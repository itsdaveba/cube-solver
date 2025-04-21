"""Cube module."""
from typing import Optional
from copy import deepcopy
import numpy as np

from cube_solver.constants import COLORS, FACES, AXES, OPPOSITE_FACE, MOVE_COUNT_STR, REPR_ORDER, NUM_CORNERS, NUM_EDGES
from cube_solver.constants import FACE_MOVES, ARRAY_MOVES, SWAP, CUBIE_IDX, COORD_MOVES, COORD_CUBIE_INDEX


class Cube:
    def __init__(self, scramble: Optional[str] = None, size: int = 3, representation="face") -> None:
        assert size > 0, "size must be greater than 0"
        self.size = size
        self.face_area = size * size
        self.representation = representation

        self.reset()
        if scramble is not None:
            self.apply_maneuver(scramble)

    def reset(self) -> None:
        if self.representation == "face":
            self.faces = np.array([np.full((self.size, self.size), color) for color in COLORS])

        elif self.representation == "array":
            self.array = np.ravel([np.full((self.size, self.size), color) for color in COLORS])

        elif self.representation in ["cubie", "coord"]:
            self.cubies = np.full((self.size,) * 4, "K")
            for face, color in zip(FACES, COLORS):
                self.cubies[CUBIE_IDX[face] + (AXES[face],)] = color

            if self.representation == "coord":
                self.orientation = np.zeros(20, dtype=int)
                self.permutation = np.arange(20, dtype=int)

    def apply_move(self, move: str, cube: "Cube" = None) -> None:
        if cube is None:
            cube = self
        else:
            cube = deepcopy(cube)

        base_move = move[0]
        shift = len(move)
        if shift > 1 and move[1] == "'":
            shift = -1

        if cube.representation == "face":
            face_move = FACE_MOVES[base_move]
            cube.faces[tuple(np.hstack(face_move))] = cube.faces[tuple(np.hstack(np.roll(face_move, shift, axis=2)))]

        elif cube.representation == "array":
            cube.array[ARRAY_MOVES[base_move]] = cube.array[np.roll(ARRAY_MOVES[base_move], shift, axis=1)]

        elif cube.representation == "cubie":
            cubies = np.rot90(cube.cubies[CUBIE_IDX[base_move]], -shift)
            if shift % 2 == 1:
                cubies = cubies[..., SWAP[AXES[base_move]]]
            cube.cubies[CUBIE_IDX[base_move]] = cubies

        elif cube.representation == "coord":
            coord_move = np.roll(COORD_MOVES[base_move], shift, axis=1)
            orientation = cube.orientation[coord_move]
            if shift % 2 == 1:
                if base_move in "FB":
                    orientation = (orientation + [[1, 2, 1, 2], [1, 1, 1, 1]]) % ([3], [2])
                elif base_move in "RL":
                    orientation[0] = (orientation[0] + [2, 1, 2, 1]) % 3
            cube.orientation[COORD_MOVES[base_move]] = orientation
            cube.permutation[COORD_MOVES[base_move]] = cube.permutation[coord_move]

        return cube

    def apply_maneuver(self, maneuver: str) -> None:
        for move in maneuver.split():
            self.apply_move(move)

    @staticmethod
    def generate_scramble(length: int = 25) -> str:
        assert length >= 1

        options = set(FACES)
        count = np.random.choice(3, size=length)
        count_strs = [MOVE_COUNT_STR[c] for c in count]

        base_move = np.random.choice([*options])
        count_str = count_strs[0]
        scramble = [base_move + count_str]

        for count_str in count_strs[1:]:
            opts = options - {base_move}
            if base_move in "DBL":
                opts -= {OPPOSITE_FACE[base_move]}
            base_move = np.random.choice([*opts])
            scramble.append(base_move + count_str)

        return " ".join(scramble)

    def update_cubies(self):  # used by "coord" representation
        cubies = self.cubies.copy()
        for face, color in zip(FACES, COLORS):
            cubies[CUBIE_IDX[face] + (AXES[face],)] = color

        # corners
        for index in range(8):
            orientation = self.orientation[index]
            permutation = self.permutation[index]
            cubie = cubies[COORD_CUBIE_INDEX[permutation]]
            if (index < 4) != (permutation < 4):
                cubie = cubie[SWAP[0]]
            self.cubies[COORD_CUBIE_INDEX[index]] = np.roll(cubie, -orientation if index < 4 else orientation)

        # edges
        def axis(x): return 2 if x < 12 else 1 if x < 16 else 0
        for index in range(8, 20):
            permutation = self.permutation[index]
            cubie = cubies[COORD_CUBIE_INDEX[permutation]]
            axes = {axis(index), axis(permutation)}
            if axes == {0, 2}:
                cubie = np.roll(cubie, 1 if axis(index) != 2 else -1)
            elif axes == {1, 2}:
                cubie = cubie[SWAP[0]]
            elif axes == {0, 1}:
                cubie = cubie[SWAP[2]]
            if self.orientation[index]:
                cubie = cubie[SWAP[axis(index)]]
            self.cubies[COORD_CUBIE_INDEX[index]] = cubie

    def get_coord(self):  # corner orientation, edge orientation, corner permutation, edge permutation
        # corner orientation index
        corner_orientation_index = 0
        for co in self.orientation[:NUM_CORNERS-1]:
            corner_orientation_index *= 3
            corner_orientation_index += co

        # edge orientation index
        edge_orientation_index = 0
        for eo in self.orientation[NUM_CORNERS:-1]:
            edge_orientation_index *= 2
            edge_orientation_index += eo

        # corner permutation index
        corner_permutation_index = 0
        factorial = 1
        for i in range(NUM_CORNERS - 2, -1, -1):
            corner_permutation_index += factorial * np.sum(self.permutation[i] > self.permutation[1+i:NUM_CORNERS])
            factorial *= NUM_CORNERS - i

        # edge permutation index
        edge_permutation_index = 0
        factorial = 1
        for i in range(NUM_CORNERS + NUM_EDGES - 3, NUM_CORNERS - 1, -1):
            edge_permutation_index += factorial * np.sum(self.permutation[i] > self.permutation[1+i:])
            factorial *= NUM_EDGES - i + NUM_CORNERS

        return corner_orientation_index, edge_orientation_index, corner_permutation_index, edge_permutation_index

    def __repr__(self) -> str:
        if self.representation == "face":
            repr = "".join(self.faces[REPR_ORDER].flatten())

        elif self.representation == "array":
            repr = "".join(np.ravel([self.array[i:i+self.face_area] for i in np.multiply(REPR_ORDER, self.face_area)]))

        elif self.representation in ["cubie", "coord"]:
            if self.representation == "coord":
                self.update_cubies()
            repr = "".join(np.ravel([self.cubies[CUBIE_IDX[face] + (AXES[face],)] for face in np.array([*FACES])[REPR_ORDER]]))

        return repr

    def __str__(self) -> str:
        repr = self.__repr__()

        # up face
        str = "  " * self.size + "  " + "--" * self.size + "---\n"
        for i in range(self.size):
            j = i * self.size
            str += "  " * self.size + "  | " + " ".join(repr[j:j+self.size]) + " |\n"

        # lateral faces
        str += "--------" * self.size + "---------\n"
        for i in range(self.size):
            js = [face * self.face_area + i * self.size for face in range(1, 5)]
            str += "| " + " | ".join([" ".join(repr[j:j+self.size]) for j in js]) + " |\n"
        str += "--------" * self.size + "---------\n"

        # down face
        for i in range(self.size):
            j = 5 * self.face_area + i * self.size
            str += "  " * self.size + "  | " + " ".join(repr[j:j+self.size]) + " |\n"
        str += "  " * self.size + "  " + "--" * self.size + "---"

        return str
