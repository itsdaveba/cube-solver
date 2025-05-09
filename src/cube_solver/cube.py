"""Cube module."""
from copy import deepcopy
import numpy as np

from cube_solver.constants import SIZE, COLORS, FACES, REPR_ORDER, NEXT_MOVES
from cube_solver.constants import FACE_MOVES, ARRAY_MOVES
from cube_solver.constants import AXES, SWAP, CUBIE_IDX
from cube_solver.constants import NUM_CORNERS, NUM_EDGES, EMPTY, NUM_AXIS_ELEMS, NUM_PERM_AXIS
from cube_solver.constants import CORNER_AXIS_OFFSET, EDGE_AXIS_OFFSET, AXIS, COMB_AXIS, COORD_MOVES, COORD_CUBIE_INDEX


class Cube:
    def __init__(self, scramble: str = None, representation="coord"):
        assert representation in ["face", "array", "cubie", "coord"], "unknown representation"
        self.representation = representation

        self.reset()
        if scramble is not None:
            self.apply_maneuver(scramble)

    def reset(self) -> None:
        if self.representation == "face":
            self.faces = np.array([np.full((SIZE, SIZE), color) for color in COLORS])

        elif self.representation == "array":
            self.array = np.ravel([np.full((SIZE, SIZE), color) for color in COLORS])

        elif self.representation in ["cubie", "coord"]:
            self.cubies = np.full((SIZE,) * 4, "K")
            for face, color in zip(FACES, COLORS):
                self.cubies[CUBIE_IDX[face] + (AXES[face],)] = color

            if self.representation == "coord":
                self.orientation = np.zeros(NUM_CORNERS + NUM_EDGES, dtype=int)
                self.permutation = np.arange(NUM_CORNERS + NUM_EDGES, dtype=int)

    def apply_move(self, move: str, cube: "Cube" = None) -> "Cube":
        cube = self if cube is None else deepcopy(cube)

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
            # orientation follows the order of axes: 0 -> UD, 1 -> FB, 2 -> RL
            coord_move = np.roll(COORD_MOVES[base_move], shift, axis=1)
            orientation = cube.orientation[coord_move]
            if shift % 2 == 1:
                if base_move in "FB":
                    orientation = (orientation + [[1, 2, 1, 2], [1, 1, 1, 1]]) % ([3], [2])  # corners and edges
                elif base_move in "RL":
                    orientation[0] = (orientation[0] + [2, 1, 2, 1]) % 3  # corners
            cube.orientation[COORD_MOVES[base_move]] = orientation
            cube.permutation[COORD_MOVES[base_move]] = cube.permutation[coord_move]

        return cube

    def apply_maneuver(self, maneuver: str) -> None:
        for move in maneuver.split():
            self.apply_move(move)

    @staticmethod
    def generate_scramble(length: int = 25) -> str:
        assert length >= 1, "scramble length must be greater or equal than 1"

        scramble = []
        move = None
        for _ in range(length):
            move = np.random.choice(NEXT_MOVES[move])
            scramble.append(move)

        return " ".join(scramble)

    def _update_cubies(self) -> None:  # used by "coord" representation
        cubies = self.cubies.copy()
        for face, color in zip(FACES, COLORS):
            cubies[CUBIE_IDX[face] + (AXES[face],)] = color

        # corners
        for index in range(NUM_CORNERS):
            orientation = self.orientation[index]
            permutation = self.permutation[index]
            cubie = cubies[COORD_CUBIE_INDEX[permutation]]
            if AXIS[index] != AXIS[permutation]:
                cubie = cubie[SWAP[0]]
            if orientation:
                cubie = np.roll(cubie, orientation if AXIS[index] else -orientation)
            self.cubies[COORD_CUBIE_INDEX[index]] = cubie

        # edges
        for index in range(NUM_CORNERS, NUM_CORNERS + NUM_EDGES):
            permutation = self.permutation[index]
            cubie = cubies[COORD_CUBIE_INDEX[permutation]]
            axes = {AXIS[index], AXIS[permutation]}
            if axes == {0, 2}:
                cubie = np.roll(cubie, 1 if AXIS[index] != 2 else -1)
            elif axes == {1, 2}:
                cubie = cubie[SWAP[0]]
            elif axes == {0, 1}:
                cubie = cubie[SWAP[2]]
            if self.orientation[index]:
                cubie = cubie[SWAP[AXIS[index]]]
            self.cubies[COORD_CUBIE_INDEX[index]] = cubie

    def get_coords(self, partial_corner: bool = False, partial_edge: bool = False) -> tuple:
        # corner orientation
        corner_orientation = 0
        for co in self.orientation[:NUM_CORNERS-1]:
            corner_orientation *= 3
            corner_orientation += co

        # edge orientation
        edge_orientation = 0
        for eo in self.orientation[NUM_CORNERS:-1]:
            edge_orientation *= 2
            edge_orientation += eo

        # corner permutation
        if not partial_corner:
            corner_permutation = 0
            for i in range(NUM_CORNERS - 1):
                corner_permutation *= NUM_CORNERS - i
                corner_permutation += np.sum(self.permutation[i] > self.permutation[i+1:NUM_CORNERS])
        else:
            corner_permutation = [0, 0]
            for cycle in range(2):
                cycle_pos = np.where(np.array([AXIS[index] for index in self.permutation[:NUM_CORNERS]]) == cycle)[0]
                cycle_permutation = self.permutation[cycle_pos]
                if len(cycle_permutation) == NUM_AXIS_ELEMS:
                    for i in range(NUM_AXIS_ELEMS - 1):
                        corner_permutation[cycle] *= NUM_AXIS_ELEMS - i
                        corner_permutation[cycle] += np.sum(cycle_permutation[i] > cycle_permutation[i+1:])
                    corner_permutation[cycle] += NUM_PERM_AXIS * np.sum(COMB_AXIS[np.arange(NUM_AXIS_ELEMS), cycle_pos])
                else:
                    corner_permutation[cycle] = EMPTY
            corner_permutation = tuple(corner_permutation)

        # edge permutation
        if not partial_edge:
            edge_permutation = 0
            for i in range(NUM_CORNERS, NUM_CORNERS + NUM_EDGES - 2):
                edge_permutation *= NUM_CORNERS + NUM_EDGES - i
                edge_permutation += np.sum(self.permutation[i] > self.permutation[i+1:])
        else:
            edge_permutation = [0, 0, 0]
            for axis in range(3):
                axis_pos = np.where(np.array([AXIS[index] for index in self.permutation[NUM_CORNERS:]]) == axis)[0]
                axis_permutation = self.permutation[axis_pos + NUM_CORNERS]
                if len(axis_permutation) == NUM_AXIS_ELEMS:
                    for i in range(NUM_AXIS_ELEMS - 1):
                        edge_permutation[axis] *= NUM_AXIS_ELEMS - i
                        edge_permutation[axis] += np.sum(axis_permutation[i] > axis_permutation[i+1:])
                    edge_permutation[axis] += NUM_PERM_AXIS * np.sum(COMB_AXIS[np.arange(NUM_AXIS_ELEMS), axis_pos])
                else:
                    edge_permutation[axis] = EMPTY
            edge_permutation = tuple(edge_permutation)

        return corner_orientation, edge_orientation, corner_permutation, edge_permutation

    def set_coords(self, coord: tuple, partial_corner: bool = False, partial_edge: bool = False) -> None:
        # corner orientation
        corner_orientation = coord[0]
        for i in range(NUM_CORNERS - 2, -1, -1):
            corner_orientation, self.orientation[i] = divmod(corner_orientation, 3)
        self.orientation[NUM_CORNERS-1] = -np.sum(self.orientation[:NUM_CORNERS-1]) % 3

        # edge orientation
        edge_orientation = coord[1]
        for i in range(NUM_CORNERS + NUM_EDGES - 2, NUM_CORNERS - 1, -1):
            edge_orientation, self.orientation[i] = divmod(edge_orientation, 2)
        self.orientation[NUM_CORNERS+NUM_EDGES-1] = -np.sum(self.orientation[NUM_CORNERS:-1]) % 2

        # corner permutation
        if not partial_corner:
            corner_parity = 0
            corner_permutation = coord[2]
            self.permutation[NUM_CORNERS-1] = 0
            for i in range(NUM_CORNERS - 2, -1, -1):
                corner_permutation, self.permutation[i] = divmod(corner_permutation, NUM_CORNERS - i)
                self.permutation[i+1:NUM_CORNERS] += self.permutation[i+1:NUM_CORNERS] >= self.permutation[i]
                corner_parity += self.permutation[i]
        else:
            corner_parity = 0  # TEST partial corner with full edge
            corner_permutation = coord[2]
            self.permutation[:NUM_CORNERS] = EMPTY
            for cycle in range(2):
                if corner_permutation[cycle] != EMPTY:
                    combination_index, permutation_index = divmod(corner_permutation[cycle], NUM_PERM_AXIS)
                    cycle_pos = np.zeros(NUM_AXIS_ELEMS, dtype=int)
                    i, pos = NUM_AXIS_ELEMS - 1, NUM_CORNERS - 1
                    while i >= 0:
                        if combination_index >= COMB_AXIS[i, pos]:
                            combination_index -= COMB_AXIS[i, pos]
                            cycle_pos[i] = pos
                            i -= 1
                        pos -= 1
                    cycle_permutation = np.zeros(NUM_AXIS_ELEMS, dtype=int)
                    for i in range(NUM_AXIS_ELEMS - 2, - 1, -1):
                        permutation_index, cycle_permutation[i] = divmod(permutation_index, NUM_AXIS_ELEMS - i)
                        cycle_permutation[i+1:] += cycle_permutation[i+1:] >= cycle_permutation[i]
                        corner_parity += cycle_permutation[i]
                    self.permutation[cycle_pos] = cycle_permutation + CORNER_AXIS_OFFSET[cycle]

        # edge permutation
        if not partial_edge:
            edge_parity = 0
            edge_permutation = coord[3]
            self.permutation[-2:] = [0, 1]
            for i in range(NUM_CORNERS + NUM_EDGES - 3, NUM_CORNERS - 1, -1):
                edge_permutation, self.permutation[i] = divmod(edge_permutation, NUM_CORNERS + NUM_EDGES - i)
                self.permutation[i+1:] += self.permutation[i+1:] >= self.permutation[i]
                edge_parity += self.permutation[i]
            self.permutation[NUM_CORNERS:] += NUM_CORNERS
            if corner_parity % 2 != edge_parity % 2:
                self.permutation[-2:] = self.permutation[[-1, -2]]
        else:
            edge_permutation = coord[3]
            self.permutation[NUM_CORNERS:] = EMPTY
            for axis in range(3):
                if edge_permutation[axis] != EMPTY:
                    combination_index, permutation_index = divmod(edge_permutation[axis], NUM_PERM_AXIS)
                    axis_pos = np.zeros(NUM_AXIS_ELEMS, dtype=int)
                    i, pos = NUM_AXIS_ELEMS - 1, NUM_EDGES - 1
                    while i >= 0:
                        if combination_index >= COMB_AXIS[i, pos]:
                            combination_index -= COMB_AXIS[i, pos]
                            axis_pos[i] = pos
                            i -= 1
                        pos -= 1
                    axis_permutation = np.zeros(NUM_AXIS_ELEMS, dtype=int)
                    for i in range(NUM_AXIS_ELEMS - 2, - 1, -1):
                        permutation_index, axis_permutation[i] = divmod(permutation_index, NUM_AXIS_ELEMS - i)
                        axis_permutation[i+1:] += axis_permutation[i+1:] >= axis_permutation[i]
                    self.permutation[axis_pos + NUM_CORNERS] = axis_permutation + EDGE_AXIS_OFFSET[axis]

    def __repr__(self) -> str:
        if self.representation == "face":
            repr = "".join(self.faces[REPR_ORDER].flatten())

        elif self.representation == "array":
            repr = "".join(np.ravel([self.array[i:i+SIZE*SIZE] for i in np.multiply(REPR_ORDER, SIZE*SIZE)]))

        elif self.representation in ["cubie", "coord"]:
            if self.representation == "coord":
                self._update_cubies()
            repr = "".join(np.ravel([self.cubies[CUBIE_IDX[face] + (AXES[face],)] for face in np.array([*FACES])[REPR_ORDER]]))

        return repr

    def __str__(self) -> str:
        repr = self.__repr__()

        # up face
        str = "  " * SIZE + "  " + "--" * SIZE + "---\n"
        for i in range(SIZE):
            j = i * SIZE
            str += "  " * SIZE + "  | " + " ".join(repr[j:j+SIZE]) + " |\n"

        # lateral faces
        str += "--------" * SIZE + "---------\n"
        for i in range(SIZE):
            js = [face * SIZE * SIZE + i * SIZE for face in range(1, 5)]
            str += "| " + " | ".join([" ".join(repr[j:j+SIZE]) for j in js]) + " |\n"
        str += "--------" * SIZE + "---------\n"

        # down face
        for i in range(SIZE):
            j = 5 * SIZE * SIZE + i * SIZE
            str += "  " * SIZE + "  | " + " ".join(repr[j:j+SIZE]) + " |\n"
        str += "  " * SIZE + "  " + "--" * SIZE + "---"

        return str
