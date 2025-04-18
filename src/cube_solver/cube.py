"""Cube module."""
from typing import Optional
import numpy as np

from cube_solver.constants import COLORS, FACES, AXES, OPPOSITE_FACE, MOVE_COUNT_STR, REPR_ORDER
from cube_solver.constants import FACE_MOVES, ARRAY_MOVES, CUBIE_IDX


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

        elif self.representation == "cubie":
            self.cubies = np.full((self.size,) * 4, "K")
            for face, color in zip(FACES, COLORS):
                self.cubies[CUBIE_IDX[face] + (AXES[face],)] = color

    def apply_move(self, move: str) -> None:
        base_move = move[0]
        shift = len(move)
        if shift > 1 and move[1] == "'":
            shift = -1

        if self.representation == "face":
            face_move = FACE_MOVES[base_move]
            self.faces[tuple(np.hstack(face_move))] = self.faces[tuple(np.hstack(np.roll(face_move, shift, axis=2)))]

        elif self.representation == "array":
            self.array[ARRAY_MOVES[base_move]] = self.array[np.roll(ARRAY_MOVES[base_move], shift, axis=1)]

        elif self.representation == "cubie":
            cubies = np.rot90(self.cubies[CUBIE_IDX[base_move]], -shift)
            if shift % 2 == 1:
                cubies = np.flip(cubies, axis=2)
                if base_move not in "FB":
                    cubies = np.roll(cubies, shift=1 if base_move in "UD" else -1, axis=2)
            self.cubies[CUBIE_IDX[base_move]] = cubies

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
            if base_move in "UFR":
                opts -= {OPPOSITE_FACE[base_move]}
            base_move = np.random.choice([*opts])
            scramble.append(base_move + count_str)

        return " ".join(scramble)

    def __repr__(self) -> str:
        if self.representation == "face":
            repr = "".join(self.faces[REPR_ORDER].flatten())

        elif self.representation == "array":
            repr = "".join(np.ravel([self.array[i:i+self.face_area] for i in np.multiply(REPR_ORDER, self.face_area)]))

        elif self.representation == "cubie":
            repr = "".join(np.ravel([self.cubies[CUBIE_IDX[face] + (AXES[face],)] for face in np.array([*FACES])[REPR_ORDER]]))

        return repr

    def __str__(self) -> str:
        if self.representation == "face":
            # up face
            str = "  " * self.size + "  " + "--" * self.size + "---\n"
            for i in range(self.size):
                str += "  " * self.size + "  | " + " ".join(self.faces[REPR_ORDER[0], i, :]) + " |\n"

            # lateral faces
            str += "--------" * self.size + "---------\n"
            for i in range(self.size):
                str += "| " + " | ".join([" ".join(row) for row in self.faces[REPR_ORDER[1:-1], i, :]]) + " |\n"
            str += "--------" * self.size + "---------\n"

            # down face
            for i in range(self.size):
                str += "  " * self.size + "  | " + " ".join(self.faces[REPR_ORDER[-1], i, :]) + " |\n"
            str += "  " * self.size + "  " + "--" * self.size + "---"

        elif self.representation == "array":
            # up face
            str = "  " * self.size + "  " + "--" * self.size + "---\n"
            for i in range(self.size):
                j = REPR_ORDER[0] * self.face_area + i * self.size
                str += "  " * self.size + "  | " + " ".join(self.array[j:j+3]) + " |\n"

            # lateral faces
            str += "--------" * self.size + "---------\n"
            for i in range(self.size):
                js = [order * self.face_area + i * self.size for order in REPR_ORDER[1:-1]]
                str += "| " + " | ".join([" ".join(self.array[j:j+3]) for j in js]) + " |\n"
            str += "--------" * self.size + "---------\n"

            # down face
            for i in range(self.size):
                j = REPR_ORDER[-1] * self.face_area + i * self.size
                str += "  " * self.size + "  | " + " ".join(self.array[j:j+3]) + " |\n"
            str += "  " * self.size + "  " + "--" * self.size + "---"

        elif self.representation == "cubie":
            # up face
            str = "  " * self.size + "  " + "--" * self.size + "---\n"
            for i in range(self.size):
                str += "  " * self.size + "  | " + " ".join(self.cubies[0, i, :, 0]) + " |\n"

            # lateral faces
            str += "--------" * self.size + "---------\n"
            for i in range(self.size):
                order = np.array([*FACES])[REPR_ORDER[1:-1]]
                row = [" ".join(self.cubies[(i,) + CUBIE_IDX[face][1:] + (AXES[face],)]) for face in order]
                str += "| " + " | ".join(row) + " |\n"
            str += "--------" * self.size + "---------\n"

            # down face
            for i in range(self.size - 1, -1, -1):
                str += "  " * self.size + "  | " + " ".join(self.cubies[2, i, :, 0]) + " |\n"
            str += "  " * self.size + "  " + "--" * self.size + "---"

        return str


if __name__ == "__main__":
    cube = Cube("F", representation="face")
    repr(cube)
