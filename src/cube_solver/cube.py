"""Cube module."""
from typing import Optional
import numpy as np

from cube_solver.constants import COLORS, FACES, AXES, OPPOSITE_FACE, MOVE_COUNT_STR, REPR_ORDER
from cube_solver.constants import FACE_MOVES, ARRAY_MOVES, CUBIE_INDEX


class Cube:
    def __init__(self, scramble: Optional[str] = None, size: int = 3, representation="face") -> None:
        assert size > 0, "size must be greater than 0"
        self.size = size
        self.representation = representation

        self.reset()
        if scramble is not None:
            self.apply_maneuver(scramble)

    def reset(self):
        if self.representation == "face":
            self.faces = np.array([[[color] * self.size for _ in range(self.size)] for color in COLORS])

        elif self.representation == "array":
            self.array = np.array([[color] * self.size * self.size for color in COLORS]).flatten()

        elif self.representation == "cubie":
            self.cubies = np.full((3, 3, 3, 3), "K")
            for face, color in zip(FACES, COLORS):
                self.cubies[CUBIE_INDEX[face] + (AXES[face],)] = color

    def apply_move(self, move: str) -> None:
        base_move = move[0]
        shift = len(move)
        if shift > 1 and move[1] == "'":
            shift = -1  # same as 3

        if self.representation == "face":
            new_faces = self.faces[tuple(np.hstack(np.roll(FACE_MOVES[base_move], shift, axis=2)))]
            self.faces[tuple(np.hstack(FACE_MOVES[base_move]))] = new_faces

        elif self.representation == "array":
            self.array[ARRAY_MOVES[base_move]] = self.array[np.roll(ARRAY_MOVES[base_move], shift, axis=1)]

        elif self.representation == "cubie":
            cubies = np.rot90(self.cubies[CUBIE_INDEX[base_move]], shift if base_move in "RDB" else -shift)
            if shift != 2:
                cubies = np.flip(cubies, axis=2)
                if base_move in "UDRL":
                    cubies = np.roll(cubies, shift=1 if base_move in "UD" else -1, axis=2)
            self.cubies[CUBIE_INDEX[base_move]] = cubies

    def apply_maneuver(self, maneuver: str) -> None:
        for move in maneuver.split():
            self.apply_move(move)

    @staticmethod
    def generate_scramble(length: int = 25) -> str:
        assert length >= 1

        options = set(FACES)
        count = np.random.choice(3, size=length)
        count_strs = [MOVE_COUNT_STR[c] for c in count]

        base_move = np.random.choice(list(options))
        count_str = count_strs[0]
        scramble = [base_move + count_str]

        for count_str in count_strs[1:]:
            opts = options - {base_move}
            if base_move in "UFR":
                opts -= {OPPOSITE_FACE[base_move]}
            base_move = np.random.choice(list(opts))
            scramble.append(base_move + count_str)

        return " ".join(scramble)

    def __repr__(self) -> str:
        if self.representation == "face":
            repr = "".join(self.faces[REPR_ORDER].flatten())

        elif self.representation == "array":
            repr = [self.array[face:face + self.size * self.size] for face in np.array(REPR_ORDER) * self.size * self.size]
            repr = "".join(np.array(repr).flatten())

        elif self.representation == "cubie":
            repr = ""
            # up
            for j in range(3):
                for k in range(3):
                    repr += self.cubies[0, j, k, 0]
            # left
            for i in range(3):
                for j in range(3):
                    repr += self.cubies[i, j, 0, 2]
            # front
            for i in range(3):
                for k in range(3):
                    repr += self.cubies[i, 2, k, 1]
            # right
            for i in range(3):
                for j in range(2, -1, -1):
                    repr += self.cubies[i, j, 2, 2]
            # back
            for i in range(3):
                for k in range(2, -1, -1):
                    repr += self.cubies[i, 0, k, 1]
            # down
            for j in range(2, -1, -1):
                for k in range(3):
                    repr += self.cubies[2, j, k, 0]
        return repr

    def __str__(self) -> str:
        if self.representation == "face":
            # up face
            str = "  " * self.size + "  "
            str += "--" * self.size + "---\n"
            for j in range(self.size):
                str += "  " * self.size + "  | "
                for k in range(self.size):
                    str += self.faces[REPR_ORDER[0], j, k] + " "
                str += "| \n"

            # lateral faces
            str += "--------" * self.size + "---------\n"
            for j in range(self.size):
                str += "| "
                for i in REPR_ORDER[1:-1]:
                    for k in range(self.size):
                        str += self.faces[i, j, k] + " "
                    str += "| "
                str += "\n"
            str += "--------" * self.size + "---------\n"

            # down face
            for j in range(self.size):
                str += "  " * self.size + "  | "
                for k in range(self.size):
                    str += self.faces[REPR_ORDER[-1], j, k] + " "
                str += "| \n"
            str += "  " * self.size + "  "
            str += "--" * self.size + "---"

        elif self.representation == "array":
            # up face
            str = "  " * self.size + "  "
            str += "--" * self.size + "---\n"
            for j in range(self.size):
                str += "  " * self.size + "  | "
                for k in range(self.size):
                    str += self.array[REPR_ORDER[0] * self.size * self.size + j * self.size + k] + " "
                str += "| \n"

            # lateral faces
            str += "--------" * self.size + "---------\n"
            for j in range(self.size):
                str += "| "
                for i in REPR_ORDER[1:-1]:
                    for k in range(self.size):
                        str += self.array[i * self.size * self.size + j * self.size + k] + " "
                    str += "| "
                str += "\n"
            str += "--------" * self.size + "---------\n"

            # down face
            for j in range(self.size):
                str += "  " * self.size + "  | "
                for k in range(self.size):
                    str += self.array[REPR_ORDER[-1] * self.size * self.size + j * self.size + k] + " "
                str += "| \n"
            str += "  " * self.size + "  "
            str += "--" * self.size + "---"

        elif self.representation == "cubie":
            # up face
            str = "  " * self.size + "  "
            str += "--" * self.size + "---\n"
            for j in range(self.size):
                str += "  " * self.size + "  | "
                for k in range(self.size):
                    str += self.cubies[0, j, k, 0] + " "
                str += "| \n"

            # lateral faces
            str += "--------" * self.size + "---------\n"
            for i in range(self.size):
                # left
                str += "| "
                for j in range(self.size):
                    str += self.cubies[i, j, 0, 2] + " "
                # front
                str += "| "
                for k in range(self.size):
                    str += self.cubies[i, 2, k, 1] + " "
                # right
                str += "| "
                for j in range(self.size - 1, -1, -1):
                    str += self.cubies[i, j, 2, 2] + " "
                # back
                str += "| "
                for k in range(self.size - 1, -1, -1):
                    str += self.cubies[i, 0, k, 1] + " "
                str += "| \n"
            str += "--------" * self.size + "---------\n"

            # down face
            for j in range(self.size - 1, -1, -1):
                str += "  " * self.size + "  | "
                for k in range(self.size):
                    str += self.cubies[2, j, k, 0] + " "
                str += "| \n"
            str += "  " * self.size + "  "
            str += "--" * self.size + "---"

        return str


if __name__ == "__main__":
    cube = Cube("F", representation="face")
    repr(cube)
