"""Cube module."""
from typing import Optional
import numpy as np

from cube_solver.constants import COLORS, FACES, OPPOSITE_FACE, MOVE_COUNT_STR, REPR_ORDER, FACE_MOVES, CUBIE_MOVES


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

        elif self.representation == "cubie":
            self.cubies = np.full((3, 3, 3, 3), "K")
            self.cubies[0, :, :, 0] = "W"
            self.cubies[:, 2, :, 1] = "G"
            self.cubies[:, :, 2, 2] = "R"
            self.cubies[2, :, :, 0] = "Y"
            self.cubies[:, 0, :, 1] = "B"
            self.cubies[:, :, 0, 2] = "O"

    def apply_move(self, move: str) -> None:
        base_move = move[0]
        shift = len(move)
        if shift > 1 and move[1] == "'":
            shift = -1  # same as 3

        if self.representation == "face":
            for indices in FACE_MOVES[base_move]:
                self.faces[indices] = self.faces[tuple(np.roll(indices, shift, axis=1))]

        elif self.representation == "cubie":
            for indices in CUBIE_MOVES[base_move]:
                cubies = self.cubies[tuple(np.roll(indices, shift, axis=1))]
                if shift != 2:  # if quarter turn
                    cubies = np.fliplr(cubies)
                    if base_move in "UDRL":
                        cubies = np.roll(cubies, shift=1 if base_move in "UD" else -1, axis=1)
                self.cubies[indices] = cubies

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
                    str += self.faces[REPR_ORDER[0]][j][k] + " "
                str += "| \n"

            # lateral faces
            str += "--------" * self.size + "---------\n"
            for j in range(self.size):
                str += "| "
                for i in REPR_ORDER[1:-1]:
                    for k in range(self.size):
                        str += self.faces[i][j][k] + " "
                    str += "| "
                str += "\n"
            str += "--------" * self.size + "---------\n"

            # down face
            for j in range(self.size):
                str += "  " * self.size + "  | "
                for k in range(self.size):
                    str += self.faces[REPR_ORDER[-1]][j][k] + " "
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
