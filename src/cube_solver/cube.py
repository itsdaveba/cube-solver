"""Cube module."""
from typing import Optional
import numpy as np

from cube_solver.constants import COLORS, FACES, OPPOSITE_FACE, MOVE_COUNT_STR, REPR_ORDER, FACE_MOVES


class Cube:
    def __init__(self, scramble: Optional[str] = None, size: int = 3) -> None:
        assert size > 0, "size must be greater than 0"

        self.size = size
        self.reset()

        if scramble is not None:
            self.apply_maneuver(scramble)

    def reset(self):
        self.faces = np.array([[[color] * self.size for _ in range(self.size)] for color in COLORS])

    def apply_move(self, move: str) -> None:
        base_move = move[0]
        shift = len(move)
        if shift > 1 and move[1] == "'":
            shift = -1  # same as 3

        for indices in FACE_MOVES[base_move]:
            self.faces[indices] = self.faces[tuple(np.roll(indices, shift, axis=1))]

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

        for move_count in count_strs[1:]:
            opts = options - {base_move}
            if base_move in "UFR":
                opts -= {OPPOSITE_FACE[base_move]}
            base_move = np.random.choice(list(opts))
            scramble.append(base_move + count_str)

        return " ".join(scramble)

    def __repr__(self) -> str:
        return "".join(self.faces[REPR_ORDER].flatten())

    def __str__(self) -> str:
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

        return str
