"""Cube module."""
from typing import Optional
import numpy as np

from cube_solver.constants import COLORS, FACES, OPPOSITE_FACE, MOVE_COUNT_STR, FACE_MOVES


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
        shift = len(move)
        if shift > 1 and move[1] == "'":
            shift = -1

        for indices in FACE_MOVES[move[0]]:
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

        move = np.random.choice(list(options)) + count_strs[0]
        scramble = [move]

        for count_str in count_strs[1:]:
            opt = options - {move[0]}
            if move[0] in FACES[:3]:
                opt -= {OPPOSITE_FACE[move[0]]}
            move = np.random.choice(list(opt)) + count_str
            scramble.append(move)

        return " ".join(scramble)

    def __repr__(self) -> str:
        return "".join(self.faces.flatten())

    def __str__(self) -> str:
        # up face
        repr = "  " * self.size + "  "
        repr += "--" * self.size + "---\n"
        for j in range(self.size):
            repr += "  " * self.size + "  | "
            for k in range(self.size):
                repr += self.faces[0][j][k] + " "
            repr += "| \n"

        # lateral faces
        repr += "--------" * self.size + "---------\n"
        for j in range(self.size):
            repr += "| "
            for i in range(1, 5):
                for k in range(self.size):
                    repr += self.faces[i][j][k] + " "
                repr += "| "
            repr += "\n"
        repr += "--------" * self.size + "---------\n"

        # down face
        for j in range(self.size):
            repr += "  " * self.size + "  | "
            for k in range(self.size):
                repr += self.faces[5][j][k] + " "
            repr += "| \n"
        repr += "  " * self.size + "  "
        repr += "--" * self.size + "---"

        return repr
