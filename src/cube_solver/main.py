"""Main module."""
from typing import Optional
import numpy as np

from cube_solver.constants import MOVE, OPPOSITE


class Cube:
    def __init__(self, scramble: Optional[str] = None, size: int = 3) -> None:
        assert size > 0, "size must be greater than 0"

        self.size = size
        self.faces = np.array([[[color] * size] * size for color in "WOGRBY"])
        if scramble is not None:
            self.apply_maneuver(scramble)

    def apply_move(self, move: str) -> None:
        shift = len(move)
        if shift > 1 and move[1] == "'":
            shift = -1
        for m in MOVE[move[0]]:
            self.faces[m] = self.faces[tuple(np.roll(m, shift, axis=1))]

    def apply_maneuver(self, maneuver: str) -> None:
        for move in maneuver.split():
            self.apply_move(move)

    @staticmethod
    def generate_scramble(length: int = 25) -> str:
        assert length >= 1

        options = set("ULFRBD")
        repetition = np.random.choice(3, size=length)
        repetition = list(map(lambda x: ["'", "", "2"][x], repetition))

        move = np.random.choice(list(options)) + repetition[0]
        scramble = [move]

        for rep in repetition[1:]:
            opt = options - {move[0]}
            if move[0] in "ULF":
                opt -= {OPPOSITE[move[0]]}
            move = np.random.choice(list(opt)) + rep
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


if __name__ == "__main__":
    cube = Cube("D F' U B2 R' B' F D2 U' L U B' R D' U F D R' F' U2 F' L F2 R2 D2")
    Cube.generate_scramble()
    print(cube)
