"""Main module."""
import numpy as np

MOVE = {
    "U": [([0, 0, 0, 0], [0, 0, 2, 2], [0, 2, 2, 0]),
          ([1, 4, 3, 2], [0, 0, 0, 0], [0, 0, 0, 0]),
          ([1, 4, 3, 2], [0, 0, 0, 0], [2, 2, 2, 2]),
          ([0, 0, 0, 0], [0, 1, 2, 1], [1, 2, 1, 0]),
          ([1, 4, 3, 2], [0, 0, 0, 0], [1, 1, 1, 1])],
    "L": [([0, 2, 5, 4], [0, 0, 0, 2], [0, 0, 0, 2]),
          ([0, 2, 5, 4], [2, 2, 2, 0], [0, 0, 0, 2]),
          ([1, 1, 1, 1], [0, 0, 2, 2], [0, 2, 2, 0]),
          ([0, 2, 5, 4], [1, 1, 1, 1], [0, 0, 0, 2]),
          ([1, 1, 1, 1], [0, 1, 2, 1], [1, 2, 1, 0])],
    "F": [([0, 3, 5, 1], [2, 0, 0, 2], [0, 0, 2, 2]),
          ([0, 3, 5, 1], [2, 2, 0, 0], [2, 0, 0, 2]),
          ([2, 2, 2, 2], [0, 0, 2, 2], [0, 2, 2, 0]),
          ([0, 3, 5, 1], [2, 1, 0, 1], [1, 0, 1, 2]),
          ([2, 2, 2, 2], [0, 1, 2, 1], [1, 2, 1, 0])],
    "R": [([0, 4, 5, 2], [0, 2, 0, 0], [2, 0, 2, 2]),
          ([0, 4, 5, 2], [2, 0, 2, 2], [2, 0, 2, 2]),
          ([3, 3, 3, 3], [0, 0, 2, 2], [0, 2, 2, 0]),
          ([0, 4, 5, 2], [1, 1, 1, 1], [2, 0, 2, 2]),
          ([3, 3, 3, 3], [0, 1, 2, 1], [1, 2, 1, 0])],
    "B": [([0, 1, 5, 3], [0, 2, 2, 0], [0, 0, 2, 2]),
          ([0, 1, 5, 3], [0, 0, 2, 2], [2, 0, 0, 2]),
          ([4, 4, 4, 4], [0, 0, 2, 2], [0, 2, 2, 0]),
          ([0, 1, 5, 3], [0, 1, 2, 1], [1, 0, 1, 2]),
          ([4, 4, 4, 4], [0, 1, 2, 1], [1, 2, 1, 0])],
    "D": [([1, 2, 3, 4], [2, 2, 2, 2], [0, 0, 0, 0]),
          ([1, 2, 3, 4], [2, 2, 2, 2], [2, 2, 2, 2]),
          ([5, 5, 5, 5], [0, 0, 2, 2], [0, 2, 2, 0]),
          ([1, 2, 3, 4], [2, 2, 2, 2], [1, 1, 1, 1]),
          ([5, 5, 5, 5], [0, 1, 2, 1], [1, 2, 1, 0])]
}


class Cube:
    def __init__(self, scramble: str | None = None, size: int = 3) -> None:
        assert size > 0, "size must be greater than 0"
        self.size = size
        self.faces = np.array([[[color] * size for _ in range(size)] for color in "WOGRBY"])  # test * size * size
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
    print(cube)
