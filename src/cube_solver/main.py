"""Main module."""
from typing import Optional, List
import numpy as np

from cube_solver.constants import COLORS, FACES, REPS, OPPOSITE, MOVE


class Cube:
    def __init__(self, scramble: Optional[str] = None, size: int = 3) -> None:
        assert size > 0, "size must be greater than 0"

        self.size = size
        self.faces = np.array([[[color] * size] * size for color in COLORS])
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

    def is_solved(self):
        return repr(self) == "".join([color * self.size * self.size for color in COLORS])

    def solve(self, max_depth: int = 3) -> str:
        solution = []
        for depth in range(max_depth + 1):
            if self._solve(depth, solution):
                break
        return " ".join(solution[::-1])

    def _solve(self, depth: int, solution: List[str]) -> bool:
        if depth == 0:
            return self.is_solved()
        for move in FACES:
            for i in range(4):
                self.apply_move(move)
                if self._solve(depth - 1, solution):
                    solution.append(move + REPS[i-2])
                    return True
        return False

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
    depth = 3
    scramble = Cube.generate_scramble(length=depth)
    print("Scramble:", scramble)
    cube = Cube(scramble)
    solution = cube.solve()
    if solution or len(scramble) == 0:
        print("Solution:", solution)
    else:
        print("No solution found")
