from typing import List

from cube_solver import Cube
from cube_solver.constants import COLORS, FACES, MOVE_COUNT_STR


class Solver:
    def __init__(self, cube: Cube = None):
        if cube is None:
            cube = Cube()
        self.cube = cube
        self.SOLVED_REPR = "".join([color * cube.size * cube.size for color in COLORS])

    def scramble(self, scramble):
        self.cube.reset()
        self.cube.apply_maneuver(scramble)

    def is_solved(self):
        return repr(self.cube) == self.SOLVED_REPR

    def solve(self, max_depth: int = 4) -> str:
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
                self.cube.apply_move(move)
                if self._solve(depth - 1, solution):
                    solution.append(move + MOVE_COUNT_STR[i-2])
                    return True
        return False
