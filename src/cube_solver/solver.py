from typing import List

from cube_solver import Cube
from cube_solver.constants import COLORS, FACES, OPPOSITE_FACE, MOVE_COUNT_STR, REPR_ORDER


class Solver:
    def __init__(self, cube: Cube = None):
        if cube is None:
            cube = Cube()
        self.cube = cube
        if cube.representation == "coord":
            self.SOLVED_COORD = (0, 0, 0, 0)  # (corner orientation, edge orientation, corner permutation, edge permutation)
        else:
            self.SOLVED_REPR = "".join([COLORS[r] * cube.face_area for r in REPR_ORDER])

    def is_solved(self, cube: Cube = None) -> bool:
        if cube.representation == "coord":
            return cube.get_coord() == self.SOLVED_COORD
        else:
            return repr(cube) == self.SOLVED_REPR

    def solve(self, max_depth: int = 5) -> str:
        solution = []
        for depth in range(max_depth + 1):
            if self._solve(depth, solution, self.cube):
                break
        return " ".join(solution[::-1])

    def _solve(self, depth: int, solution: List[str], cube: Cube, last_base_move: str = None) -> bool:
        if depth == 0:
            return self.is_solved(cube)

        base_moves = set(FACES) - {last_base_move}
        if last_base_move is not None and last_base_move in "DBL":
            base_moves -= {OPPOSITE_FACE[last_base_move]}

        for base_move in base_moves:
            for count_str in MOVE_COUNT_STR:
                next_cube = cube.apply_move(base_move + count_str, cube)
                if self._solve(depth - 1, solution, next_cube, base_move):
                    solution.append(base_move + count_str)
                    return True
        return False
