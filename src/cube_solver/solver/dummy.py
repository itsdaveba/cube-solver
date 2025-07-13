"""
DummySolver solver.

Implementation of a simple one-phase algorithm that does not use pruning tables.

Examples
--------
>>> from cube_solver import Cube, DummySolver
>>> solver = DummySolver()
>>> cube = Cube("L2 U R D'")
>>> solver.solve(cube)
"D R' U' L2"
"""
from typing import Tuple

from .solver import BaseSolver


class DummySolver(BaseSolver):
    @staticmethod
    def phase_coords(coords: Tuple[int, int], phase: int) -> Tuple[int, ...]:
        return coords
