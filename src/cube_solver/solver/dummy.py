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
from .defs import FlattenCoords
from .solver import BaseSolver


class DummySolver(BaseSolver):
    partial_corner_perm = True
    partial_edge_perm = True

    @staticmethod
    def phase_coords(coords: FlattenCoords, phase: int) -> FlattenCoords:
        return coords
