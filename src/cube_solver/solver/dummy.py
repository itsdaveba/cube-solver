"""
DummySolver solver.

Simple solver that does not use pruning tables.
"""
from .solver import BaseSolver, FlattenCoords


class DummySolver(BaseSolver):
    partial_corner_perm = True
    partial_edge_perm = True

    @staticmethod
    def phase_coords(coords: FlattenCoords, phase: int) -> FlattenCoords:
        return coords
