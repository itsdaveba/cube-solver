"""DummySolver solver."""
from .solver import BaseSolver, FlattenCoords


class DummySolver(BaseSolver):
    partial_corner_perm = True
    partial_edge_perm = True

    def phase_coords(self, coords: FlattenCoords, phase: int) -> FlattenCoords:
        return coords
