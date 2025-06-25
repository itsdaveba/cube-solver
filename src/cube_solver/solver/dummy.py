from .solver import BaseSolver, FlattenCoords


class DummySolver(BaseSolver):
    partial_corner_perm = True
    partial_edge_perm = True

    def phase_coords(self, phase: int, coords: FlattenCoords) -> FlattenCoords:
        return coords
