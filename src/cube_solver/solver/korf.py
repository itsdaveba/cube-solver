from .solver import BaseSolver, PruningDef, FlattenCoords
from ..cube.defs import CORNER_ORIENTATION_SIZE as co_size
from ..cube.defs import EDGE_ORIENTATION_SIZE as eo_size
from ..cube.defs import CORNER_PERMUTATION_SIZE as cp_size
from ..cube.defs import PARTIAL_EDGE_PERMUTATION_SIZE as pep_size


class Korf(BaseSolver):
    partial_corner_perm = False
    partial_edge_perm = True
    pruning_kwargs = [[PruningDef(name="ceo", shape=(co_size, eo_size), indexes=(0, 1)),
                       PruningDef(name="cp", shape=(cp_size,), indexes=(2,)),
                       PruningDef(name="pep1", shape=(pep_size,), indexes=(3,)),
                       PruningDef(name="pep2", shape=(pep_size,), indexes=(4,)),
                       PruningDef(name="pep3", shape=(pep_size,), indexes=(5,))]]

    def phase_coords(self, phase: int, coords: FlattenCoords) -> FlattenCoords:
        return coords
