"""Korf solver."""
from ..cube.defs import CORNER_ORIENTATION_SIZE as CO_SIZE
from ..cube.defs import EDGE_ORIENTATION_SIZE as EO_SIZE
from ..cube.defs import CORNER_PERMUTATION_SIZE as CP_SIZE
from ..cube.defs import PARTIAL_EDGE_PERMUTATION_SIZE as PEP_SIZE
from .solver import BaseSolver, PruningDef, FlattenCoords


class Korf(BaseSolver):
    partial_corner_perm = False
    partial_edge_perm = True
    pruning_kwargs = [[PruningDef(name="co_eo", shape=(CO_SIZE, EO_SIZE), indexes=(0, 1)),
                       PruningDef(name="co_cp", shape=(CO_SIZE, CP_SIZE), indexes=(0, 2)),
                       PruningDef(name="co_mep", shape=(CO_SIZE, PEP_SIZE), indexes=(0, 3)),
                       PruningDef(name="co_eep", shape=(CO_SIZE, PEP_SIZE), indexes=(0, 4)),
                       PruningDef(name="co_sep", shape=(CO_SIZE, PEP_SIZE), indexes=(0, 5)),
                       PruningDef(name="eo_cp", shape=(EO_SIZE, CP_SIZE), indexes=(1, 2)),
                       PruningDef(name="eo_mep", shape=(EO_SIZE, PEP_SIZE), indexes=(1, 3)),
                       PruningDef(name="eo_eep", shape=(EO_SIZE, PEP_SIZE), indexes=(1, 4)),
                       PruningDef(name="eo_sep", shape=(EO_SIZE, PEP_SIZE), indexes=(1, 5)),
                       PruningDef(name="cp_mep", shape=(CP_SIZE, PEP_SIZE), indexes=(2, 3)),
                       PruningDef(name="cp_eep", shape=(CP_SIZE, PEP_SIZE), indexes=(2, 4)),
                       PruningDef(name="cp_sep", shape=(CP_SIZE, PEP_SIZE), indexes=(2, 5)),
                       PruningDef(name="mep_eep", shape=(PEP_SIZE, PEP_SIZE), indexes=(3, 4)),
                       PruningDef(name="mep_sep", shape=(PEP_SIZE, PEP_SIZE), indexes=(3, 5)),
                       PruningDef(name="eep_sep", shape=(PEP_SIZE, PEP_SIZE), indexes=(4, 5))]]

    def phase_coords(self, coords: FlattenCoords, phase: int) -> FlattenCoords:
        return coords
