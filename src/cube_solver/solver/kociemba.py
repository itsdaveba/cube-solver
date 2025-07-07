"""Kociemba solver."""
from ..cube.enums import Move
from ..cube.defs import CORNER_ORIENTATION_SIZE as CO_SIZE
from ..cube.defs import EDGE_ORIENTATION_SIZE as EO_SIZE
from ..cube.defs import CORNER_PERMUTATION_SIZE as CP_SIZE
from ..cube.defs import NUM_CORNERS, NUM_EDGES, FACTORIAL, COMBINATION, NUM_ORBIT_ELEMS
from .defs import FlattenCoords, PruningDef
from .solver import BaseSolver


CC_SIZE = COMBINATION[NUM_CORNERS, NUM_ORBIT_ELEMS].item()  # corner combination
EEC_SIZE = COMBINATION[NUM_EDGES, NUM_ORBIT_ELEMS].item()  # equator edge combination
MSEP_SIZE = FACTORIAL[NUM_EDGES - NUM_ORBIT_ELEMS].item()  # middle standing edge permutation
OP_SIZE = FACTORIAL[NUM_ORBIT_ELEMS].item()  # orbit permutation

PHASE0_MOVES = [*Move.face_moves()]
RESTRICT_MOVES = [Move.F1, Move.F3, Move.B1, Move.B3, Move.R1, Move.R3, Move.L1, Move.L3]
PHASE1_MOVES = [move for move in PHASE0_MOVES if move not in RESTRICT_MOVES]


class Kociemba(BaseSolver):
    num_phases = 2
    partial_corner_perm = False
    partial_edge_perm = True
    phase_moves = [PHASE0_MOVES, PHASE1_MOVES]
    pruning_defs = [
        [PruningDef(name="co_eo", shape=(CO_SIZE, EO_SIZE), indexes=(0, 1)),
         PruningDef(name="co_eec", shape=(CO_SIZE, EEC_SIZE), indexes=(0, 2)),
         PruningDef(name="eo_eec", shape=(EO_SIZE, EEC_SIZE), indexes=(1, 2))],
        [PruningDef(name="cp_msep", shape=(CP_SIZE, MSEP_SIZE), indexes=(0, 1)),
         PruningDef(name="cp_eep", shape=(CP_SIZE, OP_SIZE), indexes=(0, 2)),
         PruningDef(name="msep_eep", shape=(MSEP_SIZE, OP_SIZE), indexes=(1, 2))]]

    @staticmethod
    def phase_coords(coords: FlattenCoords, phase: int) -> FlattenCoords:
        if phase == 0:
            corner_orientation = coords[0]
            edge_orientation = coords[1]
            equator_edge_combination = coords[4] // OP_SIZE
            return (corner_orientation, edge_orientation, equator_edge_combination)
        elif phase == 1:
            corner_permutation = coords[2]
            middle_standing_edge_permutation = coords[5] + (coords[3] - CC_SIZE + 1 + coords[3] // OP_SIZE) * OP_SIZE
            equator_edge_permutation = coords[4] % OP_SIZE
            return (corner_permutation, middle_standing_edge_permutation, equator_edge_permutation)
        raise ValueError(f"phase must be >= 0 and < {Kociemba.num_phases} (got {phase})")
