import numpy as np

from ..cube.enums import Move
from .solver import BaseSolver, PruningDef, FlattenCoords
from ..cube.defs import CORNER_ORIENTATION_SIZE as CO_SIZE
from ..cube.defs import EDGE_ORIENTATION_SIZE as EO_SIZE
from ..cube.defs import PARTIAL_CORNER_PERMUTATION_SIZE as PCP_SIZE
from ..cube.defs import PARTIAL_EDGE_PERMUTATION_SIZE as PEP_SIZE
from ..cube.defs import NONE, NUM_CORNERS, NUM_EDGES, FACTORIAL, COMBINATION, NUM_ORBIT_ELEMS


EEC_SIZE = COMBINATION[NUM_EDGES, NUM_ORBIT_ELEMS].item()
CC_SIZE = COMBINATION[NUM_CORNERS, NUM_ORBIT_ELEMS].item()
EC_SIZE = COMBINATION[NUM_EDGES - NUM_ORBIT_ELEMS, NUM_ORBIT_ELEMS].item()

PHASE0_MOVES = [*Move.face_moves()]
RESTRICT_MOVES = [Move.F1, Move.F3, Move.B1, Move.B3, Move.R1, Move.R3, Move.L1, Move.L3]
PHASE1_MOVES = [move for move in PHASE0_MOVES if move not in RESTRICT_MOVES]


class Kociemba(BaseSolver):
    num_phases = 2
    partial_corner_perm = False
    partial_edge_perm = True
    phase_moves = [PHASE0_MOVES, PHASE1_MOVES]
    pruning_kwargs = [
        [PruningDef(name="co_eec", shape=(2187, 495), indexes=(0, 2)), PruningDef(name="eo_eec", shape=(2048, 495), indexes=(1, 2))],
        [PruningDef(name="cp_eep", shape=(40320, 24), indexes=(0, 2)), PruningDef(name="ep_eep", shape=(40320, 24), indexes=(1, 2))]]

    def phase_coords(self, phase: int, coords: FlattenCoords) -> FlattenCoords:
        if phase == 0:
            edge_orientation = coords[1]
            corner_orientation = coords[0]
            equator_edge_combination = coords[4] // FACTORIAL[NUM_ORBIT_ELEMS].item()
            return (corner_orientation, edge_orientation, equator_edge_combination)
        elif phase == 1:
            corner_permutation = coords[2]
            edge_permutation = coords[3] + (coords[5] + coords[5] // 24 - 69) * 24  # TODO double-check
            equator_edge_permutation = coords[4] % 24
            return (corner_permutation, edge_permutation, equator_edge_permutation)
        raise ValueError("")  # TODO add message
