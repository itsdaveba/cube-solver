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
RESTRICT_MOVES = [Move.F1, Move.F3, Move.B1, Move.B3]
PHASE1_MOVES = [move for move in PHASE0_MOVES if move not in RESTRICT_MOVES]
RESTRICT_MOVES = [Move.R1, Move.R3, Move.L1, Move.L3]
PHASE2_MOVES = [move for move in PHASE1_MOVES if move not in RESTRICT_MOVES]
RESTRICT_MOVES = [Move.U1, Move.U3, Move.D1, Move.D3]
PHASE3_MOVES = [move for move in PHASE2_MOVES if move not in RESTRICT_MOVES]

NUM_THREADS = 6  # TODO double check
THREAD_PERM_GROUP = [1, 0, 4, 5, 2, 3]
THREAD_SELECTOR = np.full((NUM_THREADS, NUM_THREADS), NONE, dtype=int)
filter = np.eye(NUM_THREADS, dtype=bool)
for i in range(NUM_THREADS):
    THREAD_SELECTOR[filter] = i
    filter = filter[THREAD_PERM_GROUP] if i % 2 == 0 else np.rot90(filter, 2)
CORNER_THREAD = np.full((NUM_THREADS, NUM_THREADS), NONE, dtype=int)
for i, y in enumerate(THREAD_SELECTOR):
    CORNER_THREAD[y, range(NUM_THREADS)] = i
CORNER_THREAD = np.vstack((CORNER_THREAD, CORNER_THREAD[THREAD_PERM_GROUP]))
CORNER_THREAD = np.hstack((CORNER_THREAD, CORNER_THREAD[:, THREAD_PERM_GROUP]))
CORNER_THREAD = np.vstack((CORNER_THREAD, np.flipud(CORNER_THREAD)))
CORNER_THREAD = np.hstack((CORNER_THREAD, np.fliplr(CORNER_THREAD)))


class Thistlethwaite(BaseSolver):
    num_phases = 4
    partial_corner_perm = True
    partial_edge_perm = True
    phase_moves = [PHASE0_MOVES, PHASE1_MOVES, PHASE2_MOVES, PHASE3_MOVES]
    pruning_kwargs = [
        [PruningDef(name="eo", shape=(EO_SIZE,), indexes=None)],
        [PruningDef(name="co_eec", shape=(CO_SIZE, EEC_SIZE), indexes=None)],
        [PruningDef(name="cct", shape=(CC_SIZE, CC_SIZE, NUM_THREADS), indexes=None)],
        [PruningDef(name="cp", shape=(24, 4, 24, 24, 12), indexes=None)]]

    def phase_coords(self, phase: int, coords: FlattenCoords) -> FlattenCoords:
        if phase == 0:
            edge_orientation = coords[1]
            return (edge_orientation,)
        elif phase == 1:
            corner_orientation = coords[0]
            equator_edge_combination = coords[5] // FACTORIAL[NUM_ORBIT_ELEMS].item()
            return (corner_orientation, equator_edge_combination)
        elif phase == 2:
            corner_combination = coords[2] // FACTORIAL[NUM_ORBIT_ELEMS].item()
            edge_combination = coords[4] // FACTORIAL[NUM_ORBIT_ELEMS].item()
            corner_thread = CORNER_THREAD[coords[2] % FACTORIAL[NUM_ORBIT_ELEMS], coords[3] % FACTORIAL[NUM_ORBIT_ELEMS]].item()
            return (corner_combination, edge_combination, corner_thread)
        elif phase == 3:
            corner_permutation = (coords[2] % 24, (coords[3] % 24) // 6)
            edge_permutation = (coords[4] % 24, coords[5] % 24, (coords[6] % 24) // 2)
            return corner_permutation + edge_permutation
        raise ValueError("")  # TODO add message
