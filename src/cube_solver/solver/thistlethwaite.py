"""
Thistlethwaite solver.

Implementation of the four-phase algorithm proposed by Morwen B. Thistlethwaite.

For more information, see: https://www.jaapsch.net/puzzles/thistle.htm

Examples
--------
>>> from cube_solver import Cube, Thistlethwaite
>>> solver = Thistlethwaite()
>>> cube = Cube("L2 U R D' B2 D2 F B D")
>>> solver.solve(cube)
"D F B U2 B2 L2 U L' D F2 B2 R2 F2 L2 D2 F2 B2 L2 F2"

Solution divided by phases.

>>> solver.solve(cube, verbose=2)
['D F B', 'U2 B2 L2 U L', 'L2 D', 'F2 B2 R2 F2 L2 D2 F2 B2 L2 F2']
"""
import numpy as np

from ..cube.enums import Move
from ..cube.defs import CORNER_ORIENTATION_SIZE as CO_SIZE
from ..cube.defs import EDGE_ORIENTATION_SIZE as EO_SIZE
from ..cube.defs import NONE, NUM_CORNERS, NUM_EDGES, FACTORIAL, COMBINATION, NUM_ORBIT_ELEMS
from .defs import FlattenCoords, PruningDef
from .solver import BaseSolver


EEC_SIZE = COMBINATION[NUM_EDGES, NUM_ORBIT_ELEMS].item()  # equator edge combination
CC_SIZE = COMBINATION[NUM_CORNERS, NUM_ORBIT_ELEMS].item()  # corner combination
MSEC_SIZE = COMBINATION[NUM_EDGES - NUM_ORBIT_ELEMS, NUM_ORBIT_ELEMS].item()  # middle standing edge combination
OP_SIZE = FACTORIAL[NUM_ORBIT_ELEMS].item()  # orbit permutation

PHASE0_MOVES = [*Move.face_moves()]
RESTRICT_MOVES = [Move.F1, Move.F3, Move.B1, Move.B3]
PHASE1_MOVES = [move for move in PHASE0_MOVES if move not in RESTRICT_MOVES]
RESTRICT_MOVES = [Move.R1, Move.R3, Move.L1, Move.L3]
PHASE2_MOVES = [move for move in PHASE1_MOVES if move not in RESTRICT_MOVES]
RESTRICT_MOVES = [Move.U1, Move.U3, Move.D1, Move.D3]
PHASE3_MOVES = [move for move in PHASE2_MOVES if move not in RESTRICT_MOVES]

NUM_THREADS = 6
thread_perm_group = [1, 0, 4, 5, 2, 3]
thread_selector = np.full((NUM_THREADS, NUM_THREADS), NONE, dtype=int)
filter = np.eye(NUM_THREADS, dtype=bool)
for i in range(NUM_THREADS):
    thread_selector[filter] = i
    filter = filter[thread_perm_group] if i % 2 == 0 else np.rot90(filter, 2)
CORNER_THREAD = np.full((NUM_THREADS, NUM_THREADS), NONE, dtype=int)
for i, y in enumerate(thread_selector):
    CORNER_THREAD[y, range(NUM_THREADS)] = i
CORNER_THREAD = np.vstack((CORNER_THREAD, CORNER_THREAD[thread_perm_group]))
CORNER_THREAD = np.hstack((CORNER_THREAD, CORNER_THREAD[:, thread_perm_group]))
CORNER_THREAD = np.vstack((CORNER_THREAD, np.flipud(CORNER_THREAD)))
CORNER_THREAD = np.hstack((CORNER_THREAD, np.fliplr(CORNER_THREAD)))


class Thistlethwaite(BaseSolver):
    num_phases = 4
    partial_corner_perm = True
    partial_edge_perm = True
    phase_moves = [PHASE0_MOVES, PHASE1_MOVES, PHASE2_MOVES, PHASE3_MOVES]
    pruning_defs = [
        [PruningDef(name="eo", shape=EO_SIZE)], [PruningDef(name="co_eec", shape=(CO_SIZE, EEC_SIZE))],
        [PruningDef(name="cc_ct_msec", shape=(CC_SIZE, NUM_THREADS, MSEC_SIZE))],
        [PruningDef(name="cop_eop", shape=(OP_SIZE, OP_SIZE // NUM_THREADS, OP_SIZE, OP_SIZE, OP_SIZE // 2))]]

    @staticmethod
    def phase_coords(coords: FlattenCoords, phase: int) -> FlattenCoords:
        if phase == 0:
            edge_orientation = coords[1]
            return (edge_orientation,)
        if phase == 1:
            corner_orientation = coords[0]
            equator_edge_combination = coords[5] // OP_SIZE
            return (corner_orientation, equator_edge_combination)
        if phase == 2:
            corner_combination = coords[2] // OP_SIZE
            corner_thread = CORNER_THREAD[coords[2] % OP_SIZE, coords[3] % OP_SIZE].item()
            middle_standing_edge_combination = coords[4] // OP_SIZE
            return (corner_combination, corner_thread, middle_standing_edge_combination)
        if phase == 3:
            corner_orbit_permutation = (coords[2], (coords[3] % OP_SIZE) // NUM_THREADS)
            edge_orbit_permutation = (coords[4], coords[5] % OP_SIZE, (coords[6] % OP_SIZE) // 2)
            return corner_orbit_permutation + edge_orbit_permutation
        raise ValueError(f"phase must be >= 0 and < {Thistlethwaite.num_phases} (got {phase})")
