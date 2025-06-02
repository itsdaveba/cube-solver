import csolver

import time
import pickle
import numpy as np
from cube_solver.solver import BaseSolver
from cube_solver.cube import Move


PHASE1_MOVES = list(Move)
RESTRICT = [Move.F1, Move.F3, Move.B1, Move.B3, Move.R1, Move.R3, Move.L1, Move.L3]
PHASE2_MOVES = [move for move in list(Move) if move not in RESTRICT]


class Kociemba(BaseSolver):
    partial_corner_perm = False
    phase_moves = phase_moves = [PHASE1_MOVES, PHASE2_MOVES]

    def _phase_coords(self, phase: int, coords: tuple) -> tuple:
        return coords


solver = Kociemba(use_transition_tables=True)

start = time.time()
table: np.ndarray = csolver.generate_pruning_table(solver, 1, (40320, 40320, 24), [0, 1, 2])
print(f"Time: {time.time() - start}")
# print(table.shape)
# print(np.unique(table, return_counts=True))
with open("tables/kociemba_cpepsep.npy", "wb") as file:
    np.save(file, table)

# print(table.min())
# print(np.unique(array, return_counts=True))
# print(np.all(table == array))
