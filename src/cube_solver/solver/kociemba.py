"""
Kociemba solver.

Implementation of the two-phase algorithm proposed by Herbert Kociemba.

For more information, see: https://kociemba.org/cube.htm

Examples
--------
>>> from cube_solver import Cube, Kociemba
>>> solver = Kociemba()
>>> cube = Cube("L2 U R D' B2 D2 F L D")
>>> solver.solve(cube)
"F' R' U' R U F2 U' F2 R2 U' R2 U"

Solution divided by phases.

>>> solver.solve(cube, verbose=2)
["F' R' U' R", "U F2 U' F2 R2 U' R2 U"]

Find the optimal solution.

>>> solver.solve(cube, optimal=True)
'F2 R F R2 F U'
"""
from typing import Tuple

from ..cube.enums import Move
from ..cube.defs import CORNER_ORIENTATION_SIZE as CO_SIZE
from ..cube.defs import CORNER_PERMUTATION_SIZE as CP_SIZE
from .defs import MAIN_MOVES, PruningDef
from .solver import BaseSolver


PHASE0_MOVES = MAIN_MOVES
RESTRICT_MOVES = [Move.F1, Move.F3, Move.R1, Move.R3]
PHASE1_MOVES = [move for move in PHASE0_MOVES if move not in RESTRICT_MOVES]


class Kociemba(BaseSolver):
    num_phases = 2
    phase_moves = [PHASE0_MOVES, PHASE1_MOVES]
    pruning_defs = [[PruningDef(name="co", shape=CO_SIZE)], [PruningDef(name="cp", shape=CP_SIZE)]]

    @staticmethod
    def phase_coords(coords: Tuple[int, int], phase: int) -> Tuple[int, ...]:
        if phase == 0:
            corner_orientation = coords[0]
            return (corner_orientation,)
        elif phase == 1:
            corner_permutation = coords[1]
            return (corner_permutation,)
        raise ValueError(f"phase must be >= 0 and < {Kociemba.num_phases} (got {phase})")
