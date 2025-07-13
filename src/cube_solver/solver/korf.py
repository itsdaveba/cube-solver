"""
Korf solver.

Implementation of a simple one-phase algorithm that uses pruning tables proposed by Richard E. Korf.

For more information, see: https://en.wikipedia.org/wiki/Optimal_solutions_for_the_Rubik%27s_Cube#Korf's_algorithm

Examples
--------
>>> from cube_solver import Cube, Korf
>>> solver = Korf()
>>> cube = Cube("L2 U R D' B2 D2 F B")
>>> solver.solve(cube)
"F' B' D2 B2 D R' U' L2"
"""
from typing import Tuple

from ..cube.defs import CORNER_ORIENTATION_SIZE as CO_SIZE
from ..cube.defs import CORNER_PERMUTATION_SIZE as CP_SIZE
from .defs import PruningDef
from .solver import BaseSolver


class Korf(BaseSolver):
    pruning_defs = [[PruningDef(name="co_cp", shape=(CO_SIZE, CP_SIZE))]]

    @staticmethod
    def phase_coords(coords: Tuple[int, int], phase: int) -> Tuple[int, ...]:
        return coords
