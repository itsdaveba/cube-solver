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
from ..cube.defs import CORNER_ORIENTATION_SIZE as CO_SIZE
from ..cube.defs import EDGE_ORIENTATION_SIZE as EO_SIZE
from ..cube.defs import CORNER_PERMUTATION_SIZE as CP_SIZE
from ..cube.defs import PARTIAL_EDGE_PERMUTATION_SIZE as PEP_SIZE
from .defs import FlattenCoords, PruningDef
from .solver import BaseSolver


class Korf(BaseSolver):
    partial_corner_perm = False
    partial_edge_perm = True
    pruning_defs = [[
        PruningDef(name="co", shape=CO_SIZE, indexes=0), PruningDef(name="eo", shape=EO_SIZE, indexes=1),
        PruningDef(name="cp", shape=CP_SIZE, indexes=2), PruningDef(name="mep", shape=PEP_SIZE, indexes=3),
        PruningDef(name="eep", shape=PEP_SIZE, indexes=4), PruningDef(name="sep", shape=PEP_SIZE, indexes=5)]]

    @staticmethod
    def phase_coords(coords: FlattenCoords, phase: int) -> FlattenCoords:
        return coords
