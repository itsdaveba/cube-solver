import pickle
import numpy as np
from copy import deepcopy

from cube_solver import Cube
from cube_solver.constants import ALL_MOVES, ALL_MOVES_INDEX, NEXT_MOVES
from cube_solver.constants import SOLVED_COORD, SOLVED_REPR
from cube_solver.constants import CORNER_ORIENTATION_SIZE, EDGE_ORIENTATION_SIZE
from cube_solver.constants import CORNER_PERMUTATION_SIZE, EDGE_PERMUTATION_SIZE


class Solver:
    def __init__(self, cube: Cube, transition_tables: bool = False):
        if cube is None:
            cube = Cube()
        self.cube = deepcopy(cube)
        self.transition_tables = transition_tables

        if self.transition_tables and cube.representation == "coord":
            try:
                with open("co.pkl", "rb") as file:
                    self.trans_corner_orientation = pickle.load(file)
            except FileNotFoundError:
                trans_corner_orientation = np.zeros((CORNER_ORIENTATION_SIZE, len(ALL_MOVES)), dtype=np.uint16)
                for corner_orientation in range(CORNER_ORIENTATION_SIZE):
                    cube.set_coord((corner_orientation, 0, 0, (0, 0, 0)))
                    for i in range(len(ALL_MOVES)):
                        trans_corner_orientation[corner_orientation, i] = cube.apply_move(ALL_MOVES[i], cube).get_coord()[0]
                with open("co.pkl", "wb") as file:
                    pickle.dump(trans_corner_orientation, file)
                self.trans_corner_orientation = trans_corner_orientation

            try:
                with open("eo.pkl", "rb") as file:
                    self.trans_edge_orientation = pickle.load(file)
            except FileNotFoundError:
                trans_edge_orientation = np.zeros((EDGE_ORIENTATION_SIZE, len(ALL_MOVES)), dtype=np.uint16)
                for edge_orientation in range(EDGE_ORIENTATION_SIZE):
                    cube.set_coord((0, edge_orientation, 0, (0, 0, 0)))
                    for i in range(len(ALL_MOVES)):
                        trans_edge_orientation[edge_orientation, i] = cube.apply_move(ALL_MOVES[i], cube).get_coord()[1]
                with open("eo.pkl", "wb") as file:
                    pickle.dump(trans_edge_orientation, file)
                self.trans_edge_orientation = trans_edge_orientation

            try:
                with open("cp.pkl", "rb") as file:
                    self.trans_corner_permutation = pickle.load(file)
            except FileNotFoundError:
                trans_corner_permutation = np.zeros((CORNER_PERMUTATION_SIZE, len(ALL_MOVES)), dtype=np.uint16)
                for corner_permutation in range(CORNER_PERMUTATION_SIZE):
                    cube.set_coord((0, 0, corner_permutation, (0, 0, 0)))
                    for i in range(len(ALL_MOVES)):
                        trans_corner_permutation[corner_permutation, i] = cube.apply_move(ALL_MOVES[i], cube).get_coord()[2]
                with open("cp.pkl", "wb") as file:
                    pickle.dump(trans_corner_permutation, file)
                self.trans_corner_permutation = trans_corner_permutation

            try:
                with open("ep.pkl", "rb") as file:
                    self.trans_edge_permutation = pickle.load(file)
            except FileNotFoundError:
                trans_edge_permutation = np.zeros((EDGE_PERMUTATION_SIZE, len(ALL_MOVES)), dtype=np.uint16)
                for edge_permutation in range(EDGE_PERMUTATION_SIZE):
                    cube.set_coord((0, 0, 0, (0, 0, edge_permutation)))
                    for i in range(len(ALL_MOVES)):
                        trans_edge_permutation[edge_permutation, i] = cube.apply_move(ALL_MOVES[i], cube).get_coord()[3][2]
                with open("ep.pkl", "wb") as file:
                    pickle.dump(trans_edge_permutation, file)
                self.trans_edge_permutation = trans_edge_permutation

    def get_next_position(self, position: Cube | tuple, move: str) -> Cube | tuple:
        if self.transition_tables and self.cube.representation == "coord":
            return (self.trans_corner_orientation[position[0], ALL_MOVES_INDEX[move]],
                    self.trans_edge_orientation[position[1], ALL_MOVES_INDEX[move]],
                    self.trans_corner_permutation[position[2], ALL_MOVES_INDEX[move]],
                    tuple([self.trans_edge_permutation[position[3][axis], ALL_MOVES_INDEX[move]] for axis in range(3)]))
        return position.apply_move(move, position)

    def is_solved(self, position: Cube | tuple) -> bool:
        if self.cube.representation == "coord":
            return (position if self.transition_tables else position.get_coord()) == SOLVED_COORD
        return repr(position) == SOLVED_REPR

    def solve(self, max_depth: int = 5) -> str:
        solution = []
        position = self.cube.get_coord() if self.transition_tables and self.cube.representation == "coord" else self.cube
        for depth in range(max_depth + 1):
            if self._solve(depth, solution, position):
                break
        return " ".join(solution[::-1])

    def _solve(self, depth: int, solution: list[str], position: Cube | tuple, last_move: str = None) -> bool:
        if depth == 0:
            return self.is_solved(position)
        for move in NEXT_MOVES[last_move]:
            next_position = self.get_next_position(position, move)
            if self._solve(depth - 1, solution, next_position, move):
                solution.append(move)
                return True
        return False
