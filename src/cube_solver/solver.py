import pickle
import numpy as np
from copy import deepcopy
from collections import deque

from cube_solver import Cube
from cube_solver.constants import ALL_MOVES, ALL_MOVES_INDEX, NEXT_MOVES
from cube_solver.constants import SOLVED_PARTIAL_COORD, SOLVED_REPR, EMPTY
from cube_solver.constants import CORNER_ORIENTATION_SIZE, EDGE_ORIENTATION_SIZE
from cube_solver.constants import CORNER_PERMUTATION_SIZE, PARTIAL_EDGE_PERMUTATION_SIZE


class Solver:
    def __init__(self, cube: Cube, pruning_tables: bool = False, transition_tables: bool = False):
        self.cube = deepcopy(cube)
        self.pruning_tables = pruning_tables
        self.transition_tables = transition_tables

        if self.pruning_tables:
            try:
                with open("pco.pkl", "rb") as file:
                    self.prun_corner_orientation = pickle.load(file)
            except FileNotFoundError:
                cube.reset()
                corner_orientation = cube.get_coord()[0]
                prun_corner_orientation = np.full(CORNER_ORIENTATION_SIZE, EMPTY, dtype=np.int8)
                prun_corner_orientation[corner_orientation] = 0
                queue = deque([(corner_orientation, 0)])  # index, depth
                while queue:
                    corner_orientation, depth = queue.popleft()
                    cube.set_coord((corner_orientation, 0, 0, (0, 0, 0)))
                    for move in ALL_MOVES:
                        corner_orientation = cube.apply_move(move, cube).get_coord()[0]
                        if prun_corner_orientation[corner_orientation] == EMPTY:
                            prun_corner_orientation[corner_orientation] = depth + 1
                            queue.append((corner_orientation, depth + 1))
                with open("pco.pkl", "wb") as file:
                    pickle.dump(prun_corner_orientation, file)
                self.prun_corner_orientation = prun_corner_orientation

            try:
                with open("peo.pkl", "rb") as file:
                    self.prun_edge_orientation = pickle.load(file)
            except FileNotFoundError:
                cube.reset()
                edge_orientation = cube.get_coord()[1]
                prun_edge_orientation = np.full(EDGE_ORIENTATION_SIZE, EMPTY, dtype=np.int8)
                prun_edge_orientation[edge_orientation] = 0
                queue = deque([(edge_orientation, 0)])  # index, depth
                while queue:
                    edge_orientation, depth = queue.popleft()
                    cube.set_coord((0, edge_orientation, 0, (0, 0, 0)))
                    for move in ALL_MOVES:
                        edge_orientation = cube.apply_move(move, cube).get_coord()[1]
                        if prun_edge_orientation[edge_orientation] == EMPTY:
                            prun_edge_orientation[edge_orientation] = depth + 1
                            queue.append((edge_orientation, depth + 1))
                with open("peo.pkl", "wb") as file:
                    pickle.dump(prun_edge_orientation, file)
                self.prun_edge_orientation = prun_edge_orientation

            try:
                with open("pcp.pkl", "rb") as file:
                    self.prun_corner_permutation = pickle.load(file)
            except FileNotFoundError:
                cube.reset()
                corner_permutation = cube.get_coord()[2]
                prun_corner_permutation = np.full(CORNER_PERMUTATION_SIZE, EMPTY, dtype=np.int8)
                prun_corner_permutation[corner_permutation] = 0
                queue = deque([(corner_permutation, 0)])  # index, depth
                while queue:
                    corner_permutation, depth = queue.popleft()
                    cube.set_coord((0, 0, corner_permutation, (0, 0, 0)))
                    for move in ALL_MOVES:
                        corner_permutation = cube.apply_move(move, cube).get_coord()[2]
                        if prun_corner_permutation[corner_permutation] == EMPTY:
                            prun_corner_permutation[corner_permutation] = depth + 1
                            queue.append((corner_permutation, depth + 1))
                with open("pcp.pkl", "wb") as file:
                    pickle.dump(prun_corner_permutation, file)
                self.prun_corner_permutation = prun_corner_permutation

            try:
                with open("pep.pkl", "rb") as file:
                    self.prun_edge_permutation = pickle.load(file)
            except FileNotFoundError:
                prun_edge_permutation = np.full((3, PARTIAL_EDGE_PERMUTATION_SIZE), EMPTY, dtype=np.int8)
                for i in range(3):
                    cube.reset()
                    edge_permutation = cube.get_coord()[3][i]
                    prun_edge_permutation[i, edge_permutation] = 0
                    queue = deque([(edge_permutation, 0)])  # index, depth
                    while queue:
                        edge_permutation = [-1, -1, -1]
                        edge_permutation[i], depth = queue.popleft()
                        cube.set_coord((0, 0, 0, edge_permutation))
                        for move in ALL_MOVES:
                            edge_permutation = cube.apply_move(move, cube).get_coord()[3][i]
                            if prun_edge_permutation[i, edge_permutation] == EMPTY:
                                prun_edge_permutation[i, edge_permutation] = depth + 1
                                queue.append((edge_permutation, depth + 1))
                with open("pep.pkl", "wb") as file:
                    pickle.dump(prun_edge_permutation, file)
                self.prun_edge_permutation = prun_edge_permutation

        if self.transition_tables:
            try:
                with open("tco.pkl", "rb") as file:
                    self.trans_corner_orientation = pickle.load(file)
            except FileNotFoundError:
                trans_corner_orientation = np.zeros((CORNER_ORIENTATION_SIZE, len(ALL_MOVES)), dtype=np.uint16)
                for corner_orientation in range(CORNER_ORIENTATION_SIZE):
                    cube.set_coord((corner_orientation, 0, 0, (0, 0, 0)))
                    for i in range(len(ALL_MOVES)):
                        trans_corner_orientation[corner_orientation, i] = cube.apply_move(ALL_MOVES[i], cube).get_coord()[0]
                with open("tco.pkl", "wb") as file:
                    pickle.dump(trans_corner_orientation, file)
                self.trans_corner_orientation = trans_corner_orientation

            try:
                with open("teo.pkl", "rb") as file:
                    self.trans_edge_orientation = pickle.load(file)
            except FileNotFoundError:
                trans_edge_orientation = np.zeros((EDGE_ORIENTATION_SIZE, len(ALL_MOVES)), dtype=np.uint16)
                for edge_orientation in range(EDGE_ORIENTATION_SIZE):
                    cube.set_coord((0, edge_orientation, 0, (0, 0, 0)))
                    for i in range(len(ALL_MOVES)):
                        trans_edge_orientation[edge_orientation, i] = cube.apply_move(ALL_MOVES[i], cube).get_coord()[1]
                with open("teo.pkl", "wb") as file:
                    pickle.dump(trans_edge_orientation, file)
                self.trans_edge_orientation = trans_edge_orientation

            try:
                with open("tcp.pkl", "rb") as file:
                    self.trans_corner_permutation = pickle.load(file)
            except FileNotFoundError:
                trans_corner_permutation = np.zeros((CORNER_PERMUTATION_SIZE, len(ALL_MOVES)), dtype=np.uint16)
                for corner_permutation in range(CORNER_PERMUTATION_SIZE):
                    cube.set_coord((0, 0, corner_permutation, (0, 0, 0)))
                    for i in range(len(ALL_MOVES)):
                        trans_corner_permutation[corner_permutation, i] = cube.apply_move(ALL_MOVES[i], cube).get_coord()[2]
                with open("tcp.pkl", "wb") as file:
                    pickle.dump(trans_corner_permutation, file)
                self.trans_corner_permutation = trans_corner_permutation

            try:
                with open("tep.pkl", "rb") as file:
                    self.trans_edge_permutation = pickle.load(file)
            except FileNotFoundError:
                trans_edge_permutation = np.zeros((PARTIAL_EDGE_PERMUTATION_SIZE, len(ALL_MOVES)), dtype=np.uint16)
                for edge_permutation in range(PARTIAL_EDGE_PERMUTATION_SIZE):
                    cube.set_coord((0, 0, 0, (edge_permutation, EMPTY, EMPTY)))
                    for i in range(len(ALL_MOVES)):
                        trans_edge_permutation[edge_permutation, i] = cube.apply_move(ALL_MOVES[i], cube).get_coord()[3][0]
                with open("tep.pkl", "wb") as file:
                    pickle.dump(trans_edge_permutation, file)
                self.trans_edge_permutation = trans_edge_permutation

    def get_next_position(self, position: Cube | tuple, move: str) -> Cube | tuple:
        if self.transition_tables:
            return (self.trans_corner_orientation[position[0], ALL_MOVES_INDEX[move]],
                    self.trans_edge_orientation[position[1], ALL_MOVES_INDEX[move]],
                    self.trans_corner_permutation[position[2], ALL_MOVES_INDEX[move]],
                    tuple([self.trans_edge_permutation[position[3][axis], ALL_MOVES_INDEX[move]] for axis in range(3)]))
        if self.cube.representation == "coord":
            self.cube.set_coord(position)
            return self.cube.apply_move(move).get_coord()
        return position.apply_move(move, position)

    def is_solved(self, position: Cube | tuple) -> bool:
        if self.cube.representation == "coord":
            return position == SOLVED_PARTIAL_COORD
        return repr(position) == SOLVED_REPR

    def solve(self, max_depth: int = 5) -> str:
        solution = []
        position = self.cube.get_coord() if self.cube.representation == "coord" else self.cube
        for depth in range(max_depth + 1):
            if self._solve(depth, solution, position):
                break
        return " ".join(solution[::-1])

    def _solve(self, depth: int, solution: list[str], position: Cube | tuple, last_move: str = None) -> bool:
        if depth == 0:
            return self.is_solved(position)
        if self.pruning_tables:
            if self.prun_corner_orientation[position[0]] <= depth:
                if self.prun_edge_orientation[position[1]] <= depth:
                    if self.prun_corner_permutation[position[2]] <= depth:
                        if self.prun_edge_permutation[0, position[3][0]] <= depth:
                            if self.prun_edge_permutation[1, position[3][1]] <= depth:
                                if self.prun_edge_permutation[2, position[3][2]] <= depth:
                                    for move in NEXT_MOVES[last_move]:
                                        next_position = self.get_next_position(position, move)
                                        if self._solve(depth - 1, solution, next_position, move):
                                            solution.append(move)
                                            return True
            return False
        else:
            for move in NEXT_MOVES[last_move]:
                next_position = self.get_next_position(position, move)
                if self._solve(depth - 1, solution, next_position, move):
                    solution.append(move)
                    return True
            return False
