import os
import pickle
import numpy as np
from collections import deque

from cube_solver import Cube
from cube_solver.constants import ALL_MOVES, ALL_MOVES_INDEX, NEXT_MOVES, EMPTY, NUM_PERM_AXIS
from cube_solver.constants import SOLVED_PARTIAL_COORD, SOLVED_REPR, COORDS_SIZES
from cube_solver.constants import NUM_PHASES, NUM_THREADS, PHASE_MOVES, PHASE_NEXT_MOVES, PHASE_TABLE_SIZES, CORNER_THREAD


def get_pruning_table(name: str, coord_index: int) -> np.ndarray:
    try:
        with open(f"tables/{name}.pkl", "rb") as file:
            return pickle.load(file)
    except FileNotFoundError:
        cube = Cube()
        shape = COORDS_SIZES[coord_index]
        if coord_index == 3:
            shape = (3, shape)
        pruning_table = np.full(shape, EMPTY, dtype=np.int8)
        for i in range(3 if coord_index == 3 else 1):
            cube.reset()
            init_coord = cube.get_coords()[coord_index]
            if coord_index == 3:
                init_coord = init_coord[i]
            if coord_index == 3:
                pruning_table[i, init_coord] = 0
            else:
                pruning_table[init_coord] = 0
            queue = deque([(init_coord, 0)])  # index, depth
            while queue:
                coords = [0, 0, 0, [EMPTY, EMPTY, EMPTY]]
                if coord_index == 3:
                    coords[coord_index][i], depth = queue.popleft()
                else:
                    coords[coord_index], depth = queue.popleft()
                cube.set_coords(coords)
                for move in ALL_MOVES:
                    coord = cube.apply_move(move, cube).get_coords()[coord_index]
                    if coord_index == 3:
                        coord = coord[i]
                        if pruning_table[i, coord] == EMPTY:
                            pruning_table[i, coord] = depth + 1
                            queue.append((coord, depth + 1))
                    else:
                        if pruning_table[coord] == EMPTY:
                            pruning_table[coord] = depth + 1
                            queue.append((coord, depth + 1))
        with open(f"tables/{name}.pkl", "wb") as file:
            pickle.dump(pruning_table, file)
        return pruning_table


def get_transition_table(name: str, coord_index: int) -> np.ndarray:
    try:
        with open(f"tables/{name}.pkl", "rb") as file:
            return pickle.load(file)
    except FileNotFoundError:
        cube = Cube()
        transition_table = np.zeros((COORDS_SIZES[coord_index], len(ALL_MOVES)), dtype=np.uint16)
        for coord in range(COORDS_SIZES[coord_index]):
            coords = [0, 0, 0, [EMPTY, EMPTY, EMPTY]]
            if coord_index == 3:
                coords[coord_index][0] = coord
            else:
                coords[coord_index] = coord
            cube.set_coords(coords)
            for i in range(len(ALL_MOVES)):
                next_coord = cube.apply_move(ALL_MOVES[i], cube).get_coords()[coord_index]
                if coord_index == 3:
                    next_coord = next_coord[0]
                transition_table[coord, i] = next_coord
        with open(f"tables/{name}.pkl", "wb") as file:
            pickle.dump(transition_table, file)
        return transition_table


class Solver:
    def __init__(self, pruning_tables: bool = False, transition_tables: bool = False):
        self.pruning_tables = pruning_tables
        self.transition_tables = transition_tables

        if self.pruning_tables:
            os.makedirs("tables/", exist_ok=True)
            self.prun_corner_orientation = get_pruning_table("pco", 0)
            self.prun_edge_orientation = get_pruning_table("peo", 1)
            self.prun_corner_permutation = get_pruning_table("pcp", 2)
            self.prun_edge_permutation = get_pruning_table("pep", 3)

        if self.transition_tables:
            os.makedirs("tables/", exist_ok=True)
            self.trans_corner_orientation = get_transition_table("tco", 0)
            self.trans_edge_orientation = get_transition_table("teo", 1)
            self.trans_corner_permutation = get_transition_table("tcp", 2)
            self.trans_edge_permutation = get_transition_table("tep", 3)

    def _get_next_position(self, position: Cube | tuple, move: str) -> Cube | tuple:
        if self.transition_tables:
            return (self.trans_corner_orientation[position[0], ALL_MOVES_INDEX[move]],
                    self.trans_edge_orientation[position[1], ALL_MOVES_INDEX[move]],
                    self.trans_corner_permutation[position[2], ALL_MOVES_INDEX[move]],
                    tuple(self.trans_edge_permutation[position[3][axis], ALL_MOVES_INDEX[move]] for axis in range(3)))
        if isinstance(position, tuple):
            cube = Cube()
            cube.set_coords(position, partial_edge=True)
            return cube.apply_move(move).get_coords(partial_edge=True)
        return position.apply_move(move, position)

    def is_solved(self, position: Cube | tuple) -> bool:
        if isinstance(position, tuple):
            return position == SOLVED_PARTIAL_COORD
        return repr(position) == SOLVED_REPR

    def solve(self, cube: Cube, max_depth: int = 10) -> str:
        solution = []
        position = cube.get_coords(partial_edge=True) if cube.representation == "coord" else cube
        for depth in range(max_depth + 1):
            if self._solve(depth, position, solution):
                break
        return " ".join(solution[::-1])

    def _solve(self, depth: int, position: Cube | tuple, solution: list[str], last_move: str = None) -> bool:
        if depth == 0:
            return self.is_solved(position)
        if self.pruning_tables:
            if self.prun_corner_permutation[position[2]] <= depth:
                if self.prun_edge_permutation[0, position[3][0]] <= depth:
                    if self.prun_edge_permutation[1, position[3][1]] <= depth:
                        if self.prun_edge_permutation[2, position[3][2]] <= depth:
                            if self.prun_corner_orientation[position[0]] <= depth:
                                if self.prun_edge_orientation[position[1]] <= depth:
                                    for move in NEXT_MOVES[last_move]:
                                        next_position = self._get_next_position(position, move)
                                        if self._solve(depth - 1, next_position, solution, move):
                                            solution.append(move)
                                            return True
            return False
        else:
            for move in NEXT_MOVES[last_move]:
                next_position = self._get_next_position(position, move)
                if self._solve(depth - 1, next_position, solution, move):
                    solution.append(move)
                    return True
            return False

    def thistlethwaite(self, cube: Cube) -> str:
        self.phase_tables = [self._get_phase_table(phase) for phase in range(NUM_PHASES)]

        solution = []
        for phase in range(NUM_PHASES):
            phase_solution = []
            coords = cube.get_coords(partial_corner=True, partial_edge=True)
            phase_coords = self._get_phase_coords(phase, coords)
            depth = self.phase_tables[phase][phase_coords]
            if self._phase(phase, depth, cube, phase_solution):
                phase_solution = " ".join(phase_solution[::-1])
                solution.append(phase_solution)
                cube.apply_maneuver(phase_solution)
            else:
                return "NO SOLUTION FOUND"

        return " ".join(solution)

    def _get_phase_coords(self, phase: int, coords: tuple) -> tuple:
        if phase == 0:
            edge_orientation = coords[1]
            return (edge_orientation,)
        if phase == 1:
            corner_orientation = coords[0]
            edge_combination = coords[3][0] // NUM_PERM_AXIS
            return (corner_orientation, edge_combination)
        if phase == 2:
            corner_combination = coords[2][0] // NUM_PERM_AXIS
            edge_combination = coords[3][2] // NUM_PERM_AXIS
            corner_thread = CORNER_THREAD[coords[2][0] % NUM_PERM_AXIS, coords[2][1] % NUM_PERM_AXIS]
            return (corner_combination, edge_combination, corner_thread)
        if phase == 3:
            corner_permutation = [cp % NUM_PERM_AXIS for cp in coords[2]]
            corner_permutation[1] //= NUM_THREADS
            edge_permutation = [ep % NUM_PERM_AXIS for ep in coords[3]]
            edge_permutation[2] //= 2  # parity
            return tuple(corner_permutation) + tuple(edge_permutation)

    def _get_phase_table(self, phase):
        try:
            with open(f"tables/phase{phase+1}.pkl", "rb") as file:
                return pickle.load(file)
        except FileNotFoundError:
            cube = Cube()
            shape = PHASE_TABLE_SIZES[phase]
            phase_table = np.full(shape, EMPTY, dtype=np.int8)
            coords = cube.get_coords(partial_corner=True, partial_edge=True)
            phase_coords = self._get_phase_coords(phase, coords)
            phase_table[phase_coords] = 0
            queue = deque([(coords, 0)])
            while queue:
                coords, depth = queue.popleft()
                cube.set_coords(coords, partial_corner=True, partial_edge=True)
                for move in PHASE_MOVES[phase]:
                    next_cube = cube.apply_move(move, cube)
                    next_coords = next_cube.get_coords(partial_corner=True, partial_edge=True)
                    phase_coords = self._get_phase_coords(phase, next_coords)
                    if phase_table[phase_coords] == EMPTY:
                        phase_table[phase_coords] = depth + 1
                        queue.append((next_coords, depth + 1))
            with open(f"tables/phase{phase+1}.pkl", "wb") as file:
                pickle.dump(phase_table, file)
            return phase_table

    def _phase(self, phase: int, depth: int, cube: Cube, solution: list[str], last_move: str = None) -> bool:
        if depth == 0:
            return True
        for move in PHASE_NEXT_MOVES[phase][last_move]:
            next_cube = cube.apply_move(move, cube)
            coords = next_cube.get_coords(partial_corner=True, partial_edge=True)
            phase_coords = self._get_phase_coords(phase, coords)
            next_depth = self.phase_tables[phase][phase_coords]
            if next_depth < depth:
                if self._phase(phase, next_depth, next_cube, solution, move):
                    solution.append(move)
                    return True
                return False
        return False
