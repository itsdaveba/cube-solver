import numpy as np
from pathlib import Path
from copy import deepcopy
from typing import overload
from collections import deque
from abc import ABC, abstractmethod

from ..defs import CoordsType, NEXT_MOVES
from ..cube import Cube, Move, Maneuver, apply_move
from ..cube.defs import CORNER_ORIENTATION_SIZE, EDGE_ORIENTATION_SIZE
from ..cube.defs import CORNER_PERMUTATION_SIZE, EDGE_PERMUTATION_SIZE
from ..cube.defs import PARTIAL_CORNER_PERMUTATION_SIZE, PARTIAL_EDGE_PERMUTATION_SIZE
from .defs import NONE, FlattenCoords, TransitionDef, PruningDef
from . import utils

MOVE_TO_INDEX = np.zeros(max(len(Move), max(Move) + 1), dtype=int)
for cubie in NEXT_MOVES[Move.NONE]:
    MOVE_TO_INDEX[cubie] = NEXT_MOVES[Move.NONE].index(cubie)


class BaseSolver(ABC):
    num_phases: int = 1
    partial_corner_perm: bool
    partial_edge_perm: bool
    transition_kwargs: list[TransitionDef]
    pruning_kwargs: list[list[PruningDef]] = [[]]
    phase_moves: list[list[Move]] = [[*Move.face_moves()]]

    def __init_subclass__(cls):
        for attr in ("partial_corner_perm", "partial_edge_perm"):
            if not hasattr(cls, attr):
                raise AttributeError(f"'{cls.__name__}' class must define class attribute '{attr}'")

        cls.transition_kwargs = [
            TransitionDef(coord_name="co", coord_size=CORNER_ORIENTATION_SIZE),
            TransitionDef(coord_name="eo", coord_size=EDGE_ORIENTATION_SIZE),
            TransitionDef(coord_name="pcp" if cls.partial_corner_perm else "cp",
                          coord_size=PARTIAL_CORNER_PERMUTATION_SIZE if cls.partial_corner_perm else CORNER_PERMUTATION_SIZE),
            TransitionDef(coord_name="pep" if cls.partial_edge_perm else "ep",
                          coord_size=PARTIAL_EDGE_PERMUTATION_SIZE if cls.partial_edge_perm else EDGE_PERMUTATION_SIZE)]

    def __init__(self, use_transition_tables: bool = False):
        if not isinstance(use_transition_tables, bool):
            raise TypeError(f"use_transition_tables must be bool, not {type(use_transition_tables).__name__}")

        # TODO order
        self.next_moves: list[dict[Move, list[Move]]] = []
        for phase_moves in self.phase_moves:
            next_moves = {Move.NONE: phase_moves}
            next_moves.update({move: [mv for mv in NEXT_MOVES[move] if mv in phase_moves] for move in phase_moves})
            self.next_moves.append(next_moves)

        self.final_moves: list[set[Move]] = []
        for i in range(self.num_phases - 1):
            self.final_moves.append({Move.NONE} | set(self.phase_moves[i]) - set(self.phase_moves[i+1]))

        self.cube: Cube = Cube()
        init_coords = self.flatten(self.get_coords(self.cube))
        self.solved_coords: list[FlattenCoords] = [self.phase_coords(phase, init_coords) for phase in range(self.num_phases)]
        self.use_transition_tables: bool = use_transition_tables  # TODO trans tables for each phase?
        self.transition_tables: dict[str, np.ndarray] = {}
        if self.use_transition_tables:  # TODO test with different extension
            self.transition_tables = utils.get_tables("transition.npz", self.transition_kwargs,
                                                      self.generate_transition_table, accumulate=True)
        pruning_kwargs = []
        for phase, phase_kwargs in enumerate(self.pruning_kwargs):
            for kwargs in phase_kwargs:
                kwargs.phase = phase
                pruning_kwargs.append(kwargs)
        pruning_filename = f"pruning_{self.__class__.__name__.lower()}.npz"
        self.pruning_tables: dict[str, np.ndarray] = utils.get_tables(pruning_filename, pruning_kwargs,
                                                                      self.generate_pruning_table, accumulate=False)

        self.nodes: list[list[list[int]]]  # TODO needed?
        self.checks: list[list[list[int]]]  # TODO needed?

    # TODO add log: generating transition tables?, move to utils
    def generate_transition_table(self, coord_name: str, coord_size: int) -> np.ndarray:
        if coord_size - 1 > np.iinfo(np.uint16).max:
            raise ValueError("")
        transition_table = np.zeros((coord_size, len(NEXT_MOVES[Move.NONE])), dtype=np.uint16)
        for coord in range(coord_size):
            self.cube.set_coord(coord_name, coord)
            transition_table[coord] = [apply_move(self.cube, Move(mv)).get_coord(coord_name) for mv in NEXT_MOVES[Move.NONE]]
        return transition_table

    # TODO maybe move to utils when having transition tables per phase
    def generate_pruning_table(self, phase: int, shape: int | tuple[int, ...],
                               indexes: int | tuple[int, ...] | None, **_) -> np.ndarray:
        self.cube.reset()
        coords = self.get_coords(self.cube)
        phase_coords = self.phase_coords(phase, self.flatten(coords))
        prune_coords = self.prune_coords(phase_coords, indexes)
        pruning_table = np.full(shape, NONE, dtype=np.int8)
        pruning_table[prune_coords] = 0
        queue = deque([(coords, 0)])  # TODO pass the pruning coords and have transition tables for that
        while queue:
            coords, depth = queue.popleft()
            for move in self.phase_moves[phase]:
                next_coords = self.next_position(coords, Move(move))
                phase_coords = self.phase_coords(phase, self.flatten(next_coords))
                prune_coords = self.prune_coords(phase_coords, indexes)
                if pruning_table[prune_coords] == NONE:
                    pruning_table[prune_coords] = depth + 1
                    queue.append((next_coords, depth + 1))
        return pruning_table

    def get_coords(self, cube: Cube) -> CoordsType:
        return cube.get_coords(self.partial_corner_perm, self.partial_edge_perm)

    def set_coords(self, cube: Cube, coords: CoordsType):
        cube.set_coords(coords, self.partial_corner_perm, self.partial_edge_perm)

    def flatten(self, coords: CoordsType) -> FlattenCoords:
        flatten = ()
        for coord in coords:
            flatten += (coord,) if isinstance(coord, int) else coord
        return flatten

    @abstractmethod
    def phase_coords(self, phase: int, coords: FlattenCoords) -> FlattenCoords: ...

    def prune_coords(self, phase_coords: FlattenCoords, indexes: int | tuple[int, ...] | None) -> FlattenCoords:
        if indexes is None:
            return phase_coords
        if isinstance(indexes, int):
            return (phase_coords[indexes],)
        return tuple(phase_coords[index] for index in indexes)

    @overload
    def next_position(self, position: Cube, move: Move) -> Cube: ...
    @overload
    def next_position(self, position: CoordsType, move: Move) -> CoordsType: ...

    def next_position(self, position: Cube | CoordsType, move: Move) -> Cube | CoordsType:
        if isinstance(position, Cube):
            return apply_move(position, move)
        if self.use_transition_tables:
            next_position = ()
            for coord, kwargs in zip(position, self.transition_kwargs):
                if isinstance(coord, int):
                    next_position += (self.transition_tables[kwargs.coord_name][coord, MOVE_TO_INDEX[move]].item(),)
                else:
                    next_position += (tuple(self.transition_tables[kwargs.coord_name][coord, MOVE_TO_INDEX[move]].tolist()),)
            return next_position
        cube = Cube()  # TODO use self?
        self.set_coords(cube, position)
        cube.apply_move(move)
        return self.get_coords(cube)

    def is_solved(self, position: Cube | CoordsType, phase: int) -> bool:
        """
        Check whether the `cube` position is solved at the current `phase`.

        Parameters
        ----------
        cube : Cube
            Cube object to check.
        phase : int or None, optional
            Phase to check (0-indexed). If `None`, checks whether the cube is fully solved.

        Returns
        -------
        bool
            `True` if the cube position is solved, `False` otherwise.

        Examples
        --------
        >>> from cube_solver import Cube, Kociemba
        >>> cube = Cube("U F2 R2")
        >>> solver = Kociemba()
        >>> solver.is_solved(cube)
        False
        >>> solver.is_solved(cube, phase=0)
        True
        """
        if isinstance(position, Cube):
            position = self.get_coords(position)
        position = self.phase_coords(phase, self.flatten(position))
        return position == self.solved_coords[phase]

    def solve(self, cube: Cube, max_depth: int | None = None, optimal: bool = False, verbose: int = 0) -> Maneuver | None:  # TODO cube could be None for perf testing?
        """
        Solve the `cube` position.

        Parameters
        ----------
        cube : Cube
            Cube object to be solved.
        max_depth : int or None, optional
            The maximum number of moves to search. If `None`, search indefinitely.
        verbose : {0, 1, 2}, optional
            Verbosity level. If `0`, only return the solution.
            If `1`, return the solution with the move count.
            If `2`, return the solution separated by phase and move count for each phase.

        Returns
        -------
        solution : str or None
            Solution for the cube position, or `None` if no solution is found.

        Examples
        --------
        >>> from cube_solver import Cube, Kociemba
        >>> cube = Cube("U F R")
        >>> solver = Kociemba()
        >>> solver.solve(cube)
        R' F' U'
        >>> solver.solve(cube, verbose=1)
        R' F' U' (3)
        >>> solver.solve(cube, verbose=2)
        R' F' (2) | U' (1)
        """
        if not isinstance(cube, Cube):
            raise TypeError(f"cube must be Cube, not {type(cube).__name__}")
        if max_depth is not None and not isinstance(max_depth, int):
            raise TypeError(f"max_depth must be int or None, not {type(max_depth).__name__}")
        if isinstance(max_depth, int) and max_depth < 0:
            raise ValueError(f"max_depth must be >= 0 (got {max_depth})")
        if cube.permutation_parity is None:
            raise ValueError("invalid cube state")

        try:
            # self.nodes = []
            # self.checks = []
            self.optimal = optimal
            self.solution = [deque([]) for i in range(self.num_phases)]
            position = self.get_coords(cube) if self.use_transition_tables else cube
            self.skip_phase = -1
            self.max_depth = None if optimal else max_depth

            # self.nodes.append([[] for i in range(self.num_phases)])
            # self.checks.append([[] for i in range(self.num_phases)])  # TODO try not to erase after each new depth
            x = self.search(position)  # TODO change phase parameter at the end?
            if optimal:
                return Maneuver([move for phase in range(self.num_phases) for move in [*self.best_solution[phase]][::-1]])
            if x is not None:
                return Maneuver([move for phase in range(self.num_phases) for move in [*self.solution[phase]][::-1]])
        except ValueError:
            raise ValueError("invalid cube state")
        return None

    # TODO make iterative version and compare
    def search(self, position: Cube | CoordsType, phase: int = 0, current_depth: int = 0, last_move: Move = Move.NONE) -> bool | None:
        """
        Solve the `cube` position at the current `phase` and `depth`.

        Parameters
        ----------
        cube : Cube
            Cube object to be solved.
        depth : int
            Search depth.
        phase : int, optional
            Phase to solve.
        solution : deque
            If a soluition is foud, stores the sequence of moves for the solution.
        last_face : Face, optional
            Last face move performed on the cube. This helps avoid repeating the same move during the search.
            `Face.NONE` indicates that no face moves have been made yet (i.e. the first move).

        Returns
        -------
        bool
            `True` if a solution was found, `False` otherwise.
        """
        if phase == self.num_phases:
            if self.optimal:
                self.skip_phase = phase - 1
                self.max_depth = current_depth - 1  # TODO calculate solution length
                self.best_solution = deepcopy(self.solution)
            return True
        phase_depth = 0
        # self.nodes[phase].append([])
        while True if self.max_depth is None else current_depth + phase_depth <= self.max_depth:
            # self.nodes[phase][-1].append(0)
            self.solution[phase].appendleft(Move.NONE)
            length = self._search(position, phase, phase_depth, current_depth)
            if length is None:
                self.solution[phase] = deque([])
                return None
            if length:
                return True
            phase_depth += 1
        self.solution[phase] = deque([])
        return None

    def _search(self, position: Cube | CoordsType, phase: int, phase_depth: int, current_depth: int, last_move: Move = Move.NONE) -> bool | None:
        # self.nodes[phase][-1][-1] += 1
        self.solution[phase][phase_depth] = last_move
        if phase_depth == 0:  # TODO change this inside the pruning and check stats, self.num_prunes?
            if phase == self.num_phases - 1 or (last_move in self.final_moves[phase]):  # TODO move this condition inside?
                if self.is_solved(position, phase):
                    return self.search(position, phase + 1, current_depth)
            return False  # TODO simplify, maybe just return self.search?
        if not self.prune(position, phase, phase_depth):
            for move in self.next_moves[phase][last_move]:  # TODO get last move from solution
                next_position = self.next_position(position, move)
                x = self._search(next_position, phase, phase_depth - 1, current_depth + 1, move)
                if self.optimal:
                    if phase == self.skip_phase:
                        return None
                    elif self.skip_phase != -1:
                        self.skip_phase = -1
                    continue
                if x:
                    return True
        return False

    def prune(self, position: Cube | CoordsType, phase: int, depth: int) -> bool:
        """
        Prune the search tree.

        Check whether the current `depth` meets the lower bound specified by the `pruning tables`.

        Parameters
        ----------
        cube : Cube
            Cube object to check.
        depth : int
            Current search depth.
        phase : int
            Phase to prune, used to select the appropiate `pruning tables`.

        Returns
        -------
        bool
            `True` if the search tree should be pruned, `False` otherwise.
        """
        if self.pruning_tables:
            if isinstance(position, Cube):
                position = self.get_coords(position)
            phase_coords = self.phase_coords(phase, self.flatten(position))
            for kwargs in self.pruning_kwargs[phase]:
                prune_coords = self.prune_coords(phase_coords, kwargs.indexes)
                if self.pruning_tables[kwargs.name][prune_coords] > depth:
                    return True
        return False
