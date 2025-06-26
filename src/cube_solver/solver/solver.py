"""BaseSolver module."""
import numpy as np
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
    """Number of phases of the solving algorithm."""
    partial_corner_perm: bool
    """Whether the solving algorithm uses the normal or the partial corner permutation. """
    partial_edge_perm: bool
    """Whether the solving algorithm uses the normal or the partial edge permutation. """
    transition_kwargs: list[TransitionDef]
    """List of transition table definitions."""
    pruning_kwargs: list[list[PruningDef]] = []
    """List of pruning table definitions for each phase."""
    phase_moves: list[list[Move]] = []
    """List of available moves for each phase."""

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

    def __init__(self, use_pruning_tables: bool = True, use_transition_tables: bool = True):
        """
        Create :class:`.BaseSolver` object.

        Parameters
        ----------
        use_pruning_tables : bool, optional
            Whether to use pruning tables to reduce the tree search space.
            If ``True``, creates or loads the tables from the ``tables/`` directory.
            Default is ``True``.
        use_transition_tables : bool, optional
            Whether to use transition tables for cube state transitions.
            If ``True``, creates or loads the tables from the ``tables/`` directory.
            Default is ``True``.
        """
        if not isinstance(use_pruning_tables, bool):
            raise TypeError(f"use_pruning_tables must be bool, not {type(use_pruning_tables).__name__}")
        if not isinstance(use_transition_tables, bool):
            raise TypeError(f"use_transition_tables must be bool, not {type(use_transition_tables).__name__}")
        if use_pruning_tables and not self.pruning_kwargs:
            raise ValueError("cannot use pruning tables with empty class attribute 'pruning_kwargs'")

        self.cube: Cube = Cube()
        """Internal :class:`cube_solver.Cube` object. Used to generate transition and pruning tables."""
        init_coords = self.flatten(self.get_coords(self.cube))
        self.solved_coords: list[FlattenCoords] = [self.phase_coords(init_coords, phase) for phase in range(self.num_phases)]
        """Solved flatten coordinates for each phase."""

        self.use_pruning_tables: bool = use_pruning_tables
        """Whether to use pruning tables to reduce the tree search space."""
        self.use_transition_tables: bool = use_transition_tables
        """Whether to use transition tables for cube state transitions."""
        self.transition_tables: dict[str, np.ndarray] = {}
        """Transition tables used to compute cube state transitions."""
        if self.use_transition_tables:  # TODO test with different extension
            self.transition_tables = utils.get_tables("transition.npz", self.transition_kwargs,
                                                      self.generate_transition_table, accumulate=True)
        if not self.phase_moves:
            self.phase_moves = [[*Move.face_moves()] * self.num_phases]
        self.pruning_tables: dict[str, np.ndarray] = {}
        """Pruning tables used to reduce the tree search space."""
        if self.use_pruning_tables:
            pruning_kwargs = []
            for phase, phase_kwargs in enumerate(self.pruning_kwargs):
                for kwargs in phase_kwargs:
                    kwargs.phase = phase
                    pruning_kwargs.append(kwargs)
            pruning_filename = f"pruning_{self.__class__.__name__.lower()}.npz"
            self.pruning_tables = utils.get_tables(pruning_filename, pruning_kwargs,
                                                   self.generate_pruning_table, accumulate=False)

        self.final_moves: list[set[Move]] = []
        """Final allowed moves for each phase except the last."""
        for i in range(self.num_phases - 1):
            self.final_moves.append({Move.NONE} | set(self.phase_moves[i]) - set(self.phase_moves[i+1]))

        self.next_moves: list[dict[Move, list[Move]]] = []
        """Allowed next moves based on the previous move, for each phase."""
        for phase_moves in self.phase_moves:
            next_moves = {Move.NONE: phase_moves}
            next_moves.update({move: [mv for mv in NEXT_MOVES[move] if mv in phase_moves] for move in phase_moves})
            self.next_moves.append(next_moves)

        self.nodes: list[list[list[int]]]  # TODO needed?
        """Number of visited nodes during a solve."""
        self.solve_checks: list[list[list[int]]]  # TODO needed?
        """Number of solve checks during a solve."""
        self.prune_checks: list[list[list[int]]]
        """Number of prune checks during a solve."""  # TODO add more?

    # TODO get transition table for each phase
    # TODO add log: generating transition tables?, move to utils
    def generate_transition_table(self, coord_name: str, coord_size: int) -> np.ndarray:
        """
        Generate the cube coordinate transition table.

        Parameters
        ----------
        coord_name : {'co', 'eo', 'cp', 'ep', 'pcp', 'pep'}
            Cube coordinate name.
        coord_size : int
            Size of the cube coordinate.

        Returns
        -------
        transition_table : ndarray
            Cube coordinate transition table.
        """
        if coord_size - 1 > np.iinfo(np.uint16).max:
            raise ValueError(f"coord_size must be <= {np.iinfo(np.uint16).max + 1} (got {coord_size})")
        transition_table = np.zeros((coord_size, len(NEXT_MOVES[Move.NONE])), dtype=np.uint16)
        for coord in range(coord_size):
            self.cube.set_coord(coord_name, coord)
            transition_table[coord] = [apply_move(self.cube, Move(mv)).get_coord(coord_name) for mv in NEXT_MOVES[Move.NONE]]
        return transition_table

    # TODO maybe move to utils when having transition tables per phase
    def generate_pruning_table(self, phase: int, shape: int | tuple[int, ...],
                               indexes: int | tuple[int, ...] | None, **kwargs) -> np.ndarray:
        """
        Generate the phase coordinates pruning table.

        Parameters
        ----------
        phase : int
            Solver phase (0-indexed).
        shape : int or tuple of int
            Shape of the pruning table.
        indexes : int or tuple of int or None
            Index or indexes of the phase coordinates to use for the pruning table.
            If ``None``, use all the phase coordinates.

        Returns
        -------
        pruning_table : ndarray
            Phase coordinates pruning table.
        """
        self.cube.reset()
        coords = self.get_coords(self.cube)
        phase_coords = self.phase_coords(self.flatten(coords), phase)
        prune_coords = self.prune_coords(phase_coords, indexes)
        pruning_table = np.full(shape, NONE, dtype=np.int8)
        pruning_table[prune_coords] = 0
        queue = deque([(coords, 0)])  # TODO pass the pruning coords and have transition tables for that
        while queue:
            coords, depth = queue.popleft()
            for move in self.phase_moves[phase]:
                next_coords = self.next_position(coords, Move(move))
                phase_coords = self.phase_coords(self.flatten(next_coords), phase)
                prune_coords = self.prune_coords(phase_coords, indexes)
                if pruning_table[prune_coords] == NONE:
                    pruning_table[prune_coords] = depth + 1
                    queue.append((next_coords, depth + 1))
        return pruning_table

    def get_coords(self, cube: Cube) -> CoordsType:
        """
        Get cube coordinates.

        Get the `corner orientation`, `edge orientation`,
        `(partial) corner permutation` and `(partial) edge permutation` coordinates,
        according to :attr:`partial_corner_perm` and :attr:`partial_edge_perm`.

        Parameters
        ----------
        cube : Cube
            Cube object.

        Returns
        -------
        coords : tuple of (int or tuple of int)
            Cube coordinates in the following order:
            `corner orientation`, `edge orientation`, `(partial) corner permutation`, `(partial) edge permutation`,
            according to :attr:`partial_corner_perm` and :attr:`partial_edge_perm`.

        See Also
        --------
        cube_solver.Cube.get_coords
        """
        return cube.get_coords(self.partial_corner_perm, self.partial_edge_perm)

    def set_coords(self, cube: Cube, coords: CoordsType):
        """
        Set cube coordinates.

        Set the `corner orientation`, `edge orientation`,
        `(partial) corner permutation` and `(partial) edge permutation` coordinates,
        according to :attr:`partial_corner_perm` and :attr:`partial_edge_perm`.

        Parameters
        ----------
        cube : Cube
            Cube object.
        coords : tuple of (int or tuple of int)
            Cube coordinates in the following order:
            `corner orientation`, `edge orientation`, `(partial) corner permutation`, `(partial) edge permutation`,
            according to :attr:`partial_corner_perm` and :attr:`partial_edge_perm`.

        See Also
        --------
        cube_solver.Cube.set_coords
        """
        cube.set_coords(coords, self.partial_corner_perm, self.partial_edge_perm)

    def flatten(self, coords: CoordsType) -> FlattenCoords:
        """
        Get the flatten cube coordinates.

        Parameters
        ----------
        coords : tuple of (int or tuple of int)
            Cube coordinates.

        Returns
        -------
        flatten_coords : tuple of int
            Flatten cube coordinates.
        """
        flatten = ()
        for coord in coords:
            flatten += (coord,) if isinstance(coord, int) else coord
        return flatten

    @abstractmethod
    def phase_coords(self, coords: FlattenCoords, phase: int) -> FlattenCoords:
        """
        Get the coordinates for the specified phase.

        Parameters
        ----------
        coords : tuple of int
            Flatten cube coordinates.
        phase : int
            Solver phase (0-indexed).

        Returns
        -------
        phase_coords : tuple of int
            Phase coordinates.
        """

    def prune_coords(self, phase_coords: FlattenCoords, indexes: int | tuple[int, ...] | None) -> FlattenCoords:
        """
        Get the prune coordinates.

        Parameters
        ----------
        phase_coords : tuple of int
            Phase coordinates.
        indexes : int or tuple of int or None
            Index or indexes of the phase coordinates to use as the prune coordinates.
            If ``None``, use all the phase coordinates.

        Returns
        -------
        prune_coords : tuple of int
            Prune coordinates.
        """
        if indexes is None:
            return phase_coords
        if isinstance(indexes, int):
            return (phase_coords[indexes],)
        return tuple(phase_coords[index] for index in indexes)

    # TODO add examples?
    def is_solved(self, position: Cube | CoordsType, phase: int) -> bool:
        """
        Whether the position is solved at the specified phase.

        Parameters
        ----------
        position : Cube or tuple of (int or tuple of int)
            Cube object or cube coordinates to check.
        phase : int
            Solver phase (0-indexed).

        Returns
        -------
        bool
            ``True`` if the cube position is solved, ``False`` otherwise.
        """
        if isinstance(position, Cube):
            position = self.get_coords(position)
        phase_coords = self.phase_coords(self.flatten(position), phase)
        return phase_coords == self.solved_coords[phase]

    def prune(self, position: Cube | CoordsType, phase: int, depth: int) -> bool:
        """
        Whether to prune the search tree.

        Checks whether the current ``depth`` meets the lower bound specified by the :attr:`pruning_tables`.

        Parameters
        ----------
        position : Cube or tuple of (int or tuple of int)
            Cube object or cube coordinates to check.
        depth : int
            Current search depth.
        phase : int
            Solver phase (0-indexed).

        Returns
        -------
        bool
            ``True`` if the search tree should be pruned, ``False`` otherwise.
        """
        if self.pruning_tables:
            if isinstance(position, Cube):
                position = self.get_coords(position)
            phase_coords = self.phase_coords(self.flatten(position), phase)
            for kwargs in self.pruning_kwargs[phase]:
                prune_coords = self.prune_coords(phase_coords, kwargs.indexes)
                if self.pruning_tables[kwargs.name][prune_coords] > depth:
                    return True
        return False

    @overload
    def next_position(self, position: Cube, move: Move) -> Cube: ...
    @overload
    def next_position(self, position: CoordsType, move: Move) -> CoordsType: ...

    def next_position(self, position: Cube | CoordsType, move: Move) -> Cube | CoordsType:
        """
        Get the next cube position.

        Parameters
        ----------
        position : Cube or tuple of (int or tuple of int)
            Cube object or cube coordinates.
        move : Move
            Move to apply.

        Returns
        -------
        next_position : Cube or tuple of (int or tuple of int)
            Cube object or cube coordinates with the move applied.
        """
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

    # TODO add timeout for optimal or for all
    def solve(self, cube: Cube, max_length: int | None = None, optimal: bool = False, verbose: int = 0) -> Maneuver | list[Maneuver] | None:  # TODO cube could be None for perf testing?
        """
        Solve the cube position.

        Parameters
        ----------
        cube : Cube
            Cube object to be solved.
        max_length : int or None, optional
            Maximum number of moves to search. If ``None``, search indefinitely. Default is ``None``.
        optimal : bool, optimal
            If ``True``, finds the optimal solution. Default is ``False``.
        verbose : {0, 1}, optional
            Verbosity level. Default is ``0``.

            * ``0`` returns only the solution.
            * ``1`` returns the solution for each phase.
              Also logs all solutions found if ``optimal`` is ``True``.

        Returns
        -------
        solution : Maneuver or list of Maneuver or None
            Solution for the cube position,
            or list of solutions for each phase if ``verbose`` is ``1``,
            or ``None`` if no solution is found.

        Examples
        --------
        >>> from cube_solver import Cube, Kociemba
        >>> cube = Cube("U F R")
        >>> solver = Kociemba()
        >>> solver.solve(cube)
        R' F' U'
        >>> solver.solve(cube, optimal=True, verbose=1)
        R F U
        """
        if not isinstance(cube, Cube):  # TODO add test for optimal and verbose
            raise TypeError(f"cube must be Cube, not {type(cube).__name__}")
        if max_length is not None and not isinstance(max_length, int):
            raise TypeError(f"max_length must be int or None, not {type(max_length).__name__}")
        if not isinstance(optimal, bool):
            raise TypeError(f"optimal must be bool, not {type(optimal).__name__}")
        if not isinstance(verbose, int):
            raise TypeError(f"verbose must be int, not {type(verbose).__name__}")
        if isinstance(max_length, int) and max_length < 0:
            raise ValueError(f"max_length must be >= 0 (got {max_length})")
        if verbose not in (0, 1):
            raise ValueError(f"verbose must be 0 or 1 (got {verbose})")
        if cube.permutation_parity is None:
            raise ValueError("invalid cube state")
        try:
            self.get_coords(cube)
        except ValueError:
            raise ValueError("invalid cube state")

        # self.nodes = []
        # self.checks = []
        # TODO order
        self.max_length = None if optimal else max_length
        self.optimal = optimal  # TODO add to __init__
        self.verbose = verbose
        self.solution = [[] * self.num_phases]  # TODO check

        # self.nodes.append([[] for i in range(self.num_phases)])
        # self.checks.append([[] for i in range(self.num_phases)])  # TODO try not to erase after each new depth
        position = self.get_coords(cube) if self.use_transition_tables else cube
        if self.phase_search(position):
            return Maneuver([move for phase in range(self.num_phases) for move in [*self.solution[phase]][::-1]])
        # if optimal:
        #     return Maneuver([move for phase in range(self.num_phases) for move in [*self.best_solution[phase]][::-1]])
        # if x is not None:
        #     return Maneuver([move for phase in range(self.num_phases) for move in [*self.solution[phase]][::-1]])
        return None

    # TODO make iterative version and compare
    # TODO test starting with phase > 0
    def phase_search(self, position: Cube | CoordsType, phase: int = 0, current_length: int = 0) -> bool:
        """
        Solve the cube position from the specified phase.

        Parameters
        ----------
        position : Cube or tuple of (int or tuple of int)
            Cube object or cube coordinates to search.
        phase : int, optional
            Phase to solve from. Default is ``0``.
        current_length : int, optional
            Current solution length. Default is ``0``.
        last_move : Move, optional
            Last move performed on the cube. Helps avoid repeating the same move during the phase search.
            :attr:`.Move.NONE` indicates that no moves have been made yet (i.e. the start of the phase search).
            Default is :attr:`.Move.NONE`.

        Returns
        -------
        bool
            `True` if a solution was found, `False` otherwise.
        """
        if phase == self.num_phases:
            # if self.optimal:
            #     self.skip_phase = phase - 1
            #     self.max_depth = current_depth - 1
            #     self.best_solution = deepcopy(self.solution)
            return True
        depth = 0
        # self.nodes[phase].append([])
        while True if self.max_length is None else current_length + depth <= self.max_length:
            # self.nodes[phase][-1].append(0)
            self.solution[phase].append(Move.NONE)
            if self.search(position, phase, depth, current_length):
                return True
            # if length is None:
            #     self.solution[phase] = deque([])
            #     return False
            depth += 1
        self.solution[phase] = []
        return False

    def search(self, position: Cube | CoordsType, phase: int, depth: int, current_length: int, last_move: Move = Move.NONE) -> bool:
        """
        Parameters
        ----------
        last_move : Move, optional
            Last move performed on the cube. Helps avoid repeating the same move during the phase search.
            :attr:`.Move.NONE` indicates that no moves have been made yet (i.e. the start of the phase search).
            Default is :attr:`.Move.NONE`.
        """
        # self.nodes[phase][-1][-1] += 1
        self.solution[phase][depth] = last_move
        if depth == 0:  # TODO change this inside the pruning and check stats, self.num_prunes?
            if phase == self.num_phases - 1 or (last_move in self.final_moves[phase]):
                # self.solve_checks[phase][-1][-1] += 1
                if self.is_solved(position, phase):
                    return self.phase_search(position, phase + 1, current_length)
            return False
        # self.prune_checks[phase][-1][-1] += 1
        if not self.prune(position, phase, depth):
            for move in self.next_moves[phase][last_move]:  # TODO get last move from solution
                next_position = self.next_position(position, move)
                if self.search(next_position, phase, depth - 1, current_length + 1, move):
                    return True
                # if self.optimal:
                #     if phase == self.skip_phase:
                #         return None
                #     elif self.skip_phase != -1:
                #         self.skip_phase = -1
                #     continue
        return False
