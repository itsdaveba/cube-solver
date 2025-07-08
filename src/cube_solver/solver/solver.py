"""BaseSolver module."""
import time
import numpy as np
from copy import deepcopy
from typing import overload
from abc import ABC, abstractmethod
from typing import Union, List, Dict, Set

from ..logger import logger
from ..defs import CoordsType, NEXT_MOVES
from ..cube.enums import Move
from ..cube.maneuver import Maneuver
from ..cube.cube import Cube, apply_move
from ..cube.defs import CORNER_ORIENTATION_SIZE, EDGE_ORIENTATION_SIZE
from ..cube.defs import CORNER_PERMUTATION_SIZE, EDGE_PERMUTATION_SIZE
from ..cube.defs import PARTIAL_CORNER_PERMUTATION_SIZE, PARTIAL_EDGE_PERMUTATION_SIZE
from .defs import FlattenCoords, TransitionDef, PruningDef
from . import utils

MOVE_TO_INDEX = np.zeros(max(len(Move), max(Move) + 1), dtype=int)
for cubie in NEXT_MOVES[Move.NONE]:
    MOVE_TO_INDEX[cubie] = NEXT_MOVES[Move.NONE].index(cubie)


class BaseSolver(ABC):
    num_phases: int = 1
    """Number of phases of the solving algorithm."""
    partial_corner_perm: bool
    """Whether the solving algorithm uses the normal or the partial corner permutation."""
    partial_edge_perm: bool
    """Whether the solving algorithm uses the normal or the partial edge permutation."""
    phase_moves: List[List[Move]]
    """Available moves for each phase."""
    transition_defs: List[TransitionDef]
    """Transition table definitions for each phase."""
    pruning_defs: List[List[PruningDef]]
    """Pruning table definitions for each phase."""

    def __init_subclass__(cls):
        for attr in ("partial_corner_perm", "partial_edge_perm"):
            if not hasattr(cls, attr):
                raise AttributeError(f"'{cls.__name__}' class must define class attribute '{attr}'")

        if not hasattr(cls, "phase_moves"):
            cls.phase_moves = [[*Move.face_moves()] * cls.num_phases]

        cls.transition_defs = [
            TransitionDef(coord_name="co", coord_size=CORNER_ORIENTATION_SIZE),
            TransitionDef(coord_name="eo", coord_size=EDGE_ORIENTATION_SIZE),
            TransitionDef(coord_name="pcp", coord_size=PARTIAL_CORNER_PERMUTATION_SIZE) if cls.partial_corner_perm else
            TransitionDef(coord_name="cp", coord_size=CORNER_PERMUTATION_SIZE),
            TransitionDef(coord_name="pep", coord_size=PARTIAL_EDGE_PERMUTATION_SIZE) if cls.partial_edge_perm else
            TransitionDef(coord_name="ep", coord_size=EDGE_PERMUTATION_SIZE)]

        if not hasattr(cls, "pruning_defs"):
            cls.pruning_defs = [[] for _ in range(cls.num_phases)]

    def __init__(self, use_transition_tables: bool = True, use_pruning_tables: bool = True):
        """
        Create :class:`.BaseSolver` object.

        Parameters
        ----------
        use_transition_tables : bool, optional
            Whether to use transition tables for cube state transitions.
            If ``True``, creates or loads the tables from the ``tables/`` directory.
            Default is ``True``.
        use_pruning_tables : bool, optional
            Whether to use pruning tables to reduce the tree search space.
            If ``True``, creates or loads the tables from the ``tables/`` directory.
            Default is ``True``.

        See Also
        --------
        .solve
            Solve a cube position.
        """
        if not isinstance(use_transition_tables, bool):
            raise TypeError(f"use_transition_tables must be bool, not {type(use_transition_tables).__name__}")
        if not isinstance(use_pruning_tables, bool):
            raise TypeError(f"use_pruning_tables must be bool, not {type(use_pruning_tables).__name__}")

        init_coords = utils.flatten(self.get_coords(Cube()))
        self.solved_coords: List[FlattenCoords] = [self.phase_coords(init_coords, phase) for phase in range(self.num_phases)]
        """Solved flatten coordinates for each phase."""

        self.next_moves: List[Dict[Move, List[Move]]] = []
        """Allowed next moves based on the previous move, for each phase."""
        for phase_moves in self.phase_moves:
            next_moves = {Move.NONE: phase_moves}
            next_moves.update({move: [mv for mv in NEXT_MOVES[move] if mv in phase_moves] for move in phase_moves})
            self.next_moves.append(next_moves)

        self.final_moves: List[Set[Move]] = []
        """Final allowed moves for each phase, except the last."""
        for i in range(self.num_phases - 1):
            self.final_moves.append({Move.NONE} | set(self.phase_moves[i]) - set(self.phase_moves[i+1]))

        self.use_transition_tables: bool = use_transition_tables
        """Whether to use transition tables for cube state transitions."""
        self.use_pruning_tables: bool = use_pruning_tables
        """Whether to use pruning tables to reduce the tree search space."""
        self.transition_tables: Dict[str, np.ndarray] = {}
        """Transition tables used to compute cube state transitions."""
        if self.use_transition_tables:
            try:
                from .. import csolver  # type: ignore
                self.transition_tables = utils.get_tables("transition.npz", self.transition_defs,
                                                          csolver.generate_transition_table, accumulate=True)
            except Exception:  # pragma: no cover
                self.transition_tables = utils.get_tables("transition.npz", self.transition_defs,
                                                          utils.generate_transition_table, accumulate=True)
        self.pruning_tables: Dict[str, np.ndarray] = {}
        """Pruning tables used to reduce the tree search space."""
        if self.use_pruning_tables:
            pruning_defs = []
            for phase, phase_kwargs in enumerate(self.pruning_defs):
                for kwargs in phase_kwargs:
                    kwargs.solver = self
                    kwargs.phase = phase
                    pruning_defs.append(kwargs)
            if pruning_defs:
                pruning_filename = f"pruning_{self.__class__.__name__.lower()}.npz"
                try:
                    from .. import csolver  # type: ignore
                    self.pruning_tables = utils.get_tables(pruning_filename, pruning_defs, csolver.generate_pruning_table)
                except Exception:
                    self.pruning_tables = utils.get_tables(pruning_filename, pruning_defs, utils.generate_pruning_table)

        self.nodes: List[int] = [0] * self.num_phases
        """Number of visited nodes during a solve for each phase."""
        self.checks: List[int] = [0] * self.num_phases
        """Number of solve checks during a solve for each phase."""
        self.prunes: List[int] = [0] * self.num_phases
        """Number of pruned nodes during a solve for each phase."""
        self.terminated: bool = False
        """Whether the solve search was terminated due to a timeout."""

    def __repr__(self) -> str:
        """Solver string representation."""
        return self.__class__.__name__

    @staticmethod
    @abstractmethod
    def phase_coords(coords: FlattenCoords, phase: int) -> FlattenCoords:
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

        Notes
        -----
        Depending on the class attributes :attr:`partial_corner_perm` and :attr:`partial_edge_perm`,
        the ``coords`` parameter is the flattened version of the output from the :meth:`.get_coords` method.
        """

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

        Examples
        --------
        >>> from cube_solver import Cube, Kociemba
        >>> solver = Kociemba()
        >>> cube = Cube("U F2 R2")
        >>> solver.get_coords(cube)
        (0, 0, 26939, (1007, 11859, 673))
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

        Examples
        --------
        >>> from cube_solver import Cube, Kociemba
        >>> solver = Kociemba()
        >>> cube = Cube(random_state=True)
        >>> coords = (0, 0, 26939, (1007, 11859, 673))
        >>> solver.set_coords(cube, coords)
        >>> solver.get_coords(cube)
        (0, 0, 26939, (1007, 11859, 673))
        """
        cube.set_coords(coords, self.partial_corner_perm, self.partial_edge_perm)

    def is_solved(self, position: Union[Cube, CoordsType], phase: int) -> bool:
        """
        Whether the cube position is solved at the specified phase.

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

        Examples
        --------
        >>> from cube_solver import Cube, Kociemba
        >>> solver = Kociemba()
        >>> cube = Cube("U F2 R2")
        >>> solver.is_solved(cube, phase=0)
        True
        >>> solver.is_solved(cube, phase=1)
        False
        """
        if isinstance(position, Cube):
            position = self.get_coords(position)
        phase_coords = self.phase_coords(utils.flatten(position), phase)
        return phase_coords == self.solved_coords[phase]

    def prune(self, position: Union[Cube, CoordsType], phase: int, depth: int) -> bool:
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

        Examples
        --------
        >>> from cube_solver import Cube, Kociemba
        >>> solver = Kociemba()
        >>> cube = Cube("U F2 R2")
        >>> solver.prune(cube, phase=1, depth=2)
        True
        >>> solver.prune(cube, phase=1, depth=3)
        False
        """
        if self.pruning_tables:
            if isinstance(position, Cube):
                position = self.get_coords(position)
            phase_coords = self.phase_coords(utils.flatten(position), phase)
            for kwargs in self.pruning_defs[phase]:
                prune_coords = utils.select(phase_coords, kwargs.indexes)
                if self.pruning_tables[kwargs.name][prune_coords] > depth:
                    return True
        return False

    @overload
    def next_position(self, position: Cube, move: Move) -> Cube: ...
    @overload
    def next_position(self, position: CoordsType, move: Move) -> CoordsType: ...

    def next_position(self, position: Union[Cube, CoordsType], move: Move) -> Union[Cube, CoordsType]:
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

        Examples
        --------
        >>> from cube_solver import Cube, Move, Kociemba
        >>> solver = Kociemba()
        >>> cube = Cube("R2")
        >>> coords = solver.get_coords(cube)
        >>> solver.next_position(cube, Move.R2)
        WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY
        >>> solver.next_position(coords, Move.R2)
        (0, 0, 0, (0, 11856, 1656))
        """
        if isinstance(position, Cube):
            return apply_move(position, move)
        if self.use_transition_tables:
            next_position = ()
            for coord, kwargs in zip(position, self.transition_defs):
                next_coord = self.transition_tables[kwargs.name][coord, MOVE_TO_INDEX[move]]
                next_position += (next_coord.item(),) if isinstance(coord, int) else (tuple(next_coord.tolist()),)
            return next_position
        cube = Cube()
        self.set_coords(cube, position)
        cube.apply_move(move)
        return self.get_coords(cube)

    def solve(self, cube: Cube, max_length: Union[int,  None] = None, optimal: bool = False,
              timeout: Union[int, None] = None, verbose: int = 0) -> Union[Maneuver, List[Maneuver], None]:
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
        timeout : int or None, optional
            Stop the search after the specified number of seconds.
            If ``None``, no time limit is applied. Default is ``None``.
        verbose : {0, 1, 2}, optional
            Verbosity level. Default is ``0``.

            * ``0`` returns only the solution.
            * ``1`` logs all solutions found if ``optimal`` is ``True``.
            * ``2`` returns the solution for each phase.

        Returns
        -------
        solution : Maneuver or list of Maneuver or None
            Solution for the cube position,
            or list of solutions for each phase if ``verbose`` is ``2``,
            or ``None`` if no solution is found.

        Examples
        --------
        >>> from cube_solver import Cube, Kociemba
        >>> solver = Kociemba()
        >>> cube = Cube("L2 U R D' B2 D2 F B D")
        >>> solver.solve(cube)
        "D' F' B' U2 F2 D L' F2 D2 L2 F2 U D L2 B2 D L2"
        >>> solver.solve(cube, optimal=True, verbose=2)
        ["D' F' B' D2 B2 D R'", "U' L2"]
        """
        if not isinstance(cube, Cube):
            raise TypeError(f"cube must be Cube, not {type(cube).__name__}")
        if max_length is not None and not isinstance(max_length, int):
            raise TypeError(f"max_length must be int or None, not {type(max_length).__name__}")
        if not isinstance(optimal, bool):
            raise TypeError(f"optimal must be bool, not {type(optimal).__name__}")
        if timeout is not None and not isinstance(timeout, int):
            raise TypeError(f"timeout must be int or None, not {type(timeout).__name__}")
        if not isinstance(verbose, int):
            raise TypeError(f"verbose must be int, not {type(verbose).__name__}")
        if isinstance(max_length, int) and max_length < 0:
            raise ValueError(f"max_length must be >= 0 (got {max_length})")
        if isinstance(timeout, int) and timeout < 0:
            raise ValueError(f"timeout must be >= 0 (got {timeout})")
        if verbose not in (0, 1, 2):
            raise ValueError(f"verbose must be one of 0, 1, 2 (got {verbose})")
        try:
            cube.coords
        except ValueError:
            raise ValueError("invalid cube state")
        cube = Cube(repr=(repr(cube)))
        if cube.permutation_parity is None:
            raise ValueError("invalid cube state")

        for phase in range(self.num_phases):
            self.nodes[phase] = 0
            self.checks[phase] = 0
            self.prunes[phase] = 0
        self.terminated = False

        self._max_length = max_length
        self._optimal = optimal
        self._timeout = timeout
        self._verbose = verbose

        self._return_phase = False
        self._start = time.time()
        self._best_solution = None
        self._solution = [[] for _ in range(self.num_phases)]

        solution = None
        position = self.get_coords(cube) if self.use_transition_tables else cube
        if self._phase_search(position):
            solution = self._solution
        if self._optimal:
            solution = self._best_solution
        if solution:
            if verbose == 2:
                return [Maneuver(phase_solution[-2::-1]) for phase_solution in solution]
            return Maneuver([move for phase_solution in solution for move in phase_solution[-2::-1]])
        return None

    def _phase_search(self, position: Union[Cube, CoordsType], phase: int = 0, current_length: int = 0) -> bool:
        """
        Solve the cube position from the specified phase.

        Parameters
        ----------
        position : Cube or tuple of (int or tuple of int)
            Cube object or cube coordinates to search.
        phase : int, optional
            Phase to solve from (0-indexed). Default is ``0``.
        current_length : int, optional
            Current solution length. Default is ``0``.

        Returns
        -------
        bool
            ``True`` if a solution is found, ``False`` otherwise.
        """
        if phase == self.num_phases:
            if self._optimal:
                self._return_phase = True
                self._max_length = current_length - 1
                self._best_solution = deepcopy(self._solution)
                if self._verbose == 1:
                    solution = Maneuver([move for sol in self._solution for move in sol[-2::-1]], reduce=False)
                    logger.info(f"Solution: {solution} ({len(solution)})")
                elif self._verbose == 2:
                    solution = [Maneuver(phase_solution[-2::-1]) for phase_solution in self._solution]
                    logger.info(f"Solution: {' | '.join([f'{sol} ({len(sol)})' for sol in solution])}")
                return False
            return True
        depth = 0
        while True if self._max_length is None else current_length + depth <= self._max_length:
            self._solution[phase].append(Move.NONE)
            if self._search(position, phase, depth, current_length):
                return True
            elif self.terminated or (self._optimal and self._return_phase):
                self._return_phase = False
                self._solution[phase] = []
                return False
            depth += 1
        self._solution[phase] = []
        return False

    def _search(self, position: Union[Cube, CoordsType], phase: int, depth: int, current_length: int) -> bool:
        """
        Solve the cube position from the specified phase at the specified depth.

        Parameters
        ----------
        position : Cube or tuple of (int or tuple of int)
            Cube object or cube coordinates to search.
        phase : int
            Phase to solve from (0-indexed).
        depth : int
            Current search depth.
        current_length : int
            Current solution length.

        Returns
        -------
        bool
            ``True`` if a solution is found, ``False`` otherwise.
        """
        self.nodes[phase] += 1
        if self._timeout is not None and int(time.time() - self._start) >= self._timeout:
            self.terminated = True
            return False
        if depth == 0:
            if phase == self.num_phases - 1 or (self._solution[phase][0] in self.final_moves[phase]):
                self.checks[phase] += 1
                if self.is_solved(position, phase):
                    return self._phase_search(position, phase + 1, current_length)
            return False
        if not self.prune(position, phase, depth):
            for move in self.next_moves[phase][self._solution[phase][depth]]:
                self._solution[phase][depth - 1] = move
                if self._search(self.next_position(position, move), phase, depth - 1, current_length + 1):
                    return True
                elif self.terminated or (self._optimal and self._return_phase):
                    return False
            return False
        self.prunes[phase] += 1
        return False
