# import numpy as np
# from itertools import count
# from collections import deque
# from typing import overload
# from abc import ABC, abstractmethod

# import cube_solver
# from cube_solver import Cube, Move
# from cube_solver.cube.enums import Face
# from cube_solver.solver import utils
# from cube_solver.constants import EMPTY, SOLVED_PARTIAL_COORDS

# # TODO check types for each method

# OPPOSITE = {Face[face]: Face[opp] for face, opp in zip("UFRDBL", "DBLUFR")}
# FACES = list(Face)
# UFR = [Face.U, Face.F, Face.R]
# NEXT_FACES = {face: [f for f in FACES if f != face and (f != OPPOSITE[face] or face in UFR)] for face in FACES}
# NEXT_FACES[Face.NONE] = FACES
# SOLVED_COORDS = (0, 0, 0, 0)


# class BaseSolver(ABC):
#     num_phases: int = 1
#     phase_moves: list = [list(Move)]
#     solved_coords: list = [SOLVED_PARTIAL_COORDS]
#     pruning_names: list = [[]]
#     pruning_kwargs: list = [[]]
#     partial_corner_perm: bool = True
#     partial_edge_perm: bool = True
#     transition_names: list
#     transition_kwargs: list

#     def __init_subclass__(cls) -> None:
#         # for attr in ("num_phases", "next_moves", "solved_coords", "prun_tables_names", "prun_tables_kwargs"):
#         #     if not hasattr(cls, attr):
#         #         raise AttributeError(f"{cls.__name__} class must define class attribute '{attr}'")
#         cls.transition_names = [
#             "transition_eo",
#             "transition_co",
#             "transition_pcp" if cls.partial_corner_perm else "transition_cp",
#             "transition_pep" if cls.partial_edge_perm else "transition_cp"
#         ]
#         cls.transition_kwargs = [
#             dict(coord_type="co", size=2187),
#             dict(coord_type="eo", size=2048),
#             dict(coord_type="pcp", size=1680) if cls.partial_corner_perm else dict(coord_type="cp", size=40320),
#             dict(coord_type="pep", size=11880) if cls.partial_edge_perm else dict(coord_type="ep", size=479001600)
#         ]

#     def __init__(self, use_transition_tables: bool = False):
#         self._cube = Cube()
#         self.use_transition_tables = use_transition_tables

#         self.next_moves = []
#         for phase_moves in self.phase_moves:
#             next_moves = {face: [move for move in phase_moves if move.face in NEXT_FACES[face]] for face in Face}
#             next_moves[Face.NONE] = phase_moves
#             self.next_moves.append(next_moves)

#         self.transition_tables = []
#         if self.use_transition_tables:
#             self.transition_tables = self.get_tables(self.transition_names, self.transition_kwargs, self._generate_trans_table)
#         self.pruning_tables = self.get_prun_tables()

#     def get_tables(self, tables_names, tables_kwargs, generate_table_fn) -> list:
#         tables = []
#         for idx, name in enumerate(tables_names):
#             try:
#                 table = utils.load_table(name)
#             except FileNotFoundError:
#                 table = generate_table_fn(**tables_kwargs[idx])
#                 utils.save_table(name, table)
#             tables.append(table)
#         return tables

#     def get_prun_tables(self) -> list:
#         tables = []
#         for phase, (tables_names, tables_kwargs) in enumerate(zip(self.pruning_names, self.pruning_kwargs)):
#             for kwargs in tables_kwargs:
#                 kwargs["phase"] = phase
#             phase_tables = self.get_tables(tables_names, tables_kwargs, self._generate_prun_table)
#             tables.append(phase_tables)
#         return tables

#     def _generate_trans_table(self, coord_type: str, size: int) -> np.ndarray:  # TODO maybe move to utils
#         assert size < 65536  # TODO change to ValError
#         transition_table = np.zeros((size, len(Move)), dtype=np.uint16)
#         for coord in range(size):
#             self._cube.set_coord(coord_type, coord)
#             for move in Move:
#                 transition_table[coord, move] = cube_solver.apply_move(self._cube, move).get_coord(coord_type)
#         return transition_table

#     def _generate_prun_table(self, phase: int, shape: tuple, indexes: list) -> np.ndarray:
#         self._cube.reset()
#         coords = self._cube.get_coords(self.partial_corner_perm, self.partial_edge_perm)
#         phase_coords = self._phase_coords(phase, coords)
#         prune_coords = tuple(phase_coords[i] for i in indexes)
#         pruning_table = np.full(shape, EMPTY, dtype=np.int8)
#         pruning_table[prune_coords] = 0
#         queue = deque([(coords, 0)])
#         while queue:
#             coords, depth = queue.popleft()
#             for move in self.next_moves[phase][Face.NONE]:
#                 next_coords = self._next_position(coords, move)
#                 phase_coords = self._phase_coords(phase, next_coords)
#                 prune_coords = tuple(phase_coords[i] for i in indexes)
#                 if pruning_table[prune_coords] == EMPTY:
#                     pruning_table[prune_coords] = depth + 1
#                     queue.append((next_coords, depth + 1))
#         return pruning_table

#     @abstractmethod
#     def _phase_coords(self, phase: int, coords: tuple) -> tuple: ...

#     def _set_coords(self, phase: int, cube: Cube, coords: tuple):
#         cube.coords = coords

#     @overload
#     def _next_position(self, position: Cube, move: Move) -> Cube: ...
#     @overload
#     def _next_position(self, position: tuple, move: Move) -> tuple: ...

#     def _next_position(self, position: Cube | tuple, move: Move) -> Cube | tuple:
#         if isinstance(position, Cube):
#             return cube_solver.apply_move(position, move)
#         elif isinstance(position, tuple):
#             if self.use_transition_tables:
#                 return tuple(
#                     tuple(table[coord, move].tolist()) if isinstance(coord, tuple) else table[coord, move].item()
#                     for coord, table in zip(position, self.transition_tables)
#                 )
#             self._cube.set_coords(position, self.partial_corner_perm, self.partial_edge_perm)
#             self._cube.apply_move(move)
#             return self._cube.get_coords(self.partial_corner_perm, self.partial_edge_perm)
#         raise ValueError

#     def is_solved(self, position: Cube | tuple, phase: int | None = None) -> bool:
#         """
#         Check whether the `cube` position is solved at the current `phase`.

#         Parameters
#         ----------
#         cube : Cube
#             Cube object to check.
#         phase : int or None, optional
#             Phase to check (0-indexed). If `None`, checks whether the cube is fully solved.

#         Returns
#         -------
#         bool
#             `True` if the cube position is solved, `False` otherwise.

#         Examples
#         --------
#         >>> from cube_solver import Cube, Kociemba
#         >>> cube = Cube("U F2 R2")
#         >>> solver = Kociemba()
#         >>> solver.is_solved(cube)
#         False
#         >>> solver.is_solved(cube, phase=0)
#         True
#         """
#         if isinstance(position, Cube):
#             position = position.get_coords(self.partial_corner_perm, self.partial_edge_perm)
#         elif not isinstance(position, tuple):
#             raise ValueError

#         if phase is None:
#             self._cube.reset()
#             return position == self._cube.get_coords(self.partial_corner_perm, self.partial_edge_perm)
#         return self._phase_coords(phase, position) == self.solved_coords[phase]

#     def solve(self, cube: Cube, max_depth: int | None = None, verbose: int = 0, optimal=False) -> str | None:
#         """
#         Solve the `cube` position.

#         Parameters
#         ----------
#         cube : Cube
#             Cube object to be solved.
#         max_depth : int or None, optional
#             The maximum number of moves to search. If `None`, search indefinitely.
#         verbose : {0, 1, 2}, optional
#             Verbosity level. If `0`, only return the solution.
#             If `1`, return the solution with the move count.
#             If `2`, return the solution separated by phase and move count for each phase.

#         Returns
#         -------
#         solution : str or None
#             Solution for the cube position, or `None` if no solution is found.

#         Examples
#         --------
#         >>> from cube_solver import Cube, Kociemba
#         >>> cube = Cube("U F R")
#         >>> solver = Kociemba()
#         >>> solver.solve(cube)
#         R' F' U'
#         >>> solver.solve(cube, verbose=1)
#         R' F' U' (3)
#         >>> solver.solve(cube, verbose=2)
#         R' F' (2) | U' (1)
#         """
#         if not isinstance(cube, Cube):
#             raise TypeError(f"cube must be Cube, not {type(cube).__name__}")
#         if max_depth is not None and not isinstance(max_depth, int):
#             raise TypeError(f"max_depth must be int or None, not {type(max_depth).__name__}")
#         if not isinstance(verbose, int):
#             raise TypeError(f"verbose must be int, not {type(verbose).__name__}")

#         coords = cube.get_coords(self.partial_corner_perm, self.partial_edge_perm)
#         solution = self.search(coords, max_depth or 100, optimal=optimal, verbose=verbose)
#         if solution[-1][0][0] < 0:
#             return None

#         if verbose in (0, 1):
#             solution = [move.str for sol in solution[:-1] for _, move in sol[1:]]
#             if verbose == 0:
#                 return " ".join(solution)
#             return " ".join(solution) + f" ({len(solution)})"
#         if verbose == 2:
#             return " | ".join(" ".join(move.str for _, move in sol[1:]) + f" ({len(sol[1:])})" for sol in solution[:-1])
#         raise ValueError(f"verbose must be one of 0, 1, or 2 (got '{verbose}')")

#     def search(self, init_coords, max_depth, solution=[], phase=0, optimal=False, verbose=0):
#         if phase >= self.num_phases:
#             if optimal and verbose > 0:
#                 print(" | ".join(" ".join(sol) + f" ({len(sol)})" for sol in [[move[1].str for move in sol[1:]] for sol in solution]))
#             return [[(0, Move.NONE)]]

#         sol = [[(-1, Move.NONE)]]
#         phase_depth = 0
#         phase_solution = []
#         while phase_depth <= max_depth:
#             stack = [(phase_depth, init_coords, Move.NONE)]
#             while stack:
#                 depth, coords, last_move = stack.pop()
#                 while phase_solution and phase_solution[-1][0] <= depth:
#                     phase_solution.pop()
#                 phase_solution.append((depth, last_move))
#                 print(phase, phase_solution)
#                 if depth == 0:
#                     if self.is_solved(coords, phase):
#                         phase_length = phase_solution[0][0]
#                         next_phase_solution = self.search(coords, max_depth - phase_length, solution + [phase_solution], phase + 1, optimal, verbose)  # TODO restrict moves to advance e.g. F1, F3, R1, R3 for Kociemba
#                         if not optimal or phase > 0:
#                             return [phase_solution] + next_phase_solution
#                         if next_phase_solution[-1][0][0] >= 0:
#                             sol = [phase_solution[:]] + next_phase_solution
#                             next_phase_length = sum([s[0][0] for s in next_phase_solution])
#                             max_depth = phase_length + next_phase_length - 1
#                             if max_depth < optimal:
#                                 return sol
#                     continue
#                 if not self._prune(coords, depth, phase):
#                     for move in self.next_moves[phase][last_move.face]:
#                         next_coords = self._next_position(coords, move)
#                         stack.append((depth - 1, next_coords, move))
#             phase_depth += 1
#         return sol

#     def _solve(self, position: Cube | tuple, depth: int, phase: int, solution: deque, last_face: Face = Face.NONE) -> bool:
#         """
#         Solve the `cube` position at the current `phase` and `depth`.

#         Parameters
#         ----------
#         cube : Cube
#             Cube object to be solved.
#         depth : int
#             Search depth.
#         phase : int, optional
#             Phase to solve.
#         solution : deque
#             If a soluition is foud, stores the sequence of moves for the solution.
#         last_face : Face, optional
#             Last face move performed on the cube. This helps avoid repeating the same move during the search.
#             `Face.NONE` indicates that no face moves have been made yet (i.e. the first move).

#         Returns
#         -------
#         bool
#             `True` if a solution was found, `False` otherwise.
#         """
#         if depth == 0:
#             return self.is_solved(position, phase)
#         if not self._prune(position, depth, phase):
#             for move in self.next_moves[phase][last_face]:
#                 next_position = self._next_position(position, move)
#                 if self._solve(next_position, depth - 1, phase, solution, move.face):
#                     solution.appendleft(move)
#                     return True
#         return False

#     def _prune(self, position: Cube | tuple, depth: int, phase: int) -> bool:
#         """
#         Prune the search tree.

#         Check whether the current `depth` meets the lower bound specified by the `pruning tables`.

#         Parameters
#         ----------
#         cube : Cube
#             Cube object to check.
#         depth : int
#             Current search depth.
#         phase : int
#             Phase to prune, used to select the appropiate `pruning tables`.

#         Returns
#         -------
#         bool
#             `True` if the search tree should be pruned, `False` otherwise.
#         """
#         if self.pruning_tables[phase]:
#             if isinstance(position, Cube):
#                 position = position.get_coords(self.partial_corner_perm, self.partial_edge_perm)
#             phase_coords = self._phase_coords(phase, position)
#             for table, kwargs in zip(self.pruning_tables[phase], self.pruning_kwargs[phase]):
#                 prune_coords = tuple(phase_coords[i] for i in kwargs["indexes"])
#                 if table[prune_coords] > depth:
#                     return True
#         return False
