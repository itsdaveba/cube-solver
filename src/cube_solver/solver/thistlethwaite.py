# import time
# from cube_solver import Cube, Move, generate_scramble
# from cube_solver.solver import BaseSolver
# from cube_solver.constants import CORNER_THREAD

# PHASE1_MOVES = list(reversed(Move))
# RESTRICT = [Move.F1, Move.F3, Move.B1, Move.B3]
# PHASE2_MOVES = [move for move in reversed(Move) if move not in RESTRICT]
# RESTRICT += [Move.R1, Move.R3, Move.L1, Move.L3]
# PHASE3_MOVES = [move for move in reversed(Move) if move not in RESTRICT]
# RESTRICT += [Move.U1, Move.U3, Move.D1, Move.D3]
# PHASE4_MOVES = [move for move in reversed(Move) if move not in RESTRICT]


# class Thistlethwaite(BaseSolver):
#     num_phases = 4
#     phase_moves = [PHASE1_MOVES, PHASE2_MOVES, PHASE3_MOVES, PHASE4_MOVES]
#     solved_coords = [(0,), (0, 494), (0, 0, 0), (0, 0, 0, 0, 0)]
#     pruning_names = [["thistle_phase1"], ["thistle_phase2"], ["thistle_phase3"], ["thistle_phase4"]]
#     pruning_kwargs = [[dict(shape=(2048,), indexes=[0])],
#                       [dict(shape=(2187, 495), indexes=[0, 1])],
#                       [dict(shape=(70, 70, 6), indexes=[0, 1, 2])],
#                       [dict(shape=(24, 4, 24, 24, 12), indexes=[0, 1, 2, 3, 4])]]

#     def _phase_coords(self, phase: int, coords: tuple) -> tuple:
#         if phase == 0:
#             edge_orientation = coords[1]
#             return (edge_orientation,)
#         if phase == 1:
#             corner_orientation = coords[0]
#             edge_combination = coords[3][0] // 24
#             return (corner_orientation, edge_combination)
#         if phase == 2:
#             corner_combination = coords[2][0] // 24
#             edgeud_combination = coords[3][2] // 24
#             corner_thread = CORNER_THREAD[coords[2][0] % 24, coords[2][1] % 24].item()
#             return (corner_combination, edgeud_combination, corner_thread)
#         if phase == 3:
#             corner_permutation = [coord % 24 for coord in coords[2]]
#             corner_permutation[-1] //= 6
#             edge_permutation = [coord % 24 for coord in coords[3]]
#             edge_permutation[-1] //= 2
#             return tuple(corner_permutation) + tuple(edge_permutation)
#         raise ValueError


# if __name__ == "__main__":

#     cube = Cube(random_state=True)
#     solver = Thistlethwaite(use_transition_tables=True)
#     solution = solver.solve(cube, verbose=2)
#     print(solution)

#     # num_solves = 1000
#     # max_moves = 20  # God's Number

#     # avg_time = 0.0
#     # cube = Cube()
#     # solver = Thistlethwaite(use_transition_tables=True)

#     # for i in range(num_solves):
#     #     print(f"Solve Number: {i + 1}")
#     #     cube.set_random_state()
#     #     print("Cube:", repr(cube))

#     #     start = time.time()
#     #     solution = solver.solve(cube, verbose=1, optimal=max_moves)
#     #     end = time.time()
#     #     avg_time += end - start

#     #     if solution is not None:
#     #         print("Solution:", solution)
#     #     else:
#     #         print("No Solution Found")

#     #     print(f"Time: {end - start:.2f} seconds")
#     #     print(f"Average Time: {avg_time / (i + 1):.2f} seconds")
#     #     print()
