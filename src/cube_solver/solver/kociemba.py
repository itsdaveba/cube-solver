import time
from cube_solver import Cube, Move, generate_scramble
from cube_solver.solver import BaseSolver

PHASE1_MOVES = list(reversed(Move))
RESTRICT = [Move.F1, Move.F3, Move.B1, Move.B3, Move.R1, Move.R3, Move.L1, Move.L3]
PHASE2_MOVES = [move for move in list(reversed(Move)) if move not in RESTRICT]


class Kociemba(BaseSolver):
    num_phases = 2
    partial_corner_perm = False
    phase_moves = [PHASE1_MOVES, PHASE2_MOVES]
    solved_coords = [(0, 0, 494), (0, 0, 0)]
    pruning_names = [["kociemba_coeoec"], ["kociemba_cpsep", "kociemba_epsep", "kociemba_cpep"]]
    pruning_kwargs = [[dict(shape=(2187, 2048, 495), indexes=[0, 1, 2])],
                      [dict(shape=(40320, 24), indexes=[0, 2]),
                       dict(shape=(40320, 24), indexes=[1, 2]),
                       dict(shape=(40320, 40320), indexes=[0, 1])]]

    def _phase_coords(self, phase: int, coords: tuple) -> tuple:
        if phase == 0:
            corner_orientation = coords[0]
            edge_orientation = coords[1]
            edge_combination = coords[3][0] // 24
            return (corner_orientation, edge_orientation, edge_combination)
        if phase == 1:
            corner_permutation = coords[2]
            edge_permutation = coords[3][1] + (coords[3][2] + coords[3][2] // 24 - 69) * 24  # TODO double-check
            slice_edge_permutation = coords[3][0] % 24
            return (corner_permutation, edge_permutation, slice_edge_permutation)
        raise ValueError


if __name__ == "__main__":
    cube = Cube(random_state=True)
    solver = Kociemba(use_transition_tables=True)
    solution = solver.solve(cube, verbose=2)
    print(solution)

    # num_solves = 1000
    # max_moves = 20  # God's Number

    # avg_time = 0.0
    # cube = Cube()
    # solver = Kociemba(use_transition_tables=True)

    # for i in range(num_solves):
    #     print(f"Solve Number: {i + 1}")
    #     cube.set_random_state()
    #     print("Cube:", repr(cube))

    #     start = time.time()
    #     solution = solver.solve(cube, verbose=1, optimal=max_moves)
    #     end = time.time()
    #     avg_time += end - start

    #     if solution is not None:
    #         print("Solution:", solution)
    #     else:
    #         print("No Solution Found")

    #     print(f"Time: {end - start:.2f} seconds")
    #     print(f"Average Time: {avg_time / (i + 1):.2f} seconds")
    #     print()
