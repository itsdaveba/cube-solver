# from cube_solver import Cube, generate_scramble
# from cube_solver.solver import BaseSolver


# class DummySolver(BaseSolver):
#     def _phase_coords(self, phase: int, coords: tuple) -> tuple:
#         return coords


# if __name__ == "__main__":
#     scramble = generate_scramble(5)
#     print("Scramble:", scramble)
#     cube = Cube(scramble)
#     solver = DummySolver(use_transition_tables=True)
#     solution = solver.solve(cube, verbose=1)
#     print("Solution:", solution)
