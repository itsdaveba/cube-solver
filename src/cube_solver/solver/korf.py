from cube_solver import Cube, generate_scramble
from cube_solver.solver import BaseSolver


class Korf(BaseSolver):
    solved_coords = [(0, 0, 0, 1656, 11856, 1656, 0)]
    pruning_names = [["korf_co", "korf_eo", "korf_pcp1", "korf_pcp2", "korf_pep1", "korf_pep2", "korf_pep3"]]
    pruning_kwargs = [[dict(shape=(2187,), indexes=[0]),
                       dict(shape=(2048,), indexes=[1]),
                       dict(shape=(1680,), indexes=[2]),
                       dict(shape=(1680,), indexes=[3]),
                       dict(shape=(11880,), indexes=[4]),
                       dict(shape=(11880,), indexes=[5]),
                       dict(shape=(11880,), indexes=[6])]]

    def _phase_coords(self, phase: int, coords: tuple) -> tuple:
        return coords[:2] + coords[2] + coords[3]


if __name__ == "__main__":
    scramble = generate_scramble(9)
    print("Scramble:", scramble)
    cube = Cube(scramble)
    solver = Korf(use_transition_tables=True)
    solution = solver.solve(cube, verbose=1)
    print("Solution:", solution)
