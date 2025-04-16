from cube_solver import Solver


def test_solver():
    solver = Solver()

    solver.scramble("R L2 F'")
    assert solver.solve() == "F R' L2"

    solver.scramble("L F' L")
    assert solver.solve() == "L' F L'"

    solver.scramble("D L2 U'")
    assert solver.solve() == "U L2 D'"
