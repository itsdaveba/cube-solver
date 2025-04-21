from cube_solver import Cube, Solver


def test_solver():
    cube = Cube("R L2 F'", representation="face")
    assert Solver(cube).solve() == "F R' L2"

    cube = Cube("L F' L", representation="array")
    assert Solver(cube).solve() == "L' F L'"

    cube = Cube("D L2 U'", representation="cubie")
    assert Solver(cube).solve() == "U L2 D'"

    cube = Cube("B' L B2", representation="coord")
    assert Solver(cube).solve() == "B2 L' B"
