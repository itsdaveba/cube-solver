import os
from cube_solver import Cube, Solver


def test_solver():
    tables = ["pco", "peo", "pcp", "pep", "tco", "teo", "tcp", "tep"]
    for table in tables:
        if os.path.isfile(f"tables/{table}.pkl"):
            os.remove(f"tables/{table}.pkl")
    if os.path.isdir("tables/"):
        os.removedirs("tables/")

    cube = Cube("F B R' L'", representation="face")
    solver = Solver()
    assert solver.solve(cube) == "R L F' B'"

    cube = Cube("U' F B2 D", representation="array")
    solver = Solver()
    assert solver.solve(cube) == "D' F' B2 U"

    cube = Cube("B2 R L2 B2", representation="cubie")
    solver = Solver()
    assert solver.solve(cube) == "B2 R' L2 B2"

    cube = Cube("D2 U' L D'", representation="coord")
    solver = Solver()
    assert solver.solve(cube) == "D L' U D2"

    cube = Cube("F L2 D U' L F", representation="coord")
    solver = Solver(transition_tables=True)
    assert solver.solve(cube) == "F' L' U D' L2 F'"

    cube = Cube("L2 B D' F2 L' F' R", representation="coord")
    solver = Solver(pruning_tables=True)
    assert solver.solve(cube) == "R' F L F2 D B' L2"

    cube = Cube("D F' R' F B' D' F2 B2 L2 U2", representation="coord")
    solver = Solver(pruning_tables=True, transition_tables=True)
    assert solver.solve(cube) == "U2 L2 F2 B2 D F' B R F D'"
