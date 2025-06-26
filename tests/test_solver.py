import os
import time
import pytest
import numpy as np

from cube_solver import Cube, Maneuver
from cube_solver import DummySolver, Korf, Thistlethwaite, Kociemba
from cube_solver.solver import BaseSolver
from cube_solver.solver.defs import FlattenCoords, PruningDef


def num_checks(depth: int) -> int:
    if depth <= 2:
        return 104 * depth * depth - 87 * depth + 1
    return num_checks(depth - 1) * 12 + num_checks(depth - 2) * 18


def solve_cubes(solver: BaseSolver, scramble_length: int, num_cubes: int):
    print()
    print("Solve Test")
    print("----------")
    print(f"Solver: {solver.__class__.__name__}")
    print(f"Transition Tables: {solver.use_transition_tables}")
    print(f"Scramble Length: {scramble_length}")
    print(f"Num Cubes: {num_cubes}")
    print("Solving Cubes", end="")
    times = []
    for i in range(num_cubes):
        print(".", end="")
        scramble = Maneuver.random(scramble_length)
        cube = Cube(scramble)
        times.append(time.perf_counter())
        solution = solver.solve(cube)
        times[i] = time.perf_counter() - times[i]
        assert solution == scramble.inverse
    print()
    print(f"Min Time: {np.min(times)}")
    print(f"Max Time: {np.max(times)}")
    print(f"Avg Time: {np.mean(times)}")
    print(f"Std Time: {np.std(times)}")


def depth_test(solver: BaseSolver, max_depth: int, num_cubes: int):  # TODO add stats about nodes and checks
    print()
    print("Depth Test")
    print("----------")
    print(f"Solver: {solver.__class__.__name__}")
    print(f"Transition Tables: {solver.use_transition_tables}")
    print(f"Max Depth: {max_depth}")
    print(f"Num Cubes: {num_cubes}")
    print("\nDepth\tSolving\t\tMin\t\tMax\t\tMean\t\tStd")
    print("----------------------------------------------------------------------------------")
    for depth in range(max_depth + 1):
        times = []
        print(f"{depth}", end="\t")
        for i in range(num_cubes):
            scramble = Maneuver.random(depth)
            cube = Cube(scramble)
            times.append(time.perf_counter())
            solution = solver.solve(cube)
            times[i] = time.perf_counter() - times[i]
            assert solution == scramble.inverse
            print("#", end="")
        print(f"\t{np.min(times):.4f} sec\t{np.max(times):.4f} sec", end="\t")
        print(f"{np.mean(times):.4f} sec\t{np.std(times):.4f} sec")


def test_solver():
    cube = Cube()
    cube.set_coord("pcp", 0)
    with pytest.raises(AttributeError, match=r"'TestSolver' class must define class attribute 'partial_corner_perm'"):
        class TestSolver(BaseSolver):
            pass
    with pytest.raises(AttributeError, match=r"'TestSolver' class must define class attribute 'partial_edge_perm'"):
        class TestSolver(BaseSolver):
            partial_corner_perm = True
    class TestSolver(BaseSolver):
        partial_corner_perm = True
        partial_edge_perm = True
    with pytest.raises(TypeError):
        TestSolver(None)
    class TestSolver(BaseSolver):
        partial_corner_perm = True
        partial_edge_perm = True
        def phase_coords(self, coords: FlattenCoords, phase: int) -> FlattenCoords:
            return coords
    with pytest.raises(TypeError, match=r"use_pruning_tables must be bool, not NoneType"):
        TestSolver(None, None)
    with pytest.raises(TypeError, match=r"use_transition_tables must be bool, not NoneType"):
        TestSolver(True, None)
    with pytest.raises(ValueError, match=r"cannot use pruning tables with empty class attribute 'pruning_kwargs'"):
        TestSolver(True, False)
    solver = TestSolver(False, False)
    assert solver.pruning_tables == {}
    assert solver.transition_tables == {}
    with pytest.raises(TypeError, match=r"cube must be Cube, not NoneType"):
        solver.solve(None, "", None, None)
    with pytest.raises(TypeError, match=r"max_length must be int or None, not str"):
        solver.solve(cube, "", None, None)
    with pytest.raises(TypeError, match=r"optimal must be bool, not NoneType"):
        solver.solve(cube, -1, None, None)
    with pytest.raises(TypeError, match=r"verbose must be int, not NoneType"):
        solver.solve(cube, -1, False, None)
    with pytest.raises(ValueError, match=r"max_length must be >= 0 \(got -1\)"):
        solver.solve(cube, -1, False, -1)
    with pytest.raises(ValueError, match=r"verbose must be 0 or 1 \(got -1\)"):
        solver.solve(cube, 0, False, -1)
    with pytest.raises(ValueError, match=r"invalid cube state"):
        solver.solve(cube, 0, False, 0)
    with pytest.warns(UserWarning, match=r"invalid corner orientation"):
        cube = Cube(repr="OGWYYOBWWYGRGRRWGBYBGBGBRRYRYOOOBGOGGRBYBWYRBWWROWWOYO")
    with pytest.raises(ValueError, match=r"invalid cube state"):
        solver.solve(cube, 0, False, 0)
    scramble = Maneuver("U F2 R'")
    cube = Cube(scramble)
    assert solver.solve(cube, 0, False, 0) is None
    assert solver.solve(cube, None, True, 1) == scramble.inverse
    if os.path.exists("tables/transition.npz"):
        os.remove("tables/transition.npz")
    if os.path.exists("tables/pruning_testsolver.npz"):
        os.remove("tables/pruning_testsolver.npz")
    class TestSolver(BaseSolver):
        partial_corner_perm = True
        partial_edge_perm = True
        pruning_kwargs = [[PruningDef(name="pcp", shape=1680, indexes=2)]]
        def phase_coords(self, coords: FlattenCoords, phase: int) -> FlattenCoords:
            return coords
    with pytest.raises(ValueError, match=r"coord_size must be <= 65536 \(got 239500800\)"):
        solver.generate_transition_table("ep", 239500800)
    solver = TestSolver(True, False)
    assert solver.pruning_tables != {}
    assert solver.transition_tables == {}
    assert solver.solve(cube) == scramble.inverse
    solver = TestSolver(False, True)
    assert solver.pruning_tables == {}
    assert solver.transition_tables != {}
    assert solver.solve(cube) == scramble.inverse
    solver = TestSolver(True, True)
    assert solver.pruning_tables != {}
    assert solver.transition_tables != {}
    assert solver.solve(cube) == scramble.inverse
    os.remove("tables/pruning_testsolver.npz")


def test_dummy():
    num_cubes = 12

    max_depth = 3
    solver = DummySolver()
    depth_test(solver, max_depth, num_cubes)
    # assert solver.checks == [num_checks(i) for i in range(max_depth + 1)]
    # assert all(np.cumsum(solver.checks) == solver.nodes)

    max_depth = 4
    solver = DummySolver(use_transition_tables=True)
    depth_test(solver, max_depth, num_cubes)
    # assert solver.checks == [num_checks(i) for i in range(max_depth + 1)]
    # assert all(np.cumsum(solver.checks) == solver.nodes)


def test_korf():
    num_cubes = 12

    max_depth = 8
    solver = Korf()
    depth_test(solver, max_depth, num_cubes)

    max_depth = 10
    solver = Korf(use_transition_tables=True)
    depth_test(solver, max_depth, num_cubes)


def test_thistlethwaite():
    num_cubes = 12

    max_depth = 30
    solver = Thistlethwaite()
    depth_test(solver, max_depth, num_cubes)

    max_depth = 30
    solver = Thistlethwaite(use_transition_tables=True)
    depth_test(solver, max_depth, num_cubes)


def test_kociemba():
    num_cubes = 12

    # max_depth = 30
    # solver = Kociemba()
    # depth_test(solver, max_depth, num_cubes)

    # max_depth = 30
    # solver = Kociemba(use_transition_tables=True)
    # depth_test(solver, max_depth, num_cubes)
