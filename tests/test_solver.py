import os
import time
import pytest
import numpy as np

from cube_solver import Cube, Maneuver
from cube_solver import DummySolver, Korf, Thistlethwaite, Kociemba
from cube_solver.solver import BaseSolver, utils
from cube_solver.solver.defs import FlattenCoords, PruningDef


def num_nodes(depth: int) -> int:
    if depth <= 2:
        return 104 * depth * depth - 87 * depth + 1
    return num_nodes(depth - 1) * 12 + num_nodes(depth - 2) * 18


def solve_test(solver: BaseSolver, num_cubes: int):
    print()
    print("Solve Test")
    print("----------")
    print(f"Solver: {solver.__class__.__name__}")
    print(f"Transition Tables: {solver.use_transition_tables}")
    print(f"Pruning Tables: {solver.use_pruning_tables}")
    print(f"Num Cubes: {num_cubes}")
    times = []
    nodes = []
    for i in range(num_cubes):
        cube = Cube(random_state=True)
        times.append(time.perf_counter())
        solution = solver.solve(cube)
        times[i] = time.perf_counter() - times[i]
        nodes.append([sum(nds) for phase_nodes in solver.nodes for nds in phase_nodes])
        assert solution is not None
        print("#", end="")
    nodes = np.mean(nodes, axis=0)  # TODO add solve length
    print("\nMin Time\tMax Time\tMean Time\tStd Time\tMean Phase Nodes")
    print("--------------------------------------------------------------------------------")
    print(f"{np.min(times):.4f} sec\t{np.max(times):.4f} sec", end="\t")
    print(f"{np.mean(times):.4f} sec\t{np.std(times):.4f} sec", end="\t")
    for phase_nodes in nodes:
        print(f"{phase_nodes:.0f}", end=" ")
    print()


def depth_test(solver: BaseSolver, max_depth: int, num_cubes: int):  # TODO add stats about nodes and checks
    print()
    print("Depth Test")
    print("----------")
    print(f"Solver: {solver.__class__.__name__}")
    print(f"Transition Tables: {solver.use_transition_tables}")
    print(f"Pruning Tables: {solver.use_pruning_tables}")
    print(f"Max Depth: {max_depth}")
    print(f"Num Cubes: {num_cubes}")
    print("\nDepth\tSolving\t\tMin Time\tMax Time\tMean Time\tStd Time\tMean Phase Nodes")
    print("--------------------------------------------------------------------------------------------------------")
    for depth in range(max_depth + 1):
        times = []
        nodes = []
        print(f"{depth}", end="\t")
        for i in range(num_cubes):
            scramble = Maneuver.random(depth)
            cube = Cube(scramble)
            times.append(time.perf_counter())
            solution = solver.solve(cube)
            times[i] = time.perf_counter() - times[i]
            nodes.append([sum(nds) for phase_nodes in solver.nodes for nds in phase_nodes])
            assert solution == scramble.inverse
            print("#", end="")
        nodes = np.mean(nodes, axis=0)
        print(f"\t{np.min(times):.4f} sec\t{np.max(times):.4f} sec", end="\t")
        print(f"{np.mean(times):.4f} sec\t{np.std(times):.4f} sec", end="\t")
        for phase_nodes in nodes:
            print(f"{phase_nodes:.0f}", end=" ")
        print()


def test_solver():
    cube = Cube()
    cube.set_coord("pep", 0)
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
    with pytest.raises(TypeError, match=r"use_transition_tables must be bool, not NoneType"):
        TestSolver(None, None)
    with pytest.raises(TypeError, match=r"use_pruning_tables must be bool, not NoneType"):
        TestSolver(False, None)
    solver = TestSolver(False, False)
    assert solver.transition_tables == {}
    assert solver.pruning_tables == {}
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
    with pytest.raises(ValueError, match=r"verbose must be one of 0, 1, 2 \(got -1\)"):
        solver.solve(cube, 0, False, -1)
    with pytest.raises(ValueError, match=r"invalid cube state"):
        solver.solve(cube, 0, False, 0)
    with pytest.warns(UserWarning, match=r"invalid cube parity"):
        cube.set_coord("pep", (0, 11857, 1656))
    with pytest.warns(UserWarning, match=r"invalid cube parity"):
        with pytest.raises(ValueError, match=r"invalid cube state"):
            solver.solve(cube, 0, False, 0)
    with pytest.raises(ValueError, match=r"coord_size must be <= 65536 \(got 239500800\)"):
        solver.generate_transition_table("ep", 239500800)
    scramble = Maneuver("U F2 R'")
    cube = Cube(scramble)
    assert solver.solve(cube, 0, False, 0) is None
    assert solver.solve(cube, None, True, 2) == [scramble.inverse]
    if os.path.exists("tables/transition.npz"):
        tables = utils.load_tables("tables/transition.npz")
        if "cp" in tables:
            tables = {"cp": tables["cp"]}
            utils.save_tables("tables/transition.npz", tables)
        else:
            os.remove("tables/transition.npz")  # TODO just remove when creating tables with C++
    solver = TestSolver(True, False)
    assert solver.transition_tables != {}
    assert solver.pruning_tables == {}
    assert solver.solve(cube) == scramble.inverse
    class TestSolver(BaseSolver):
        partial_corner_perm = True
        partial_edge_perm = True
        pruning_kwargs = [[PruningDef(name="partial_corner_perm", shape=1680, indexes=2)]]
        def phase_coords(self, coords: FlattenCoords, phase: int) -> FlattenCoords:
            return coords
    solver = TestSolver(False, True)
    assert solver.transition_tables == {}
    assert solver.pruning_tables != {}
    assert solver.solve(cube) == scramble.inverse
    class TestSolver(BaseSolver):
        partial_corner_perm = True
        partial_edge_perm = True
        pruning_kwargs = [[PruningDef(name="pcp", shape=1680, indexes=2)]]
        def phase_coords(self, coords: FlattenCoords, phase: int) -> FlattenCoords:
            return coords
    solver = TestSolver(True, True)
    assert solver.transition_tables != {}
    assert solver.pruning_tables != {}
    assert solver.solve(cube) == scramble.inverse
    os.remove("tables/pruning_testsolver.npz")


def test_dummy():
    num_cubes = 12

    max_depth = 4
    solver = DummySolver()
    depth_test(solver, max_depth, num_cubes)


def test_korf():
    num_cubes = 12

    max_depth = 4
    solver = Korf(use_pruning_tables=False)
    depth_test(solver, max_depth, num_cubes)

    max_depth = 8
    solver = Korf()
    depth_test(solver, max_depth, num_cubes)


def test_thistlethwaite():
    num_cubes = 100
    solver = Thistlethwaite()
    with pytest.raises(ValueError, match=r"phase must be < 4 \(got 4\)"):
        solver.phase_coords((), 4)
    solve_test(solver, num_cubes)


def test_kociemba():
    num_cubes = 10
    solver = Kociemba()
    with pytest.raises(ValueError, match=r"phase must be < 2 \(got 2\)"):
        solver.phase_coords((), 2)
    solve_test(solver, num_cubes)
