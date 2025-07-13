import os
import time
import pytest
import numpy as np
from typing import Tuple

from cube_solver import Cube, Maneuver, apply_maneuver
from cube_solver import BaseSolver, DummySolver, Korf, Kociemba
from cube_solver.solver.defs import TransitionDef, PruningDef
from cube_solver.solver import utils


def solve_test(solver: BaseSolver, num_cubes: int):
    print()
    print("\nSolve Test")
    print("----------")
    print(f"Solver: {solver.__class__.__name__}")
    print(f"Transition Tables: {solver.use_transition_tables}")
    print(f"Pruning Tables: {solver.use_pruning_tables}")
    print(f"Num Cubes: {num_cubes}\n")
    times = []
    nodes = []
    lengths = []
    for i in range(num_cubes):
        cube = Cube(random_state=True)
        times.append(time.perf_counter())
        solution = solver.solve(cube)
        times[i] = time.perf_counter() - times[i]
        nodes.append(solver.nodes)
        assert isinstance(solution, Maneuver)
        assert apply_maneuver(cube, solution).is_solved
        lengths.append(len(solution))
        if i and i % 50 == 0:
            print(f" {i}")
        print("#", end="")
    print(f" {num_cubes}")
    nodes = np.mean(nodes, axis=0)
    print("\nMean Time\tMean Solution Length\tMean Phase Nodes")
    print("--------------------------------------------------------")
    print(f"{np.mean(times):.4f} sec", end="\t")
    print(f"{np.mean(lengths)} moves", end="\t\t")
    for phase_nodes in nodes:
        print(f"{phase_nodes:.0f}", end=" ")
    print()


def depth_test(solver: BaseSolver, max_depth: int, num_cubes: int):
    print()
    print("\nDepth Test")
    print("----------")
    print(f"Solver: {solver.__class__.__name__}")
    print(f"Transition Tables: {solver.use_transition_tables}")
    print(f"Pruning Tables: {solver.use_pruning_tables}")
    print(f"Max Depth: {max_depth}")
    print(f"Num Cubes: {num_cubes}")
    print("\nDepth\tSolving\t\tMean Time\tMean Phase Nodes")
    print("--------------------------------------------------------")
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
            nodes.append(solver.nodes)
            assert isinstance(solution, Maneuver)
            assert apply_maneuver(cube, solution).is_solved
            print("#", end="")
        nodes = np.mean(nodes, axis=0)
        print(f"\t{np.mean(times):.4f} sec", end="\t")
        for phase_nodes in nodes:
            print(f"{phase_nodes:.0f}", end=" ")
        print()


def test_utils():
    # select
    cube = Cube()
    cube.coords = (67, 2053)
    assert utils.select(cube.coords, None) == (67, 2053)
    assert utils.select(cube.coords, 0) == (67,)
    assert utils.select(cube.coords, (0, 1)) == (67, 2053)

    # load and save
    with pytest.raises(TypeError, match=r"path must be str or Path, not NoneType"):
        utils.save_tables(None, None)
    with pytest.raises(TypeError, match=r"tables must be dict, not NoneType"):
        utils.save_tables("tables/test.npz", None)
    utils.save_tables("tables/test.npz", {"test": np.array([0])})
    with pytest.raises(TypeError, match=r"path must be str or Path, not NoneType"):
        utils.load_tables(None)
    assert np.all(utils.load_tables("tables/test.npz")["test"] == [0])

    # get tables
    with pytest.raises(TypeError, match=r"filename must be str, not NoneType"):
        utils.get_tables(None, None, None, None)
    with pytest.raises(TypeError, match=r"tables_defs must be list, not NoneType"):
        utils.get_tables("", None, None, None)
    with pytest.raises(TypeError, match=r"generate_table_fn must be Callable, not NoneType"):
        utils.get_tables("", [None], None, None)
    def fn(coord_size, **kwargs):
        return np.zeros(coord_size, dtype=int)
    with pytest.raises(TypeError, match=r"accumulate must be bool, not NoneType"):
        utils.get_tables("", [None], fn, None)
    with pytest.raises(TypeError, match=r"tables_defs elements must be TableDef, not NoneType"):
        utils.get_tables("", [None], fn)
    tables = utils.get_tables("test", [TransitionDef("test", 1)], fn)
    assert np.all(tables["test"] == [0])
    tables = utils.get_tables("test", [TransitionDef("test", 2), TransitionDef("TEST", 2)], fn)
    assert np.all(tables["test"] == [0])
    assert np.all(tables["TEST"] == [0, 0])
    tables = utils.get_tables("test", [TransitionDef("TEST", 1)], fn, True)
    assert np.all(tables["test"] == [0])
    assert np.all(tables["TEST"] == [0, 0])
    tables = utils.get_tables("test", [TransitionDef("TEST", 1)], fn)
    assert "test" not in tables.keys()
    assert np.all(tables["TEST"] == [0, 0])
    os.remove("tables/test")
    os.remove("tables/test.npz")

    # generate transition table
    class TestSolver(BaseSolver):
        pruning_defs = [[PruningDef("co", 729, 0)]]
        @staticmethod
        def phase_coords(coords: Tuple[int, int], phase: int) -> Tuple[int, ...]:
            return coords
    solver = TestSolver()
    with pytest.raises(TypeError, match=r"coord_name must be str, not NoneType"):
        utils.generate_transition_table(None, None)
    with pytest.raises(TypeError, match=r"coord_size must be int, not NoneType"):
        utils.generate_transition_table("co", None)
    with pytest.raises(ValueError, match=r"coord_size must be > 0 and <= 65536 \(got 0\)"):
        utils.generate_transition_table("co", 0)
    with pytest.raises(ValueError, match=r"coord_size must be > 0 and <= 65536 \(got 65537\)"):
        utils.generate_transition_table("co", 65537)
    assert np.all(utils.generate_transition_table("co", 729) == solver.transition_tables["co"])

    # generate pruning table
    with pytest.raises(TypeError, match=r"phase must be int, not NoneType"):
        utils.generate_pruning_table(solver, None, None, "")
    with pytest.raises(TypeError, match=r"shape must be int or tuple, not NoneType"):
        utils.generate_pruning_table(solver, -1, None, "")
    with pytest.raises(TypeError, match=r"indexes must be int or tuple or None, not str"):
        utils.generate_pruning_table(solver, -1, (None,), "")
    with pytest.raises(ValueError, match=r"phase must be >= 0 and < 1 \(got -1\)"):
        utils.generate_pruning_table(solver, -1, (None,), (None,))
    with pytest.raises(ValueError, match=r"phase must be >= 0 and < 1 \(got 1\)"):
        utils.generate_pruning_table(solver, 1, (None,), (None,))
    with pytest.raises(TypeError, match=r"shape elements must be int, not NoneType"):
        utils.generate_pruning_table(solver, 0, (None,), (None,))
    with pytest.raises(TypeError, match=r"indexes elements must be int, not NoneType"):
        utils.generate_pruning_table(solver, 0, (729,), (None,))
    assert np.all(utils.generate_pruning_table(solver, 0, (729,), (0,)) == solver.pruning_tables["co"])
    assert np.all(utils.generate_pruning_table(solver, 0, 729, 0) == solver.pruning_tables["co"])
    os.remove("tables/pruning_testsolver.npz")


def test_solver():
    match = r"invalid string representation, setting undefined orientation and permutation values with -1"
    with pytest.warns(UserWarning, match=match):
        cube = Cube(repr="OYBRBYBWOGORYGGOOYWWBWRG")
    class TestSolver(BaseSolver):
        pass
    with pytest.raises(TypeError):
        TestSolver(None)
    class TestSolver(BaseSolver):
        @staticmethod
        def phase_coords(coords: Tuple[int, int], phase: int) -> Tuple[int, ...]:
            return coords
    with pytest.raises(TypeError, match=r"use_transition_tables must be bool, not NoneType"):
        TestSolver(None, None)
    with pytest.raises(TypeError, match=r"use_pruning_tables must be bool, not NoneType"):
        TestSolver(False, None)
    solver = TestSolver(False, False)
    assert repr(solver) == "TestSolver"
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
    with pytest.warns(UserWarning, match=r"invalid corner permutation"):
        cube = Cube(repr="BYBRRYBWOGORYGGOOWWWBWRG")
    with pytest.raises(ValueError, match=r"invalid cube state"):
        solver.solve(cube, 0, False, 0)
    scramble = Maneuver("U F2 R'")
    cube = Cube(scramble)
    assert solver.solve(cube, 0, False, 0) is None
    solution = solver.solve(cube, None, False, 0)
    assert isinstance(solution, Maneuver)
    assert apply_maneuver(cube, solution).coords == (0, 0)
    solution = solver.solve(cube, None, True, 0)
    assert isinstance(solution, Maneuver)
    assert apply_maneuver(cube, solution).coords == (0, 0)
    solution = solver.solve(cube, None, True, 1)
    assert isinstance(solution, Maneuver)
    assert apply_maneuver(cube, solution).coords == (0, 0)
    solution = solver.solve(cube, None, True, 2)
    assert isinstance(solution, list)
    assert apply_maneuver(cube, solution[0]).coords == (0, 0)
    solver = TestSolver(True, False)
    assert solver.transition_tables != {}
    assert solver.pruning_tables == {}
    solution = solver.solve(cube)
    assert isinstance(solution, Maneuver)
    assert apply_maneuver(cube, solution).coords == (0, 0)
    solver = TestSolver(False, True)
    assert solver.transition_tables == {}
    assert solver.pruning_tables == {}
    solution = solver.solve(cube)
    assert isinstance(solution, Maneuver)
    assert apply_maneuver(cube, solution).coords == (0, 0)
    class TestSolver(BaseSolver):
        pruning_defs = [[PruningDef("co", 729, 0)]]
        @staticmethod
        def phase_coords(coords: Tuple[int, int], phase: int) -> Tuple[int, ...]:
            return coords
    solver = TestSolver(False, True)
    assert solver.transition_tables == {}
    assert solver.pruning_tables != {}
    solution = solver.solve(cube)
    assert isinstance(solution, Maneuver)
    assert apply_maneuver(cube, solution).coords == (0, 0)
    solver = TestSolver(True, True)
    assert solver.transition_tables != {}
    assert solver.pruning_tables != {}
    solution = solver.solve(cube)
    assert isinstance(solution, Maneuver)
    assert apply_maneuver(cube, solution).coords == (0, 0)
    os.remove("tables/pruning_testsolver.npz")


def test_dummy():
    DummySolver()
    num_cubes = 12

    max_depth = 3
    solver = DummySolver(False, False)
    depth_test(solver, max_depth, num_cubes)

    max_depth = 5
    solver = DummySolver(True, False)
    depth_test(solver, max_depth, num_cubes)


def test_kociemba():
    num_cubes = 100

    solver = Kociemba(True)
    with pytest.raises(ValueError, match=r"phase must be >= 0 and < 2 \(got -1\)"):
        solver.phase_coords((), -1)
    solve_test(solver, num_cubes)


def test_korf():
    num_cubes = 100

    solver = Korf(True)
    solve_test(solver, num_cubes)
