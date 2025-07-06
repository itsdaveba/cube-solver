import os
import time
import pytest
import numpy as np

from cube_solver import Cube, Maneuver, apply_maneuver
from cube_solver import BaseSolver, DummySolver, Korf, Thistlethwaite, Kociemba
from cube_solver.solver.defs import FlattenCoords, PruningDef
from cube_solver.solver import utils


def num_nodes(depth: int) -> int:
    if depth <= 2:
        return 104 * depth * depth - 87 * depth + 1
    return num_nodes(depth - 1) * 12 + num_nodes(depth - 2) * 18


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
        assert solution is not None and not isinstance(solution, list)
        lengths.append(len(solution))
        assert apply_maneuver(cube, solution).is_solved
        if i and i % 50 == 0:
            print(f" {i}")
        print("#", end="")
    print(f" {num_cubes}")
    nodes = np.mean(nodes, axis=0)  # TODO add solve length
    print("\nMean Time\tMean Solution Length\tMean Phase Nodes")
    print("--------------------------------------------------------")
    print(f"{np.mean(times):.4f} sec", end="\t")
    print(f"{np.mean(lengths)} moves", end="\t\t")
    for phase_nodes in nodes:
        print(f"{phase_nodes:.0f}", end=" ")
    print()


def depth_test(solver: BaseSolver, max_depth: int, num_cubes: int):  # TODO add stats about nodes and checks
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
            assert solution == scramble.inverse
            print("#", end="")
        nodes = np.mean(nodes, axis=0)
        print(f"\t{np.mean(times):.4f} sec", end="\t")
        for phase_nodes in nodes:
            print(f"{phase_nodes:.0f}", end=" ")
        print()


def test_utils():
    # flatten and unflatten
    cube = Cube()
    coords = cube.get_coords(False, False)
    flatten_coords = utils.flatten(coords)
    assert utils.unflatten(flatten_coords, False, False) == coords
    coords = cube.get_coords(False, True)
    flatten_coords = utils.flatten(coords)
    assert utils.unflatten(flatten_coords, False, True) == coords
    coords = cube.get_coords(True, False)
    flatten_coords = utils.flatten(coords)
    assert utils.unflatten(flatten_coords, True, False) == coords
    coords = cube.get_coords(True, True)
    flatten_coords = utils.flatten(coords)
    assert flatten_coords == (0, 0, 0, 1656, 0, 11856, 1656)
    assert utils.unflatten(flatten_coords, True, True) == coords

    # select
    assert utils.select(flatten_coords, None) == flatten_coords
    assert utils.select(flatten_coords, 5) == (11856,)
    assert utils.select(flatten_coords, (3, 5, 6)) == (1656, 11856, 1656)

    # load and save
    with pytest.raises(TypeError, match=r"path must be str or Path, not NoneType"):
        utils.save_tables(None, None)
    with pytest.raises(TypeError, match=r"tables must be dict, not NoneType"):
        utils.save_tables("tables/test.npz", None)
    utils.save_tables("tables/test.npz", {"test": np.array([0])})
    with pytest.raises(TypeError, match=r"path must be str or Path, not NoneType"):
        utils.load_tables(None)
    assert utils.load_tables("tables/test.npz") == {"test": [0]}

    # get tables
    with pytest.raises(TypeError, match=r"filename must be str, not NoneType"):
        utils.get_tables(None, None, None, None)
    with pytest.raises(TypeError, match=r"tables_defs must be list, not NoneType"):
        utils.get_tables("", None, None, None)
    with pytest.raises(TypeError, match=r"generate_table_fn must be Callable, not NoneType"):
        utils.get_tables("", [None], None, None)
    def fn(**kwargs):
        return np.array([0])
    with pytest.raises(TypeError, match=r"accumulate must be bool, not NoneType"):
        utils.get_tables("", [None], fn, None)
    with pytest.raises(TypeError, match=r"tables_defs elements must be TableDef, not NoneType"):
        utils.get_tables("", [None], fn)
    tables = utils.get_tables("test", [PruningDef("test")], fn)
    assert tables == {"test": [0]}
    tables = utils.get_tables("test", [PruningDef("test"), PruningDef("TEST")], fn)
    assert tables == {"test": [0], "TEST": [0]}
    tables = utils.get_tables("test", [PruningDef("TEST")], fn, True)
    assert tables == {"test": [0], "TEST": [0]}
    tables = utils.get_tables("test", [PruningDef("TEST")], fn)
    assert tables == {"TEST": [0]}
    tables = utils.get_tables("test", [PruningDef("TEST")], fn)
    assert tables == {"TEST": [0]}
    os.remove("tables/test")
    os.remove("tables/test.npz")

    # generate transition table
    class TestSolver(BaseSolver):
        partial_corner_perm = True
        partial_edge_perm = False
        @staticmethod
        def get_phase_coords(coords: FlattenCoords, phase: int) -> FlattenCoords:  # TODO change name
            return coords
        @staticmethod  # TODO check coords type and pahse_coords type
        def get_coords_from_phase_coords(phase_coords: FlattenCoords, phase: int) -> FlattenCoords:  # TODO change name
            return phase_coords
    with pytest.raises(ValueError, match=r"size must be <= 65536 \(got 239500800\)"):
        solver = TestSolver()
    class TestSolver(BaseSolver):
        partial_corner_perm = True
        partial_edge_perm = True
        pruning_defs = [[PruningDef("eo", indexes=1)]]
        @staticmethod
        def get_phase_coords(coords: FlattenCoords, phase: int) -> FlattenCoords:  # TODO change name
            return coords
        @staticmethod  # TODO check coords type and pahse_coords type
        def get_coords_from_phase_coords(phase_coords: FlattenCoords, phase: int) -> FlattenCoords:  # TODO change name
            return phase_coords
    solver = TestSolver()
    with pytest.raises(TypeError, match=r"phase must be int, not NoneType"):
        utils.generate_transition_table(solver, None, None)
    with pytest.raises(TypeError, match=r"index must be int, not NoneType"):
        utils.generate_transition_table(solver, -2, None)
    with pytest.raises(ValueError, match=r"phase must be >= -1 and < 1 \(got -2\)"):
        utils.generate_transition_table(solver, -2, -1)
    with pytest.raises(ValueError, match=r"phase must be >= -1 and < 1 \(got 1\)"):
        utils.generate_transition_table(solver, 1, -1)
    with pytest.raises(ValueError, match=r"index must be >= 0 and < 4 \(got -1\)"):
        utils.generate_transition_table(solver, -1, -1)
    with pytest.raises(ValueError, match=r"index must be >= 0 and < 4 \(got 4\)"):
        utils.generate_transition_table(solver, -1, 4)
    with pytest.raises(ValueError, match=r"missing required keyword argument: 'name'"):
        utils.generate_transition_table(solver, -1, 2)
    assert np.all(utils.generate_transition_table(solver, -1, 2, name="pcp") == solver.transition_tables["pcp"])
    assert np.all(utils.generate_transition_table(solver, 0, 1) == solver.transition_tables["eo"])

    # generate pruning table
    with pytest.raises(TypeError, match=r"phase must be int, not NoneType"):
        utils.generate_pruning_table(solver, None, "")
    with pytest.raises(TypeError, match=r"indexes must be int, tuple, or None, not str"):
        utils.generate_pruning_table(solver, -1, "")
    with pytest.raises(ValueError, match=r"phase must be >= 0 and < 1 \(got -1\)"):
        utils.generate_pruning_table(solver, -1, (None,))
    with pytest.raises(ValueError, match=r"phase must be >= 0 and < 1 \(got 1\)"):
        utils.generate_pruning_table(solver, 1, (None,))
    with pytest.raises(TypeError, match=r"indexes elements must be int, not NoneType"):
        utils.generate_pruning_table(solver, 0, (None,))
    with pytest.raises(ValueError, match=r"indexes elements must be >= 0 and < 4 \(got \(-1,\)\)"):
        utils.generate_pruning_table(solver, 0, (-1,))
    with pytest.raises(ValueError, match=r"indexes elements must be >= 0 and < 4 \(got \(4,\)\)"):
        utils.generate_pruning_table(solver, 0, (4,))
    with pytest.raises(ValueError, match=r"indexes must be >= 0 and < 4 \(got -1\)"):
        utils.generate_pruning_table(solver, 0, -1)
    with pytest.raises(ValueError, match=r"indexes must be >= 0 and < 4 \(got 4\)"):
        utils.generate_pruning_table(solver, 0, 4)
    assert np.all(utils.generate_pruning_table(solver, 0, (1,)) == solver.pruning_tables["eo"])
    assert np.all(utils.generate_pruning_table(solver, 0, 1) == solver.pruning_tables["eo"])
    os.remove("tables/pruning_testsolver.npz")


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
        solver.solve(None, "", None, "", None)
    with pytest.raises(TypeError, match=r"max_length must be int or None, not str"):
        solver.solve(cube, "", None, "", None)
    with pytest.raises(TypeError, match=r"optimal must be bool, not NoneType"):
        solver.solve(cube, -1, None, "", None)
    with pytest.raises(TypeError, match=r"timeout must be int or None, not str"):
        solver.solve(cube, -1, False, "", None)
    with pytest.raises(TypeError, match=r"verbose must be int, not NoneType"):
        solver.solve(cube, -1, False, -1, None)
    with pytest.raises(ValueError, match=r"max_length must be >= 0 \(got -1\)"):
        solver.solve(cube, -1, False, -1, -1)
    with pytest.raises(ValueError, match=r"timeout must be >= 0 \(got -1\)"):
        solver.solve(cube, 0, False, -1, -1)
    with pytest.raises(ValueError, match=r"verbose must be one of 0, 1, 2 \(got -1\)"):
        solver.solve(cube, 0, False, 0, -1)
    with pytest.raises(ValueError, match=r"invalid cube state"):
        solver.solve(cube, 0, False, 0, 0)
    with pytest.warns(UserWarning, match=r"invalid cube parity"):
        cube.set_coord("pep", (0, 11857, 1656))
    with pytest.warns(UserWarning, match=r"invalid cube parity"):
        with pytest.raises(ValueError, match=r"invalid cube state"):
            solver.solve(cube, 0, False, 0, 0)
    with pytest.raises(ValueError, match=r"coord_size must be <= 65536 \(got 239500800\)"):
        solver.generate_transition_table("ep", 239500800)
    scramble = Maneuver("U F2 R'")
    cube = Cube(scramble)
    assert solver.solve(cube, 0, False, 0, 0) is None
    assert solver.terminate
    assert solver.solve(cube, 0, False, None, 0) is None
    assert not solver.terminate
    assert solver.solve(cube, None, True, None, 0) == scramble.inverse
    assert not solver.terminate
    assert solver.solve(cube, None, True, None, 1) == scramble.inverse
    assert not solver.terminate
    assert solver.solve(cube, None, True, None, 2) == [scramble.inverse]
    assert not solver.terminate
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
    assert np.all(solver.generate_transition_table("pcp", 1680) == solver.transition_tables["pcp"])
    tables = utils.load_tables("tables/transition.npz")
    utils.save_tables("tables/transition.npz", tables)
    os.remove("tables/pruning_testsolver.npz")


def test_dummy():
    num_cubes = 12

    max_depth = 3
    solver = DummySolver(use_transition_tables=False)
    depth_test(solver, max_depth, num_cubes)

    max_depth = 4
    solver = DummySolver()
    depth_test(solver, max_depth, num_cubes)


def test_korf():
    num_cubes = 12

    max_depth = 4
    solver = Korf(use_pruning_tables=False)
    depth_test(solver, max_depth, num_cubes)

    max_depth = 11
    solver = Korf()
    depth_test(solver, max_depth, num_cubes)


def test_thistlethwaite():
    num_cubes = 100
    solver = Thistlethwaite()
    with pytest.raises(ValueError, match=r"phase must be < 4 \(got 4\)"):
        solver.phase_coords((), 4)
    solve_test(solver, num_cubes)


def test_kociemba():
    num_cubes = 100
    solver = Kociemba()
    with pytest.raises(ValueError, match=r"phase must be < 2 \(got 2\)"):
        solver.phase_coords((), 2)
    solve_test(solver, num_cubes)
