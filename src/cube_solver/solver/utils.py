"""Solver utils module."""
from __future__ import annotations

import numpy as np
from pathlib import Path
from collections import deque
from dataclasses import asdict
from typing import Sequence, Callable
from typing_extensions import TYPE_CHECKING

from ..logger import logger
from ..defs import CoordsType
from ..cube import Cube, apply_move
from .defs import NONE, FlattenCoords, TableDef

if TYPE_CHECKING:
    from .solver import BaseSolver


def flatten(coords: CoordsType) -> FlattenCoords:
    """
    Get the flattened cube coordinates.

    Parameters
    ----------
    coords : tuple of (int or tuple of int)
        Cube coordinates.

    Returns
    -------
    flatten_coords : tuple of int
        Flattened cube coordinates.
    """
    flatten_coords = ()
    for coord in coords:
        flatten_coords += (coord,) if isinstance(coord, int) else coord
    return flatten_coords


def unflatten(coords: FlattenCoords, partial_corner_perm: bool, partial_edge_perm: bool) -> CoordsType:
    """
    Get the unflattened cube coordinates.

    Parameters
    ----------
    coords : tuple of int
        Flattened cube coordinates.
    partial_corner_perm : bool
        If ``True``, unflatten the `partial corner permutation` coordinate,
        otherwise unflatten the normal `corner permutation` coordinate.
    partial_edge_perm : bool
        If ``True``, unflatten the `partial edge permutation` coordinate,
        otherwise unflatten the normal `edge permutation` coordinate.

    Returns
    -------
    unflatten_coords : tuple of (int or tuple of int)
        Unflattened cube coordinates.
    """
    unflatten_coords = tuple(coords[:2])
    unflatten_coords += (tuple(coords[2:4]),) if partial_corner_perm else (coords[2],)
    i = 4 if partial_corner_perm else 3
    unflatten_coords += (tuple(coords[i:]),) if partial_edge_perm else (coords[i],)
    return unflatten_coords


def select(coords: FlattenCoords, indexes: int | tuple[int, ...] | None) -> FlattenCoords:
    """
    Select coordinates.

    Parameters
    ----------
    coords : tuple of int
        Coordinates.
    indexes : int or tuple of int or None
        Index or indexes of the coordinates to select.
        If ``None``, all coordinates are selected.

    Returns
    -------
    selected_coords : tuple of int
        Selected coordinates.
    """
    if indexes is None:
        return coords
    if isinstance(indexes, int):
        return (coords[indexes],)
    return tuple(coords[index] for index in indexes)


def load_tables(path: str | Path) -> dict[str, np.ndarray]:
    """
    Load the tables from a file.

    Parameters
    ----------
    path : str or Path
        Path of the file.

    Returns
    -------
    tables : dict
        Dictionary containig the tables.
    """
    if not isinstance(path, (str, Path)):
        raise TypeError(f"path must be str or Path, not {type(path).__name__}")

    if isinstance(path, str):
        path = Path(path)
    with np.load(path, allow_pickle=False) as data:
        tables = dict(data)
    return tables


def save_tables(path: str | Path, tables: dict[str, np.ndarray]):
    """
    Save the tables into a single file.

    Parameters
    ----------
    path : str or Path
        Path of the file.
    tables : dict
        Dictionary containig the tables.
    """
    if not isinstance(path, (str, Path)):
        raise TypeError(f"path must be str or Path, not {type(path).__name__}")
    if not isinstance(tables, dict):
        raise TypeError(f"tables must be dict, not {type(tables).__name__}")

    if isinstance(path, str):
        path = Path(path)
    path.parent.mkdir(exist_ok=True)
    with path.open("wb") as file:
        np.savez(file, allow_pickle=False, **tables)


def get_tables(filename: str, tables_defs: Sequence[TableDef],
               generate_table_fn: Callable[..., np.ndarray], accumulate: bool = False) -> dict[str, np.ndarray]:
    """
    Create or load tables from the ``tables/`` directory according to the ``tables_defs``.

    If the file does not exist, or if it exists but is missing some tables,
    create the missing tables from ``tables_defs`` and update the file.

    Parameters
    ----------
    filename : str
        Name of the file in the ``tables/`` directory.
    tables_defs : list of TableDef
        Table definitions.
    generate_table_fn : Callable
        Function to generate the table.
        It must accept the TableDef keyword arguments.
    accumulate : bool, optional.
        Whether to keep the tables not included in ``tables_defs``.
        Default is ``False``.

    Returns
    -------
    tables : dict
        Dictionary containig the tables.
        The keys represent the name of the table from :attr:`TableDef.name`.
    """
    if not isinstance(filename, str):
        raise TypeError(f"filename must be str, not {type(filename).__name__}")
    if not isinstance(tables_defs, list):
        raise TypeError(f"tables_defs must be list, not {type(tables_defs).__name__}")
    if not isinstance(generate_table_fn, Callable):
        raise TypeError(f"generate_table_fn must be Callable, not {type(generate_table_fn).__name__}")
    if not isinstance(accumulate, bool):
        raise TypeError(f"accumulate must be bool, not {type(accumulate).__name__}")
    for kwargs in tables_defs:
        if not isinstance(kwargs, TableDef):
            raise TypeError(f"tables_defs elements must be TableDef, not {type(kwargs).__name__}")

    path = Path(f"tables/{filename}")
    try:
        tables = load_tables(path)
        save = False
        for kwargs in tables_defs:
            if kwargs.name not in tables:
                logger.info(f"Updating {path}")
                tables[kwargs.name] = generate_table_fn(**asdict(kwargs))
                save = True
        if not accumulate:
            names = {kwargs.name for kwargs in tables_defs}
            for name in tables.keys() - names:
                logger.info(f"Updating {path}")
                del tables[name]
                save = True
        if save:
            save_tables(path, tables)
    except FileNotFoundError:
        logger.info(f"Creating {path}")
        tables = {kwargs.name: generate_table_fn(**asdict(kwargs)) for kwargs in tables_defs}
        save_tables(path, tables)
    return tables


# TODO double check when phase == -1
def generate_transition_table(solver: BaseSolver, phase: int, index: int, name: str | None = None) -> np.ndarray:
    """
    Generate the phase coordinate transition table.

    Parameters
    ----------
    solver : BaseSolver
        Solver object.
    phase : int
        Solver phase (0-indexed).
        If ``-1``, use the cube coordinates instead of the phase coordinates.
    index : int
        Index of the phase coordinates to use for the transition table.
    name : str or None, optional
        Cube coordinate name if ``phase`` is ``-1``, ignored otherwise.
        Default is ``None``.

    Returns
    -------
    transition_table : ndarray
        Phase coordinate transition table.
    """
    if not isinstance(phase, int):
        raise TypeError(f"phase must be int, not {type(phase).__name__}")
    if not isinstance(index, int):
        raise TypeError(f"index must be int, not {type(index).__name__}")
    if phase < NONE or phase >= solver.num_phases:
        raise ValueError(f"phase must be >= -1 and < {solver.num_phases} (got {phase})")
    if index < 0 or index >= len(solver.phase_coords_sizes[phase]):
        raise ValueError(f"index must be >= 0 and < {len(solver.phase_coords_sizes[phase])} (got {index})")
    if phase == NONE and name is None:
        raise ValueError("missing required keyword argument: 'name'")

    size = solver.phase_coords_sizes[phase][index]
    if size - 1 > np.iinfo(np.uint16).max:
        raise ValueError(f"size must be <= {np.iinfo(np.uint16).max + 1} (got {size})")
    transition_table = np.zeros((size, len(solver.phase_moves[phase])), dtype=np.uint16)

    cube = Cube()
    if phase == NONE:
        assert name is not None
        for coord in range(size):
            cube.set_coord(name, coord)
            transition_table[coord] = [apply_move(cube, move).get_coord(name) for move in solver.phase_moves[phase]]
        return transition_table
    coords = flatten(cube.get_coords(solver.partial_corner_perm, solver.partial_edge_perm))
    phase_coords = [*solver.get_phase_coords(coords, phase)]
    for phase_coords[index] in range(size):
        coords = solver.get_coords_from_phase_coords(tuple(phase_coords), phase)
        coords = unflatten(coords, solver.partial_corner_perm, solver.partial_edge_perm)
        cube.set_coords(coords, solver.partial_corner_perm, solver.partial_edge_perm)
        for m, move in enumerate(solver.phase_moves[phase]):
            coords = apply_move(cube, move).get_coords(solver.partial_corner_perm, solver.partial_edge_perm)
            transition_table[phase_coords[index], m] = solver.get_phase_coords(flatten(coords), phase)[index]
    return transition_table


def generate_pruning_table(solver: BaseSolver, phase: int, indexes: int | tuple[int, ...] | None, **kwargs) -> np.ndarray:
    """
    Generate the phase coordinates pruning table.

    Parameters
    ----------
    solver : BaseSolver
        Solver object.
    phase : int
        Solver phase (0-indexed).
    shape : int or tuple of int
        Shape of the pruning table.
    indexes : int or tuple of int or None
        Index or indexes of the phase coordinates to use for the pruning table.
        If ``None``, use all the phase coordinates.

    Returns
    -------
    pruning_table : ndarray
        Phase coordinates pruning table.
    """
    if not isinstance(phase, int):
        raise TypeError(f"phase must be int, not {type(phase).__name__}")
    if indexes is not None and not isinstance(indexes, (int, tuple)):
        raise TypeError(f"indexes must be int, tuple, or None, not {type(indexes).__name__}")
    if phase < 0 or phase >= solver.num_phases:
        raise ValueError(f"phase must be >= 0 and < {solver.num_phases} (got {phase})")
    n = len(solver.phase_coords_sizes[phase])
    if isinstance(indexes, tuple):
        for index in indexes:
            if not isinstance(index, int):
                raise TypeError(f"indexes elements must be int, not {type(index).__name__}")
            if index < 0 or index >= n:
                raise ValueError(f"indexes elements must be >= 0 and < {n} (got {indexes})")
    elif indexes is not None and (indexes < 0 or indexes >= n):
        raise ValueError(f"indexes must be >= 0 and < {n} (got {indexes})")

    shape = select(solver.phase_coords_sizes[phase], indexes)
    pruning_table = np.full(shape, NONE, dtype=np.int8)

    cube = Cube()
    coords = cube.get_coords(solver.partial_corner_perm, solver.partial_edge_perm)
    phase_coords = solver.get_phase_coords(flatten(coords), phase)
    prune_coords = select(phase_coords, indexes)
    pruning_table[prune_coords] = 0
    queue = deque([(phase_coords, 0)])
    while queue:
        phase_coords, depth = queue.popleft()
        for move in solver.phase_moves[phase]:
            next_coords = flatten(solver.next_position(phase_coords, phase, move))
            prune_coords = select(next_coords, indexes)
            if pruning_table[prune_coords] == NONE:
                pruning_table[prune_coords] = depth + 1
                queue.append((next_coords, depth + 1))
    return pruning_table
