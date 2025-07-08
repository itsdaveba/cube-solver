"""Solver utils module."""
from __future__ import annotations

import numpy as np
from pathlib import Path
from collections import deque
from dataclasses import asdict
from typing_extensions import TYPE_CHECKING
from typing import Union, Tuple, Sequence, Dict, Callable

from ..logger import logger
from ..defs import CoordsType, NEXT_MOVES
from ..cube.enums import Move
from ..cube.cube import Cube, apply_move
from .defs import NONE, FlattenCoords, TransitionDef, PruningDef, TableDef

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

    Examples
    --------
    >>> from cube_solver import Cube
    >>> from cube_solver.solver import utils
    >>> cube = Cube()
    >>> coords = cube.get_coords(partial_corner_perm=True, partial_edge_perm=True)
    >>> coords
    (0, 0, (0, 1656), (0, 11856, 1656))
    >>> utils.flatten(coords)
    (0, 0, 0, 1656, 0, 11856, 1656)
    """
    flatten_coords = ()
    for coord in coords:
        flatten_coords += (coord,) if isinstance(coord, int) else coord
    return flatten_coords


def select(coords: FlattenCoords, indexes: Union[int, Tuple[int, ...], None]) -> FlattenCoords:
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

    Examples
    --------
    >>> from cube_solver import Cube
    >>> from cube_solver.solver import utils
    >>> cube = Cube()
    >>> coords = cube.get_coords(partial_corner_perm=True, partial_edge_perm=True)
    >>> flatten_coords = utils.flatten(coords)
    >>> utils.select(flatten_coords, None)
    (0, 0, 0, 1656, 0, 11856, 1656)
    >>> utils.select(flatten_coords, 5)
    (11856,)
    >>> utils.select(flatten_coords, (3, 5, 6))
    (1656, 11856, 1656)
    """
    if indexes is None:
        return coords
    if isinstance(indexes, int):
        return (coords[indexes],)
    return tuple(coords[index] for index in indexes)


def load_tables(path: Union[str, Path]) -> Dict[str, np.ndarray]:
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


def save_tables(path: Union[str, Path], tables: Dict[str, np.ndarray]):
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
        np.savez(file, **tables)  # type: ignore


def get_tables(filename: str, tables_defs: Sequence[TableDef],
               generate_table_fn: Callable[..., np.ndarray], accumulate: bool = False) -> Dict[str, np.ndarray]:
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
    accumulate : bool, optional
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
        if not isinstance(kwargs, (TransitionDef, PruningDef)):
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


def generate_transition_table(coord_name: str, coord_size: int) -> np.ndarray:
    """
    Generate the cube coordinate transition table.

    Parameters
    ----------
    coord_name : str
        Cube coordinate name.
    coord_size : int
        Cube coordinate size.

    Returns
    -------
    transition_table : ndarray
        Cube coordinate transition table.
    """
    if not isinstance(coord_name, str):
        raise TypeError(f"coord_name must be str, not {type(coord_name).__name__}")
    if not isinstance(coord_size, int):
        raise TypeError(f"coord_size must be int, not {type(coord_size).__name__}")

    if coord_size <= 0 or coord_size - 1 > np.iinfo(np.uint16).max:
        raise ValueError(f"coord_size must be > 0 and <= {np.iinfo(np.uint16).max + 1} (got {coord_size})")
    transition_table = np.zeros((coord_size, len(NEXT_MOVES[Move.NONE])), dtype=np.uint16)

    cube = Cube()
    for coord in range(coord_size):
        cube.set_coord(coord_name, coord)
        transition_table[coord] = [apply_move(cube, move).get_coord(coord_name) for move in NEXT_MOVES[Move.NONE]]
    return transition_table


def generate_pruning_table(solver: BaseSolver, phase: int, shape: Union[int, Tuple[int, ...]],
                           indexes: Union[int, Tuple[int, ...], None], **kwargs) -> np.ndarray:
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
    if not isinstance(shape, (int, tuple)):
        raise TypeError(f"shape must be int or tuple, not {type(shape).__name__}")
    if indexes is not None and not isinstance(indexes, (int, tuple)):
        raise TypeError(f"indexes must be int or tuple or None, not {type(indexes).__name__}")
    if phase < 0 or phase >= solver.num_phases:
        raise ValueError(f"phase must be >= 0 and < {solver.num_phases} (got {phase})")
    if isinstance(shape, tuple):
        for size in shape:
            if not isinstance(size, int):
                raise TypeError(f"shape elements must be int, not {type(size).__name__}")
    if isinstance(indexes, tuple):
        for index in indexes:
            if not isinstance(index, int):
                raise TypeError(f"indexes elements must be int, not {type(index).__name__}")

    pruning_table = np.full(shape, NONE, dtype=np.int8)
    init_coords = solver.get_coords(Cube())
    phase_coords = solver.phase_coords(flatten(init_coords), phase)
    prune_coords = select(phase_coords, indexes)
    pruning_table[prune_coords] = 0
    queue = deque([(init_coords, 0)])
    while queue:
        coords, depth = queue.popleft()
        for move in solver.phase_moves[phase]:
            next_coords = solver.next_position(coords, move)
            phase_coords = solver.phase_coords(flatten(next_coords), phase)
            prune_coords = select(phase_coords, indexes)
            if pruning_table[prune_coords] == NONE:
                pruning_table[prune_coords] = depth + 1
                queue.append((next_coords, depth + 1))
    return pruning_table
