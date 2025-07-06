"""Solver utils module."""
import numpy as np
from pathlib import Path
from dataclasses import asdict
from typing import Sequence, Callable

from ..logger import logger
from .defs import TableDef


# TODO document

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


def get_tables(filename: str, tables_kwargs: Sequence[TableDef],
               generate_table_fn: Callable[..., np.ndarray], accumulate: bool) -> dict[str, np.ndarray]:
    path = Path(f"tables/{filename}")
    try:
        tables = load_tables(path)
        save = False
        for kwargs in tables_kwargs:
            if kwargs.name not in tables:
                logger.info(f"Updating {path}")
                tables[kwargs.name] = generate_table_fn(**asdict(kwargs))
                save = True
        if not accumulate:
            names = {kwargs.name for kwargs in tables_kwargs}
            for name in tables.keys() - names:
                del tables[name]
                save = True
        if save:
            save_tables(path, tables)
    except FileNotFoundError:
        logger.info(f"Creating {path}")
        tables = {kwargs.name: generate_table_fn(**asdict(kwargs)) for kwargs in tables_kwargs}
        save_tables(path, tables)
    return tables
