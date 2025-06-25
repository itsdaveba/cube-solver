import numpy as np
from pathlib import Path
from dataclasses import asdict
from typing import Sequence, Callable

from .defs import TableDef

# TODO document


def load_tables(path: Path) -> dict[str, np.ndarray]:
    with path.open("rb") as file:
        tables = np.load(file, allow_pickle=False)
        tables = dict(tables)
    return tables


def save_tables(path: Path, tables: dict[str, np.ndarray]):
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
        tables = {kwargs.name: generate_table_fn(**asdict(kwargs)) for kwargs in tables_kwargs}
        save_tables(path, tables)
    return tables
