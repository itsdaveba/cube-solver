import numpy as np
from pathlib import Path


def load_tables(path: Path) -> dict[str, np.ndarray]:
    with path.open("rb") as file:
        tables = np.load(file, allow_pickle=False)
        tables = dict(tables)
    return tables


def save_tables(path: Path, tables: dict[str, np.ndarray]):
    path.parent.mkdir(exist_ok=True)
    with path.open("wb") as file:
        np.savez(file, allow_pickle=False, **tables)


# TODO def get_tables()
