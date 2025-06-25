from typing import TypedDict

NONE = -1

FlattenCoords = tuple[int, ...]


class TransitionDef(TypedDict):
    coord_name: str
    coord_size: int


class PruningDef(TypedDict):
    name: str
    shape: tuple[int, ...]
    indexes: tuple[int, ...] | None
