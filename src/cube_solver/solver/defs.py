from dataclasses import dataclass

NONE = -1

FlattenCoords = tuple[int, ...]


@dataclass
class TransitionDef:
    coord_name: str
    coord_size: int

    @property
    def name(self) -> str:
        return self.coord_name


@dataclass
class PruningDef:
    name: str
    shape: int | tuple[int, ...]
    indexes: int | tuple[int, ...] | None
    phase: int | None = None


TableDef = TransitionDef | PruningDef
