"""Solver definitions."""
from __future__ import annotations

from typing import Union, Tuple
from dataclasses import dataclass
from typing_extensions import TYPE_CHECKING

if TYPE_CHECKING:
    from .solver import BaseSolver

NONE = -1

FlattenCoords = Tuple[int, ...]


@dataclass
class TransitionDef:
    """
    Transition table definition.

    Parameters
    ----------
    coord_name : str
        Cube coordinate name.
    coord_size : int
        Cube coordinate size.
    """
    coord_name: str  #: :meta private:
    coord_size: int  #: :meta private:

    @property
    def name(self) -> str:
        """Transition table name (same as ``coord_name``)."""
        return self.coord_name


@dataclass
class PruningDef:
    """
    Pruning table definition.

    Parameters
    ----------
    name : str
        Pruning table name.
    shape : int or tuple of int
        Pruning table shape.
    indexes : int or tuple of int or None, optional
        Index or indexes of the phase coordinates to use for the pruning table.
        If ``None``, use all the phase coordinates.
    solver : BaseSolver or None, optional
        Solver object. Default is ``None``.
    phase : int or None, optional
        Solver phase (0-indexed). Default is ``None``.
    """
    name: str  #: :meta private:
    shape: Union[int, Tuple[int, ...]]  #: :meta private:
    indexes: Union[int, Tuple[int, ...], None] = None  #: :meta private:
    solver: Union[BaseSolver, None] = None  #: :meta private:
    phase: Union[int, None] = None  #: :meta private:


TableDef = Union[TransitionDef, PruningDef]
"""Table definition."""
