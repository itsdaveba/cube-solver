"""Common definitions."""
from typing import Union, Tuple

from .cube.enums import Move, Layer

CoordType = Union[int, Tuple[int, ...]]
"""Cube coordinate type."""
CoordsType = Tuple[CoordType, ...]
"""Cube coordinates type."""

face_moves = [*Move.face_moves()]
main_layers = (Layer.UP, Layer.FRONT, Layer.RIGHT)
NEXT_MOVES = {Move.NONE: face_moves}
NEXT_MOVES.update({move: [next for next in face_moves if move.axis != next.axis or
                          (move.layers[0] != next.layers[0] and move.layers[0] in main_layers)] for move in face_moves})
