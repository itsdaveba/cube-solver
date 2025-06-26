"""Common definitions."""
from .cube.enums import Move, Layer

CoordType = int | tuple[int, ...]
CoordsType = tuple[CoordType, ...]

face_moves = [*Move.face_moves()]
face_moves.reverse()
main_layers = (Layer.UP, Layer.FRONT, Layer.RIGHT)
NEXT_MOVES = {Move.NONE: face_moves}
NEXT_MOVES.update({move: [next for next in face_moves if move.axis != next.axis or
                          (move.layers[0] != next.layers[0] and move.layers[0] in main_layers)] for move in face_moves})
