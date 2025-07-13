"""Common definitions."""
from .cube.enums import Move


face_moves = [*Move.face_moves()]
NEXT_MOVES = {Move.NONE: face_moves}
NEXT_MOVES.update({move: [next for next in face_moves if move.axis != next.axis] for move in face_moves})
