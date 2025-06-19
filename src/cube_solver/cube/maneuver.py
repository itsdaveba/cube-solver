from __future__ import annotations

import random
from typing import Iterator, overload
from collections import Counter
from itertools import combinations

from .cube import Cube
from .enums import Move, Layer

face_moves = [*Move.face_moves()]
face_moves.reverse()
main_layers = (Layer.UP, Layer.FRONT, Layer.RIGHT)
NEXT_MOVES = {move: [next for next in face_moves if move.axis != next.axis or
                     (move.layers[0] != next.layers[0] and move.layers[0] in main_layers)] for move in face_moves}
NEXT_MOVES[Move.NONE] = face_moves

MOVES_DICT = {frozenset(): []}
for first, second in combinations(Move, 2):
    if first.axis == second.axis and first.name[:-1] != second.name[:-1]:
        counter = Counter(dict(zip(first.layers, first.shifts)))
        counter.update(dict(zip(second.layers, second.shifts)))
        key = frozenset((layer, shift % 4) for layer, shift in counter.items() if shift % 4)
        if key not in MOVES_DICT:
            MOVES_DICT[key] = [first, second]
    elif first == Move.NONE:
        counter = Counter(dict(zip(second.layers, second.shifts)))
        key = frozenset((layer, shift % 4) for layer, shift in counter.items())
        MOVES_DICT[key] = [second]


def reduce_moves(moves: list[Move]) -> list[Move]:
    n = len(moves)
    for i in range(n - 1):
        if moves[i].axis == moves[i+1].axis:
            counter = Counter(dict(zip(moves[i].layers, moves[i].shifts)))
            counter.update(dict(zip(moves[i+1].layers, moves[i+1].shifts)))
            key = frozenset((layer, shift % 4) for layer, shift in counter.items() if shift % 4)
            if key in MOVES_DICT and len(MOVES_DICT[key]) < 2:
                return reduce_moves(moves[:i] + MOVES_DICT[key] + moves[i+2:])
            if i + 2 < n and moves[i].axis == moves[i+2].axis:
                counter.update(dict(zip(moves[i+2].layers, moves[i+2].shifts)))
                key = frozenset((layer, shift % 4) for layer, shift in counter.items() if shift % 4)
                if key in MOVES_DICT:
                    return reduce_moves(moves[:i] + MOVES_DICT[key] + moves[i+3:])
                if i + 3 < n and moves[i].axis == moves[i+3].axis:
                    counter.update(dict(zip(moves[i+3].layers, moves[i+3].shifts)))
                    key = frozenset((layer, shift % 4) for layer, shift in counter.items() if shift % 4)
                    return reduce_moves(moves[:i] + MOVES_DICT[key] + moves[i+4:])
    return moves


class Maneuver(str):
    def __new__(cls, moves: str | list[Move]) -> Maneuver:
        """
        Maneuver.

        Parameters
        ---------
        moves : str or list of Move
            Maneuver moves.
        """
        cls.moves: tuple[Move, ...]

        if not isinstance(moves, (str, list)):
            raise TypeError(f"moves must be str or list, not {type(moves).__name__}")

        if isinstance(moves, str):
            moves = [Move.from_string(move_str) for move_str in moves.split()]
        else:
            for move in moves:
                if not isinstance(move, Move):
                    raise TypeError(f"moves list elements must be Move, not {type(move).__name__}")

        moves = reduce_moves([move for move in moves if move != Move.NONE])
        obj = super().__new__(cls, " ".join(move.string for move in moves))
        obj.moves = tuple(moves)
        return obj

    def __ne__(self, other: str | list[Move]) -> bool:
        return not self.__eq__(other)

    def __eq__(self, other: str | list[Move]) -> bool:
        """Maneuver equality comparison."""
        if not isinstance(other, Maneuver):
            try:
                other = Maneuver(other)
            except Exception:
                return False
        return repr(Cube(self)) == repr(Cube(other))

    def __len__(self) -> int:
        return len(self.moves)

    @overload
    def __getitem__(self, key: int) -> Move: ...
    @overload
    def __getitem__(self, key: slice) -> Maneuver: ...

    def __getitem__(self, key: int | slice) -> Move | Maneuver:
        if not isinstance(key, (int, slice)):
            raise TypeError(f"Maneuver indices must be int or slice, not {type(key).__name__}")
        if isinstance(key, int):
            try:
                return self.moves[key]
            except IndexError:
                raise IndexError("Maneuver index out of range")
        return Maneuver([*self.moves[key]])

    def __iter__(self) -> Iterator[Move]:
        for move in self.moves:
            yield move

    def __contains__(self, item: str | Move) -> bool:
        if isinstance(item, str):
            item = Move.from_string(item)
        return item in self.moves

    def __neg__(self) -> Maneuver:
        return self.inverse

    def __add__(self, other: str | list[Move]) -> Maneuver:
        if not isinstance(other, Maneuver):
            other = Maneuver(other)
        return Maneuver([*self.moves] + [*other.moves])

    def __radd__(self, other: str | list[Move]) -> Maneuver:
        other = Maneuver(other)
        return other.__add__(self)

    def __sub__(self, other: str | list[Move]) -> Maneuver:
        if not isinstance(other, Maneuver):
            other = Maneuver(other)
        return self.__add__(other.__neg__())

    def __rsub__(self, other: str | list[Move]) -> Maneuver:
        other = Maneuver(other)
        return other.__sub__(self)

    def __mul__(self, other: int | str | list[Move]) -> Maneuver:
        if isinstance(other, int):
            return Maneuver([*self.moves] * other)
        if not isinstance(other, Maneuver):
            other = Maneuver(other)
        return self.__add__(other).__sub__(self)

    def __rmul__(self, other: int | str | list[Move]) -> Maneuver:
        if isinstance(other, int):
            return Maneuver([*self.moves] * other)
        other = Maneuver(other)
        return other.__mul__(self)

    def __matmul__(self, other: str | list[Move]) -> Maneuver:
        if not isinstance(other, Maneuver):
            other = Maneuver(other)
        return self.__mul__(other).__sub__(other)

    def __rmatmul__(self, other: str | list[Move]) -> Maneuver:
        other = Maneuver(other)
        return other.__matmul__(self)

    @classmethod
    def random(cls, length: int = 25) -> Maneuver:  # TODO include other move types?
        if not isinstance(length, int):
            raise TypeError(f"length must be int, not {type(length).__name__}")
        if length < 0:
            raise ValueError(f"length must be >= 0 (got {length})")

        moves = [Move.NONE]
        for i in range(length):
            moves.append(random.choice(NEXT_MOVES[moves[-1]]))
        return cls(moves[1:])

    @property
    def inverse(self) -> Maneuver:
        """Inverse maneuver."""
        return Maneuver([move.inverse for move in self.moves[::-1]])
