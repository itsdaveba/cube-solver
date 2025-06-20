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
NEXT_MOVES = {Move.NONE: face_moves}
NEXT_MOVES.update({move: [next for next in face_moves if move.axis != next.axis or
                          (move.layers[0] != next.layers[0] and move.layers[0] in main_layers)] for move in face_moves})

REDUCED_MOVES = {frozenset(): []}
for first, second in combinations(Move, 2):
    if first.axis == second.axis and first.name[:-1] != second.name[:-1]:
        counter = Counter(dict(zip(first.layers, first.shifts)))
        counter.update(dict(zip(second.layers, second.shifts)))
        key = frozenset((layer, shift % 4) for layer, shift in counter.items() if shift % 4)
        if key not in REDUCED_MOVES:
            REDUCED_MOVES[key] = [first, second]
    elif first == Move.NONE:
        counter = Counter(dict(zip(second.layers, second.shifts)))
        key = frozenset((layer, shift % 4) for layer, shift in counter.items())
        REDUCED_MOVES[key] = [second]


def _reduce_moves(moves: list[Move]) -> list[Move]:
    """Try to reduce consecutive moves along the same axis."""
    for i in range(len(moves) - 1):
        if moves[i].axis == moves[i+1].axis:
            counter = Counter(dict(zip(moves[i].layers, moves[i].shifts)))
            reduced_moves = _reduced(counter, moves, i, 1)
            if reduced_moves is not None:
                return reduced_moves
    return moves


def _reduced(counter: Counter, moves: list[Move], i: int, j: int) -> list[Move] | None:
    """Update counter and reduce if possible."""
    counter.update(dict(zip(moves[i+j].layers, moves[i+j].shifts)))
    key = frozenset((layer, shift % 4) for layer, shift in counter.items() if shift % 4)
    if key in REDUCED_MOVES and len(REDUCED_MOVES[key]) <= j:
        return _reduce_moves(moves[:i] + REDUCED_MOVES[key] + moves[i+j+1:])
    if i + j + 1 < len(moves) and moves[i].axis == moves[i+j+1].axis:
        return _reduced(counter, moves, i, j + 1)


class Maneuver(str):
    def __new__(cls, moves: str | list[Move]):
        """
        Create :class:`Maneuver` object.

        A :class:`Maneuver` is a subclass of :class:`str` that represents
        the sequence of moves that can be applied to a :class:`Cube` object.

        Parameters
        ----------
        moves : str or list of Move.
            Sequence of moves.

        Examples
        --------
        >>> from cube_solver import Cube, Move, Maneuver

        Apply a maneuver to a cube.

        >>> Maneuver("U F2 R'")
        "U F2 R'"
        >>> Maneuver([Move.U1, Move.F2, Move.R3])
        "U F2 R'"
        >>> Cube(Maneuver("U F2 R'")) == Cube("U F2 R'")
        True

        Compare equvalent maneuvers.

        >>> Maneuver("R L") == "L R"
        True
        >>> Maneuver("R L F2 B2 R' L' D R L F2 B2 R' L'") == "U"
        True

        Reduce consecutive moves along the same axis.

        >>> Maneuver("U U")
        'U2'
        >>> Maneuver("R L R")
        'R2 L'
        >>> Maneuver("F S' B2 F")
        'S z2'

        Inverse maneuver.

        >>> Maneuver("R U R' U'").inverse
        "U R U' R'"

        Container operations.

        >>> maneuver = Maneuver("U F2 R'")
        >>> maneuver[0]
        Move.U1
        >>> maneuver[:2]
        'U F2'
        >>> [move for move in maneuver]
        [Move.U1, Move.F2, Move.R3]
        >>> Move.U1 in maneuver
        True

        Numeric operations.

        >>> A = Maneuver("U R'")
        >>> B = Maneuver("F'")
        >>> -A     # '-' negation, same as A'
        "R U'"
        >>> A + B  # '+' addition, same as A B
        "U R' F'"
        >>> A - B  # '-' subtraction, same as A B'
        "U R' F"
        >>> 2 * A  # '*' multiplication, same as A A
        "U R' U R'"
        >>> A * B  # '*' conjugation, same as A B A'
        "U R' F' R U'"
        >>> A @ B  # '@' commutator, same as A B A' B'
        "U R' F' R U' F"
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

        moves = _reduce_moves([move for move in moves if move != Move.NONE])
        obj = super().__new__(cls, " ".join(move.string for move in moves))
        obj.moves = tuple(moves)
        return obj

    def __ne__(self, other: object) -> bool:
        return not self.__eq__(other)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Maneuver):
            try:
                if not isinstance(other, (str, list)):
                    return False
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

    @property
    def inverse(self) -> Maneuver:
        """
        Inverse maneuver.

        Examples
        --------
        >>> from cube_solver import Maneuver
        >>> Maneuver("R U R' U'").inverse
        "U R U' R'"
        """
        return Maneuver([move.inverse for move in self.moves[::-1]])

    @classmethod
    def random(cls, length: int = 25) -> Maneuver:  # TODO include other move types?
        """
        Generate a random maneuver.

        Parameters
        ----------
        length : int, optional
            Maneuver length. Default is ``25``.

        Returns
        -------
        maneuver : Maneuver
            Random maneuver of the specified length.

        Examples
        --------
        >>> from cube_solver import Maneuver
        >>> Maneuver.random(10)  # result might differ # doctest: +SKIP
        "D2 R2 L' B2 D L D F L D2"
        """
        if not isinstance(length, int):
            raise TypeError(f"length must be int, not {type(length).__name__}")
        if length < 0:
            raise ValueError(f"length must be >= 0 (got {length})")

        moves = [Move.NONE]
        for i in range(length):
            moves.append(random.choice(NEXT_MOVES[moves[-1]]))
        return cls(moves[1:])
