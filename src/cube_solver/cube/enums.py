"""Enums module."""
# TODO document
from typing import Iterator
from enum import Enum


class IntEnum(int, Enum):
    def __repr__(self) -> str:
        return self.name


class Color(IntEnum):  # up, front, right, down, back, left
    WHITE = 0
    GREEN = 1
    ORANGE = 2
    YELLOW = 3
    BLUE = 4
    RED = 5
    BLACK = 6

    def __str__(self) -> str:
        return self.name[0]


class Face(IntEnum):
    U = 0
    F = 1
    R = 2
    D = 3
    B = 4
    L = 5

    @property
    def cubie_slice(self) -> tuple:
        return cubie_slice[self]

    @property
    def perm(self) -> list:
        return face_perm[self]


cubie_slice = {
    Face.U: (0, slice(None), slice(None), 0),
    Face.F: (slice(None), 2, slice(None), 1),
    Face.R: (slice(None), slice(None, None, -1), 2, 2),
    Face.D: (2, slice(None, None, -1), slice(None), 0),
    Face.B: (slice(None), 0, slice(None, None, -1), 1),
    Face.L: (slice(None), slice(None), 0, 2),
}

face_perm = {
    Face.U: [[0, 4, 1, 5], [8, 13, 9, 12]],
    Face.F: [[1, 7, 3, 5], [9, 19, 11, 18]],
    Face.R: [[1, 4, 2, 7], [13, 17, 15, 19]],
    Face.D: [[2, 6, 3, 7], [10, 14, 11, 15]],
    Face.B: [[0, 6, 2, 4], [8, 16, 10, 17]],
    Face.L: [[0, 5, 3, 6], [12, 18, 14, 16]]
}


class Move(IntEnum):
    U1 = 0
    U2 = 1
    U3 = 2
    F1 = 3
    F2 = 4
    F3 = 5
    R1 = 6
    R2 = 7
    R3 = 8
    D1 = 9
    D2 = 10
    D3 = 11
    B1 = 12
    B2 = 13
    B3 = 14
    L1 = 15
    L2 = 16
    L3 = 17

    @property
    def face(self) -> Face:
        return move_face[self]

    @property
    def shift(self) -> int:
        return move_shift[self]

    # @property
    # def str(self) -> str:
    #     str = self.name
    #     if self.name[1] == "1":
    #         str = str[0]
    #     elif self.name[1] == "3":
    #         str = str[0] + "'"
    #     return str


move_face = {move: Face[move.name[0]] for move in Move}
# face_map[Move.NONE] = Face.NONE
move_shift = {move: int(move.name[1]) if move.name[1] != 3 else -1 for move in Move}
# shift_map[Move.NONE] = 0


class Cubie(IntEnum):
    # corners
    UBL = 0
    UFR = 1
    DBR = 2
    DFL = 3
    UBR = 4
    UFL = 5
    DBL = 6
    DFR = 7
    # edges
    UB = 8
    UF = 9
    DB = 10
    DF = 11
    UL = 12
    UR = 13
    DL = 14
    DR = 15
    BL = 16
    BR = 17
    FL = 18
    FR = 19
    # centers
    U = 20
    F = 21
    R = 22
    D = 23
    B = 24
    L = 25
    # core
    CORE = 26

    @classmethod
    def corners(cls) -> Iterator:
        for i in range(8):
            yield cls(i)

    @classmethod
    def edges(cls) -> Iterator:
        for i in range(8, 20):
            yield cls(i)

    @classmethod
    def centers(cls) -> Iterator:
        for i in range(20, 26):
            yield cls(i)

    @property
    def axis(self) -> int:
        return cubie_axis[self]

    @property
    def index(self) -> tuple:
        return cubie_index[self]


# TODO+ [EMPTY]
cubie_axis = {
    # corners
    Cubie.UBL: 0,
    Cubie.UFR: 0,
    Cubie.DBR: 0,
    Cubie.DFL: 0,
    Cubie.UBR: 1,
    Cubie.UFL: 1,
    Cubie.DBL: 1,
    Cubie.DFR: 1,
    # edges
    Cubie.UB: 2,
    Cubie.UF: 2,
    Cubie.DB: 2,
    Cubie.DF: 2,
    Cubie.UL: 1,
    Cubie.UR: 1,
    Cubie.DL: 1,
    Cubie.DR: 1,
    Cubie.BL: 0,
    Cubie.BR: 0,
    Cubie.FL: 0,
    Cubie.FR: 0,
    # centers
    Cubie.U: 0,
    Cubie.F: 1,
    Cubie.R: 2,
    Cubie.D: 0,
    Cubie.B: 1,
    Cubie.L: 2,
    # core
    Cubie.CORE: -1
}


cubie_index = {
    # corners
    Cubie.UBL: (0, 0, 0),
    Cubie.UFR: (0, 2, 2),
    Cubie.DBR: (2, 0, 2),
    Cubie.DFL: (2, 2, 0),
    Cubie.UBR: (0, 0, 2),
    Cubie.UFL: (0, 2, 0),
    Cubie.DBL: (2, 0, 0),
    Cubie.DFR: (2, 2, 2),
    # edges
    Cubie.UB: (0, 0, 1),
    Cubie.UF: (0, 2, 1),
    Cubie.DB: (2, 0, 1),
    Cubie.DF: (2, 2, 1),
    Cubie.UL: (0, 1, 0),
    Cubie.UR: (0, 1, 2),
    Cubie.DL: (2, 1, 0),
    Cubie.DR: (2, 1, 2),
    Cubie.BL: (1, 0, 0),
    Cubie.BR: (1, 0, 2),
    Cubie.FL: (1, 2, 0),
    Cubie.FR: (1, 2, 2),
    # centers
    Cubie.U: (0, 1, 1),
    Cubie.F: (1, 2, 1),
    Cubie.R: (1, 1, 2),
    Cubie.D: (2, 1, 1),
    Cubie.B: (1, 0, 1),
    Cubie.L: (1, 1, 0),
    # core
    Cubie.CORE: (1, 1, 1)
}


if __name__ == "__main__":
    print(list(Cubie.corners()))
