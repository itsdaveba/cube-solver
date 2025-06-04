"""Enums module."""
from typing import Iterator
from enum import Enum, auto


class IntEnum(int, Enum):
    """Integer enumeration."""
    def __repr__(self) -> str:
        return self.name


class Color(IntEnum):  # up, front, right, down, back, left
    """Color enumeration."""
    NONE = -1  #: No color.
    WHITE = auto()  #: White color (Up face).
    GREEN = auto()  #: Green color (Front face).
    RED = auto()  #: Red color (Right face).
    YELLOW = auto()  #: Yellow color (Down face).
    BLUE = auto()  #: Blue color (Back face).
    ORANGE = auto()  #: Orange color (Left face).

    @classmethod
    def colors(cls) -> Iterator:
        """Iterate over valid colors."""
        for i in range(6):
            yield cls(i)

    def __str__(self) -> str:
        return self.name[0]


class Face(IntEnum):
    """Face enumeration."""
    NONE = -1  #: No face.
    UP = auto()  #: Up face.
    FRONT = auto()  #: Front face.
    RIGHT = auto()  #: Right face.
    DOWN = auto()  #: Down face.
    BACK = auto()  #: Back face.
    LEFT = auto()  #: Left face.

    @classmethod
    def faces(cls) -> Iterator:
        """Iterate over valid faces."""
        for i in range(6):
            yield cls(i)


cubie_slice = {
    Face.NONE: (),
    Face.UP: (0, slice(None), slice(None), 0),
    Face.FRONT: (slice(None), 2, slice(None), 1),
    Face.RIGHT: (slice(None), slice(None, None, -1), 2, 2),
    Face.DOWN: (2, slice(None, None, -1), slice(None), 0),
    Face.BACK: (slice(None), 0, slice(None, None, -1), 1),
    Face.LEFT: (slice(None), slice(None), 0, 2),
}

face_perm = {
    Face.NONE: [],
    Face.UP: [[0, 4, 1, 5], [8, 13, 9, 12]],
    Face.FRONT: [[1, 7, 3, 5], [9, 19, 11, 18]],
    Face.RIGHT: [[1, 4, 2, 7], [13, 17, 15, 19]],
    Face.DOWN: [[2, 6, 3, 7], [10, 14, 11, 15]],
    Face.BACK: [[0, 6, 2, 4], [8, 16, 10, 17]],
    Face.LEFT: [[0, 5, 3, 6], [12, 18, 14, 16]]
}


class Move(IntEnum):
    """Move enumeration."""
    NONE = -1  #: No move.
    # face moves
    U1 = auto()  #: `U` face move.
    U2 = auto()  #: `U2` face move.
    U3 = auto()  #: `U'` face move.
    F1 = auto()  #: `F` face move.
    F2 = auto()  #: `F2` face move.
    F3 = auto()  #: `F'` face move.
    R1 = auto()  #: `R` face move.
    R2 = auto()  #: `R2` face move.
    R3 = auto()  #: `R'` face move.
    D1 = auto()  #: `D` face move.
    D2 = auto()  #: `D2` face move.
    D3 = auto()  #: `D'` face move.
    B1 = auto()  #: `B` face move.
    B2 = auto()  #: `B2` face move.
    B3 = auto()  #: `B'` face move.
    L1 = auto()  #: `L` face move.
    L2 = auto()  #: `L2` face move.
    L3 = auto()  #: `L'` face move.
    # slice moves
    M1 = auto()  #: `M` slice move.
    M2 = auto()  #: `M2` slice move.
    M3 = auto()  #: `M'` slice move.
    E1 = auto()  #: `E` slice move.
    E2 = auto()  #: `E2` slice move.
    E3 = auto()  #: `E'` slice move.
    S1 = auto()  #: `S` slice move.
    S2 = auto()  #: `S2` slice move.
    S3 = auto()  #: `S'` slice move.
    # wide moves
    UW1 = auto()  #: `Uw` or `u` wide move.
    UW2 = auto()  #: `Uw2` or `u2` wide move.
    UW3 = auto()  #: `Uw'` or `u'` wide move.
    FW1 = auto()  #: `Fw` or `f` wide move.
    FW2 = auto()  #: `Fw2` or `f2` wide move.
    FW3 = auto()  #: `Fw'` or `f'` wide move.
    RW1 = auto()  #: `Rw` or `r` wide move.
    RW2 = auto()  #: `Rw2` or `r2` wide move.
    RW3 = auto()  #: `Rw'` or `r'` wide move.
    DW1 = auto()  #: `Dw` or `d` wide move.
    DW2 = auto()  #: `Dw2` or `d2` wide move.
    DW3 = auto()  #: `Dw'` or `d'` wide move.
    BW1 = auto()  #: `Bw` or `b` wide move.
    BW2 = auto()  #: `Bw2` or `b2` wide move.
    BW3 = auto()  #: `Bw'` or `b'` wide move.
    LW1 = auto()  #: `Lw` or `l` wide move.
    LW2 = auto()  #: `Lw2` or `l2` wide move.
    LW3 = auto()  #: `Lw'` or `l'` wide move.
    # cube rotations
    X1 = auto()  #: `x` rotation.
    X2 = auto()  #: `x2` rotation.
    X3 = auto()  #: `x'` rotation.
    Y1 = auto()  #: `y` rotation.
    Y2 = auto()  #: `y2` rotation.
    Y3 = auto()  #: `y'` rotation.
    Z1 = auto()  #: `z` rotation.
    Z2 = auto()  #: `z2` rotation.
    Z3 = auto()  #: `z'` rotation.

    @classmethod
    def face_moves(cls) -> Iterator:
        """Iterate over face moves."""
        for i in range(18):
            yield cls(i)

    @classmethod
    def slice_moves(cls) -> Iterator:
        """Iterate over slice moves."""
        for i in range(18, 27):
            yield cls(i)

    @classmethod
    def wide_moves(cls) -> Iterator:
        """Iterate over wide moves."""
        for i in range(27, 45):
            yield cls(i)

    @classmethod
    def rotations(cls) -> Iterator:
        """Iterate over rotations."""
        for i in range(45, 54):
            yield cls(i)


class Cubie(IntEnum):
    """Cubie enumeration."""
    NONE = -1  #: No cubie.
    # corners
    UBL = auto()  #: Up-Back-Left corner.
    UFR = auto()  #: Up-Front-Right corner.
    DBR = auto()  #: Down-Back-Right corner.
    DFL = auto()  #: Down-Front-Left corner.
    UBR = auto()  #: Up-Back-Right corner.
    UFL = auto()  #: Up-Front-Left corner.
    DBL = auto()  #: Down-Back-Left corner.
    DFR = auto()  #: Down-Front-Right corner.
    # edges
    UB = auto()  #: Up-Back edge.
    UF = auto()  #: Up-Front edge.
    DB = auto()  #: Down-Back edge.
    DF = auto()  #: Down-Front edge.
    UL = auto()  #: Up-Left edge.
    UR = auto()  #: Up-Right edge.
    DL = auto()  #: Down-Left edge.
    DR = auto()  #: Down-Right edge.
    BL = auto()  #: Back-Left edge.
    BR = auto()  #: Back-Right edge.
    FL = auto()  #: Front-Left edge.
    FR = auto()  #: Front-Right edge.
    # centers
    U = auto()  #: Up center.
    F = auto()  #: Front center.
    R = auto()  #: Right center.
    D = auto()  #: Down center.
    B = auto()  #: Back center.
    L = auto()  #: Left center.
    # core
    CORE = auto()  #: Core.

    @classmethod
    def corners(cls) -> Iterator:
        """Iterate over corner cubies."""
        for i in range(8):
            yield cls(i)

    @classmethod
    def edges(cls) -> Iterator:
        """Iterate over edge cubies."""
        for i in range(8, 20):
            yield cls(i)

    @classmethod
    def centers(cls) -> Iterator:
        """Iterate over center cubies."""
        for i in range(20, 26):
            yield cls(i)

    @property
    def axis(self) -> int:  # TODO make private
        """
        Cubie axis.

        For edges and centers:

        * `0` along `Up` and `Down` faces.
        * `1` along `Front` and `Back` faces.
        * `2` along `Right` and `Left` faces.

        For corners:

        * `0` along `Up-Front-Left` and `Down-Back-Right` corners.
        * `1` along `Up-Front-Right` and `Down-Back-Left` corners.
        """
        return cubie_axis[self]

    @property
    def index(self) -> tuple:  # TODO make private
        """
        Cubie index.

        Used as the index of the cubie representation array.
        """
        return cubie_index[self]


cubie_axis = {
    Cubie.NONE: -1,
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
    Cubie.NONE: (),
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
