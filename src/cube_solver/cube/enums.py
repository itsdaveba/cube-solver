"""Enums module."""
from typing import Iterator
from enum import Enum, auto


class IntEnum(int, Enum):
    """Integer enumeration."""
    __repr__ = Enum.__str__


class Axis(IntEnum):
    """Axis enumeration."""
    NONE = -1  #: no axis.
    Y = auto()  #: `Y` axis along :attr:`Face.UP` and :attr:`Face.DOWN`.
    Z = auto()  #: `Z` axis along :attr:`Face.FRONT` and :attr:`Face.BACK`.
    X = auto()  #: `X` axis along :attr:`Face.RIGHT` and :attr:`Face.LEFT`.
    DIAG_M11 = auto()  #: Diagonal axis along :attr:`Cubie.UFL` and :attr:`Cubie.DBR` corners.
    DIAG_111 = auto()  #: Diagonal axis along :attr:`Cubie.UFR` and :attr:`Cubie.DBL` corners.


class Layer(IntEnum):
    """Layer enumeration."""
    NONE = -1  #: no layer.
    U = auto()  #: `U` layer.
    F = auto()  #: `F` layer.
    R = auto()  #: `R` layer.
    D = auto()  #: `D` layer.
    B = auto()  #: `B` layer.
    L = auto()  #: `L` layer.
    M = auto()  #: `M` layer.
    E = auto()  #: `E` layer.
    S = auto()  #: `S` layer.

    @property
    def axis(self) -> Axis:
        """Layer axis."""
        return layer_axis[self]

    @property
    def perm(self) -> list[list["Cubie"]]:
        """Layer permutation."""
        return layer_perm[self]


class Color(IntEnum):
    """Color enumeration."""
    NONE = -1  #: No color.
    WHITE = auto()  #: `White` color.
    GREEN = auto()  #: `Green` color.
    RED = auto()  #: `Red` color.
    YELLOW = auto()  #: `Yellow` color.
    BLUE = auto()  #: `Blue` color.
    ORANGE = auto()  #: `Orange` color.

    @property
    def char(self) -> str:
        """Character representation of the color."""
        return self.name[0]

    @classmethod
    def from_char(cls, char: str) -> "Color":
        """
        Return the corresponding :class:`Color` enum.

        Parameters
        ----------
        char : {'N', 'W', 'G', 'R', 'Y', 'B', 'O'}
            Character representing the color.

            * 'N' means :attr:`Color.NONE`.
            * 'W' means :attr:`Color.WHITE`.
            * 'G' means :attr:`Color.GREEN`.
            * 'R' means :attr:`Color.RED`.
            * 'Y' means :attr:`Color.YELLOW`.
            * 'B' means :attr:`Color.BLUE`.
            * 'O' means :attr:`Color.ORANGE`.
        """
        if not isinstance(char, str):
            raise TypeError(f"char must be str, not {type(char).__name__}")
        try:
            return char_color[char]
        except KeyError:
            raise ValueError(f"invalid color character, got '{char}'")

    @classmethod
    def colors(cls) -> Iterator["Color"]:
        """Iterate over valid colors."""
        for i in range(6):
            yield cls(i)


class Face(IntEnum):
    """Face enumeration."""
    NONE = -1  #: No face.
    UP = auto()  #: `Up` face.
    FRONT = auto()  #: `Front` face.
    RIGHT = auto()  #: `Right` face.
    DOWN = auto()  #: `Down` face.
    BACK = auto()  #: `Back` face.
    LEFT = auto()  #: `Left` face.

    @property
    def char(self) -> str:
        """Character representation of the face."""
        return self.name[0]

    @property
    def axis(self) -> Axis:
        """Face axis."""
        if self == Face.NONE:
            return Axis.NONE
        return Layer[self.char].axis

    @property
    def opposite(self) -> "Face":
        """Opposite face."""
        return face_opposite[self]

    @property
    def _cubie_slice(self) -> tuple[int | slice]:
        """Face slice for the cubie representation."""
        return face_cubie_slice[self]

    @classmethod
    def from_char(cls, char: str) -> "Face":
        """
        Return the corresponding :class:`Face` enum.

        Parameters
        ----------
        char : {'N', 'U', 'F', 'R', 'D', 'B', 'L'}
            Character representing the color.

            * 'N' means :attr:`Face.NONE`.
            * 'U' means :attr:`Face.UP`.
            * 'F' means :attr:`Face.FRONT`.
            * 'R' means :attr:`Face.RIGHT`.
            * 'D' means :attr:`Face.DOWN`.
            * 'B' means :attr:`Face.BACK`.
            * 'L' means :attr:`Face.LEFT`.
        """
        if not isinstance(char, str):
            raise TypeError(f"char must be str, not {type(char).__name__}")
        try:
            return char_face[char]
        except KeyError:
            raise ValueError(f"invalid face character, got '{char}'")

    @classmethod
    def faces(cls) -> Iterator["Face"]:
        """Iterate over valid faces."""
        for i in range(6):
            yield cls(i)


class Cubie(IntEnum):
    """Cubie enumeration."""
    NONE = -1  #: No cubie.
    # corners
    UBL = auto()  #: `Up-Back-Left` corner.
    UFR = auto()  #: `Up-Front-Right` corner.
    DBR = auto()  #: `Down-Back-Right` corner.
    DFL = auto()  #: `Down-Front-Left` corner.
    UBR = auto()  #: `Up-Back-Right` corner.
    UFL = auto()  #: `Up-Front-Left` corner.
    DBL = auto()  #: `Down-Back-Left` corner.
    DFR = auto()  #: `Down-Front-Right` corner.
    # edges
    UB = auto()  #: `Up-Back` edge.
    UF = auto()  #: `Up-Front` edge.
    DB = auto()  #: `Down-Back` edge.
    DF = auto()  #: `Down-Front` edge.
    UL = auto()  #: `Up-Left` edge.
    UR = auto()  #: `Up-Right` edge.
    DL = auto()  #: `Down-Left` edge.
    DR = auto()  #: `Down-Right` edge.
    BL = auto()  #: `Back-Left` edge.
    BR = auto()  #: `Back-Right` edge.
    FL = auto()  #: `Front-Left` edge.
    FR = auto()  #: `Front-Right` edge.
    # centers
    U = auto()  #: `Up` center.
    F = auto()  #: `Front` center.
    R = auto()  #: `Right` center.
    D = auto()  #: `Down` center.
    B = auto()  #: `Back` center.
    L = auto()  #: `Left` center.
    # core
    CORE = auto()  #: `Core`.

    @property
    def axis(self) -> Axis:
        """
        Cubie axis.

        * :attr:`Axis.NONE`: :attr:`Cubie.NONE`, :attr:`Cubie.CORE`
        * :attr:`Axis.X`: :attr:`Cubie.R`, :attr:`Cubie.L`, :attr:`Cubie.UB`,
          :attr:`Cubie.UF`, :attr:`Cubie.DB`, :attr:`Cubie.DF`
        * :attr:`Axis.Y`: :attr:`Cubie.U`, :attr:`Cubie.D`, :attr:`Cubie.BL`,
          :attr:`Cubie.BR`, :attr:`Cubie.FL`, :attr:`Cubie.FR`
        * :attr:`Axis.Z`: :attr:`Cubie.F`, :attr:`Cubie.B`, :attr:`Cubie.UL`,
          :attr:`Cubie.UR`, :attr:`Cubie.DL`, :attr:`Cubie.DR`
        * :attr:`Axis.DIAG_111`: :attr:`Cubie.UBR`, :attr:`Cubie.UFL`, :attr:`Cubie.DBL`, :attr:`Cubie.DFR`
        * :attr:`Axis.DIAG_M11`: :attr:`Cubie.UBL`, :attr:`Cubie.UFR`, :attr:`Cubie.DBR`, :attr:`Cubie.DFL`
        """
        return cubie_axis[self]

    @property
    def _index(self) -> tuple[int, int, int]:
        """
        Cubie index.

        Used as the index of the cubie representation array.
        """
        return cubie_index[self]

    @property
    def faces(self) -> list[Face]:
        """Cubie visible faces."""
        if self == Cubie.NONE:
            return [Face.NONE]
        if self == Cubie.CORE:
            return []
        faces = [Face.from_char(char) for char in self.name]
        if self.axis == Axis.DIAG_M11:
            faces[1:] = faces[2], faces[1]
        return faces

    @property
    def is_corner(self) -> bool:
        """Whether the cubie is a corner."""
        return 0 <= self < 8

    @property
    def is_edge(self) -> bool:
        """Whether the cubie is an edge."""
        return 8 <= self < 20

    @property
    def is_center(self) -> bool:
        """Whether the cubie is a center."""
        return 20 <= self < 26

    @classmethod
    def from_faces(cls, faces: list[Face]) -> "Cubie":
        """
        Return the corresponding :class:`Cubie` enum.

        Parameters
        ----------
        faces : list[Face]
            Visible faces of the cubie. If `faces` represent a corner,
            they must be provided in clockwise order to verify that it's a valid corner.
        """
        if not isinstance(faces, list):
            raise TypeError(f"faces must be list, not {type(faces).__name__}")
        if len(faces) > 3:
            raise ValueError(f"faces length must be at most 3, got {len(faces)}")
        for face in faces:
            if not isinstance(face, Face):
                raise TypeError(f"faces elements must be Face, not {type(faces).__name__}")
        min_faces = min([faces[i:] + faces[:i] for i in range(len(faces))]) if faces else []
        try:
            return faces_cubie[tuple(min_faces)]
        except KeyError:
            raise ValueError(f"invalid cubie faces, got {faces}")

    @classmethod
    def corners(cls) -> Iterator["Cubie"]:
        """Iterate over corner cubies."""
        for i in range(8):
            yield cls(i)

    @classmethod
    def edges(cls) -> Iterator["Cubie"]:
        """Iterate over edge cubies."""
        for i in range(8, 20):
            yield cls(i)

    @classmethod
    def centers(cls) -> Iterator["Cubie"]:
        """Iterate over center cubies."""
        for i in range(20, 26):
            yield cls(i)


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

    @property
    def string(self) -> str:
        """String representation of the move."""
        str = self.name
        if str[0] in "XYZ":
            str = str[0].lower() + str[1:]
        if str[1] == "W":
            str = self.name[0] + "w" + self.name[2:]
        if str[-1] == "1":
            str = str[:-1]
        elif str[-1] == "3":
            str = str[:-1] + "'"
        return str

    @property
    def axis(self) -> Axis:
        """Move axis."""
        if self == Move.NONE:
            return Axis.NONE
        if self.is_rotation:
            return Axis[self.name[0]]
        return Layer[self.name[0]].axis

    @property
    def layers(self) -> list[Layer]:
        """Move layers."""
        if self == Move.NONE:
            return [Layer.NONE]
        if self.is_rotation:
            if self.axis == Axis.X:
                return [Layer.R, Layer.M, Layer.L]
            if self.axis == Axis.Y:
                return [Layer.U, Layer.E, Layer.D]
            if self.axis == Axis.Z:
                return [Layer.F, Layer.S, Layer.B]
        layers = [Layer[self.name[0]]]
        if self.is_wide:
            if self.axis == Axis.X:
                layers += [Layer.M]
            if self.axis == Axis.Y:
                layers += [Layer.E]
            if self.axis == Axis.Z:
                layers += [Layer.S]
        return layers

    @property
    def shifts(self) -> list[int]:
        """Permutation shift for each layer."""
        if self == Move.NONE:
            return [0]
        shift = int(self.name[-1]) if self.name[-1] != "3" else -1
        if self.is_rotation:
            return [shift, shift if self.axis == Axis.Z else -shift, -shift]
        shifts = [shift]
        if self.is_wide:
            if self.axis == Axis.Z:
                shift = -shift
            shifts += [-shift if Move[self.axis.name + "1"].layers[0].name == self.name[0] else shift]
        return shifts

    @property
    def is_face(self) -> bool:
        """Whether this is a face move."""
        return 0 <= self < 18

    @property
    def is_slice(self) -> bool:
        """Whether this is a slice move."""
        return 18 <= self < 27

    @property
    def is_wide(self) -> bool:
        """Whether this is a wide move."""
        return 27 <= self < 45

    @property
    def is_rotation(self) -> bool:
        """Whether the move is a rotation."""
        return 45 <= self < 54

    @classmethod
    def from_str(cls, string: str) -> "Move":
        """
        Return the corresponding :class:`Move` enum.

        Parameters
        ----------
        string : str
            String representation of the move.
        """
        if not isinstance(string, str):
            raise TypeError(f"string must be str, not {type(string).__name__}")
        try:
            return str_move[string]
        except KeyError:
            raise ValueError(f"invalid move string, got '{string}'")

    @classmethod
    def face_moves(cls) -> Iterator["Move"]:
        """Iterate over face moves."""
        for i in range(18):
            yield cls(i)

    @classmethod
    def slice_moves(cls) -> Iterator["Move"]:
        """Iterate over slice moves."""
        for i in range(18, 27):
            yield cls(i)

    @classmethod
    def wide_moves(cls) -> Iterator["Move"]:
        """Iterate over wide moves."""
        for i in range(27, 45):
            yield cls(i)

    @classmethod
    def rotations(cls) -> Iterator["Move"]:
        """Iterate over rotations."""
        for i in range(45, 54):
            yield cls(i)


layer_axis = {
    Layer.NONE: Axis.NONE,
    Layer.U: Axis.Y,
    Layer.F: Axis.Z,
    Layer.R: Axis.X,
    Layer.D: Axis.Y,
    Layer.B: Axis.Z,
    Layer.L: Axis.X,
    Layer.M: Axis.X,
    Layer.E: Axis.Y,
    Layer.S: Axis.Z
}

layer_perm = {
    Layer.NONE: [[Cubie.NONE]],
    Layer.U: [[Cubie.UBL, Cubie.UBR, Cubie.UFR, Cubie.UFL], [Cubie.UB, Cubie.UR, Cubie.UF, Cubie.UL]],
    Layer.F: [[Cubie.UFR, Cubie.DFR, Cubie.DFL, Cubie.UFL], [Cubie.UF, Cubie.FR, Cubie.DF, Cubie.FL]],
    Layer.R: [[Cubie.UFR, Cubie.UBR, Cubie.DBR, Cubie.DFR], [Cubie.UR, Cubie.BR, Cubie.DR, Cubie.FR]],
    Layer.D: [[Cubie.DBR, Cubie.DBL, Cubie.DFL, Cubie.DFR], [Cubie.DB, Cubie.DL, Cubie.DF, Cubie.DR]],
    Layer.B: [[Cubie.UBL, Cubie.DBL, Cubie.DBR, Cubie.UBR], [Cubie.UB, Cubie.BL, Cubie.DB, Cubie.BR]],
    Layer.L: [[Cubie.UBL, Cubie.UFL, Cubie.DFL, Cubie.DBL], [Cubie.UL, Cubie.FL, Cubie.DL, Cubie.BL]],
    Layer.M: [[Cubie.U, Cubie.F, Cubie.D, Cubie.B], [Cubie.UB, Cubie.UF, Cubie.DF, Cubie.DB]],
    Layer.E: [[Cubie.F, Cubie.R, Cubie.B, Cubie.L], [Cubie.BL, Cubie.FL, Cubie.FR, Cubie.BR]],
    Layer.S: [[Cubie.U, Cubie.R, Cubie.D, Cubie.L], [Cubie.UL, Cubie.UR, Cubie.DR, Cubie.DL]]
}

face_opposite = {
    Face.NONE: Face.NONE,
    Face.UP: Face.DOWN,
    Face.FRONT: Face.BACK,
    Face.RIGHT: Face.LEFT,
    Face.DOWN: Face.UP,
    Face.BACK: Face.FRONT,
    Face.LEFT: Face.RIGHT
}

face_cubie_slice = {
    Face.NONE: (slice(None), slice(None), slice(None), slice(None)),
    Face.UP: (0, slice(None), slice(None), Face.UP.axis),
    Face.FRONT: (slice(None), 2, slice(None), Face.FRONT.axis),
    Face.RIGHT: (slice(None), slice(None, None, -1), 2, Face.RIGHT.axis),
    Face.DOWN: (2, slice(None, None, -1), slice(None), Face.DOWN.axis),
    Face.BACK: (slice(None), 0, slice(None, None, -1), Face.BACK.axis),
    Face.LEFT: (slice(None), slice(None), 0, Face.LEFT.axis)
}

cubie_axis = {
    Cubie.NONE: Axis.NONE,
    # corners
    Cubie.UBL: Axis.DIAG_M11,
    Cubie.UFR: Axis.DIAG_M11,
    Cubie.DBR: Axis.DIAG_M11,
    Cubie.DFL: Axis.DIAG_M11,
    Cubie.UBR: Axis.DIAG_111,
    Cubie.UFL: Axis.DIAG_111,
    Cubie.DBL: Axis.DIAG_111,
    Cubie.DFR: Axis.DIAG_111,
    # edges
    Cubie.UB: Axis.X,
    Cubie.UF: Axis.X,
    Cubie.DB: Axis.X,
    Cubie.DF: Axis.X,
    Cubie.UL: Axis.Z,
    Cubie.UR: Axis.Z,
    Cubie.DL: Axis.Z,
    Cubie.DR: Axis.Z,
    Cubie.BL: Axis.Y,
    Cubie.BR: Axis.Y,
    Cubie.FL: Axis.Y,
    Cubie.FR: Axis.Y,
    # centers
    Cubie.U: Axis.Y,
    Cubie.F: Axis.Z,
    Cubie.R: Axis.X,
    Cubie.D: Axis.Y,
    Cubie.B: Axis.Z,
    Cubie.L: Axis.X,
    # core
    Cubie.CORE: Axis.NONE
}

cubie_index = {
    Cubie.NONE: (slice(None), slice(None), slice(None)),
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

char_color = {color.char: color for color in Color}
char_face = {face.char: face for face in Face}
faces_cubie = {tuple(min([cb.faces[i:] + cb.faces[:i] for i in range(len(cb.faces))])if cb.faces else []): cb for cb in Cubie}
str_move = {move.string: move for move in Move}
str_move.update({move.string[0].lower() + move.string[2:]: move for move in Move.wide_moves()})
