"""Enums module."""
from __future__ import annotations

from enum import Enum, auto
from typing import Union, List, Tuple, Iterator

from .defs import NONE


class IntEnum(int, Enum):
    """Integer enumeration."""
    __repr__ = Enum.__str__


class Axis(IntEnum):
    """Axis enumeration."""
    NONE = NONE  #: No axis.
    X = auto()  #: `X` axis along :attr:`Cubie.R` and :attr:`Cubie.L` centers.
    Y = auto()  #: `Y` axis along :attr:`Cubie.U` and :attr:`Cubie.D` centers.
    Z = auto()  #: `Z` axis along :attr:`Cubie.F` and :attr:`Cubie.B` centers.
    DIAG_111 = auto()  #: Diagonal axis along :attr:`Cubie.UFR` and :attr:`Cubie.DBL` corners.
    DIAG_M11 = auto()  #: Diagonal axis along :attr:`Cubie.UFL` and :attr:`Cubie.DBR` corners.
    DIAG_1M1 = auto()  #: Diagonal axis along :attr:`Cubie.UBL` and :attr:`Cubie.DFR` corners.
    DIAG_11M = auto()  #: Diagonal axis along :attr:`Cubie.UBR` and :attr:`Cubie.DFL` corners.
    EDGE_011 = auto()  #: Edge axis along :attr:`Cubie.UF` and :attr:`Cubie.DB` edges.
    EDGE_0M1 = auto()  #: Edge axis along :attr:`Cubie.UB` and :attr:`Cubie.DF` edges.
    EDGE_101 = auto()  #: Edge axis along :attr:`Cubie.FR` and :attr:`Cubie.BL` edges.
    EDGE_M01 = auto()  #: Edge axis along :attr:`Cubie.FL` and :attr:`Cubie.BR` edges.
    EDGE_110 = auto()  #: Edge axis along :attr:`Cubie.UR` and :attr:`Cubie.DL` edges.
    EDGE_M10 = auto()  #: Edge axis along :attr:`Cubie.UL` and :attr:`Cubie.DR` edges.

    @property
    def is_cartesian(self) -> bool:
        """Whether this is a cartesian axis."""
        return 0 <= self < 3

    @property
    def is_diagonal(self) -> bool:
        """Whether this is a diagonal axis."""
        return 3 <= self < 7

    @property
    def is_edge(self) -> bool:
        """Whether this is an edge axis."""
        return 7 <= self < 13

    @classmethod
    def axes(cls) -> Iterator[Axis]:
        """Iterate over valid axes."""
        for i in range(13):
            yield cls(i)

    @classmethod
    def cartesian_axes(cls) -> Iterator[Axis]:
        """Iterate over cartesian axes."""
        for i in range(3):
            yield cls(i)

    @classmethod
    def diagonal_axes(cls) -> Iterator[Axis]:
        """Iterate over diagonal axes."""
        for i in range(3, 7):
            yield cls(i)

    @classmethod
    def edge_axes(cls) -> Iterator[Axis]:
        """Iterate over edge axes."""
        for i in range(7, 13):
            yield cls(i)


class Orbit(IntEnum):
    """Orbit enumeration."""
    NONE = NONE  #: No orbit.
    SLICE_MIDDLE = auto()  #: `Middle` slice orbit.
    SLICE_EQUATOR = auto()  #: `Equator` slice orbit.
    SLICE_STANDING = auto()  #: `Standing` slice orbit.
    TETRAD_111 = auto()
    """`Tetrad` orbit containing the :attr:`Cubie.UBL`, :attr:`Cubie.UFR`, :attr:`Cubie.DBR`, and :attr:`Cubie.DFL` corners."""
    TETRAD_M11 = auto()
    """`Tetrad` orbit containing the :attr:`Cubie.UBR`, :attr:`Cubie.UFL`, :attr:`Cubie.DBL`, and :attr:`Cubie.DFR` corners."""

    @property
    def is_slice(self) -> bool:
        """Whether this is a slice orbit."""
        return 0 <= self < 3

    @property
    def is_tetrad(self) -> bool:
        """Whether this is a tetrad orbit."""
        return 3 <= self < 5

    @classmethod
    def orbits(cls) -> Iterator[Orbit]:
        """Iterate over valid orbits."""
        for i in range(5):
            yield cls(i)

    @classmethod
    def slices(cls) -> Iterator[Orbit]:
        """Iterate over slice orbits."""
        for i in range(3):
            yield cls(i)

    @classmethod
    def tetrads(cls) -> Iterator[Orbit]:
        """Iterate over tetrad orbits."""
        for i in range(3, 5):
            yield cls(i)


class Layer(IntEnum):
    """Layer enumeration."""
    NONE = NONE  #: No layer.
    UP = auto()  #: `Up` layer.
    FRONT = auto()  #: `Front` layer.
    RIGHT = auto()  #: `Right` layer.
    DOWN = auto()  #: `Down` layer.
    BACK = auto()  #: `Back` layer.
    LEFT = auto()  #: `Left` layer.
    MIDDLE = auto()  #: `Middle` layer.
    EQUATOR = auto()  #: `Equator` layer.
    STANDING = auto()  #: `Standing` layer.

    @property
    def char(self) -> str:
        """Character representation of the layer."""
        return self.name[0]

    @property
    def axis(self) -> Axis:
        """Layer axis."""
        return layer_axis[self]

    @property
    def perm(self) -> List[List[Cubie]]:
        """Layer permutation."""
        return layer_perm[self]

    @property
    def is_outer(self) -> bool:
        """Whether this is an outer layer."""
        return 0 <= self < 6

    @property
    def is_inner(self) -> bool:
        """Whether this is an inner layer."""
        return 6 <= self < 9

    @classmethod
    def from_char(cls, char: str) -> Layer:
        """
        Return the corresponding :class:`Layer` enum.

        Parameters
        ----------
        char : {'N', 'U', 'F', 'R', 'D', 'B', 'L', 'M', 'E', 'S'}
            Character representing the layer.

            * `'N'` means :attr:`Layer.NONE`.
            * `'U'` means :attr:`Layer.UP`.
            * `'F'` means :attr:`Layer.FRONT`.
            * `'R'` means :attr:`Layer.RIGHT`.
            * `'D'` means :attr:`Layer.DOWN`.
            * `'B'` means :attr:`Layer.BACK`.
            * `'L'` means :attr:`Layer.LEFT`.
            * `'M'` means :attr:`Layer.MIDDLE`.
            * `'E'` means :attr:`Layer.EQUATOR`.
            * `'S'` means :attr:`Layer.STANDING`.

        Returns
        -------
        layer : Layer
            :class:`Layer` enum.
        """
        if not isinstance(char, str):
            raise TypeError(f"char must be str, not {type(char).__name__}")
        try:
            return char_layer[char]
        except KeyError:
            raise ValueError(f"invalid face character (got '{char}')")

    @classmethod
    def layers(cls) -> Iterator[Layer]:
        """Iterate over valid layers."""
        for i in range(9):
            yield cls(i)

    @classmethod
    def outers(cls) -> Iterator[Layer]:
        """Iterate over outer layers."""
        for i in range(6):
            yield cls(i)

    @classmethod
    def inners(cls) -> Iterator[Layer]:
        """Iterate over inner layers."""
        for i in range(6, 9):
            yield cls(i)


class Color(IntEnum):
    """Color enumeration."""
    NONE = NONE  #: No color.
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
    def from_char(cls, char: str) -> Color:
        """
        Return the corresponding :class:`Color` enum.

        Parameters
        ----------
        char : {'N', 'W', 'G', 'R', 'Y', 'B', 'O'}
            Character representing the color.

            * `'N'` means :attr:`Color.NONE`.
            * `'W'` means :attr:`Color.WHITE`.
            * `'G'` means :attr:`Color.GREEN`.
            * `'R'` means :attr:`Color.RED`.
            * `'Y'` means :attr:`Color.YELLOW`.
            * `'B'` means :attr:`Color.BLUE`.
            * `'O'` means :attr:`Color.ORANGE`.

        Returns
        -------
        color : Color
            :class:`Color` enum.
        """
        if not isinstance(char, str):
            raise TypeError(f"char must be str, not {type(char).__name__}")
        try:
            return char_color[char]
        except KeyError:
            raise ValueError(f"invalid color character (got '{char}')")

    @classmethod
    def colors(cls) -> Iterator[Color]:
        """Iterate over valid colors."""
        for i in range(6):
            yield cls(i)


class Face(IntEnum):
    """Face enumeration."""
    NONE = NONE  #: No face.
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
        return Layer[self.name].axis

    @property
    def opposite(self) -> Face:
        """Opposite face."""
        return face_opposite[self]

    @property
    def _index(self) -> Tuple[Union[int, slice], ...]:
        """Face index of the color representation array."""
        return face_cubie_index[self]

    @classmethod
    def from_char(cls, char: str) -> Face:
        """
        Return the corresponding :class:`Face` enum.

        Parameters
        ----------
        char : {'N', 'U', 'F', 'R', 'D', 'B', 'L'}
            Character representing the face.

            * `'N'` means :attr:`Face.NONE`.
            * `'U'` means :attr:`Face.UP`.
            * `'F'` means :attr:`Face.FRONT`.
            * `'R'` means :attr:`Face.RIGHT`.
            * `'D'` means :attr:`Face.DOWN`.
            * `'B'` means :attr:`Face.BACK`.
            * `'L'` means :attr:`Face.LEFT`.

        Returns
        -------
        face : Face
            :class:`Face` enum.
        """
        if not isinstance(char, str):
            raise TypeError(f"char must be str, not {type(char).__name__}")
        try:
            return char_face[char]
        except KeyError:
            raise ValueError(f"invalid face character (got '{char}')")

    @classmethod
    def faces(cls) -> Iterator[Face]:
        """Iterate over valid faces."""
        for i in range(6):
            yield cls(i)


class Cubie(IntEnum):
    """Cubie enumeration."""
    NONE = NONE  #: No cubie.
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
        """Cubie axis."""
        if self.is_center:
            return Layer.from_char(self.name).axis
        return cubie_axis[self]

    @property
    def orbit(self) -> Orbit:
        """
        Cubie orbit.

        * :attr:`Orbit.NONE`: :attr:`Cubie.NONE`, :attr:`Cubie.CORE`, :attr:`Cubie.U`, :attr:`Cubie.F`,
          :attr:`Cubie.R`, :attr:`Cubie.D`, :attr:`Cubie.B`, :attr:`Cubie.L`
        * :attr:`Orbit.SLICE_MIDDLE`: :attr:`Cubie.UB`, :attr:`Cubie.UF`, :attr:`Cubie.DB`, :attr:`Cubie.DF`
        * :attr:`Orbit.SLICE_EQUATOR`: :attr:`Cubie.BL`, :attr:`Cubie.BR`, :attr:`Cubie.FL`, :attr:`Cubie.FR`
        * :attr:`Orbit.SLICE_STANDING`: :attr:`Cubie.UL`, :attr:`Cubie.UR`, :attr:`Cubie.DL`, :attr:`Cubie.DR`
        * :attr:`Orbit.TETRAD_111`: :attr:`Cubie.UBL`, :attr:`Cubie.UFR`, :attr:`Cubie.DBR`, :attr:`Cubie.DFL`
        * :attr:`Orbit.TETRAD_M11`: :attr:`Cubie.UBR`, :attr:`Cubie.UFL`, :attr:`Cubie.DBL`, :attr:`Cubie.DFR`
        """
        if self.is_center:
            return Orbit.NONE
        return cubie_orbit[self]

    @property
    def faces(self) -> List[Face]:
        """Cubie visible faces."""
        if self == Cubie.NONE:
            return [Face.NONE]
        if self == Cubie.CORE:
            return []
        faces = [Face.from_char(char) for char in self.name]
        if self.orbit == Orbit.TETRAD_111:
            faces[1:] = faces[2], faces[1]
        return faces

    @property
    def _index(self) -> Tuple[int, ...]:
        """Cubie index of the color representation array."""
        return cubie_index[self]

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
    def from_faces(cls, faces: List[Face]) -> Cubie:
        """
        Return the corresponding :class:`Cubie` enum.

        Parameters
        ----------
        faces : list of Face
            Visible faces of the cubie. If `faces` represent a corner,
            they must be provided in clockwise order to verify that it's a valid corner.

        Returns
        -------
        cubie : Cubie
            :class:`Cubie` enum.
        """
        if not isinstance(faces, list):
            raise TypeError(f"faces must be list, not {type(faces).__name__}")
        if len(faces) > 3:
            raise ValueError(f"faces length must be at most 3 (got {len(faces)})")
        min_faces = []
        for i, face in enumerate(faces):
            if not isinstance(face, Face):
                raise TypeError(f"faces elements must be Face, not {type(faces[0]).__name__}")
            min_faces.append(faces[i:] + faces[:i])
        if min_faces:
            min_faces = min(min_faces)
        try:
            return faces_cubie[tuple(min_faces)]
        except KeyError:
            raise ValueError(f"invalid cubie faces (got {faces})")

    @classmethod
    def cubies(cls) -> Iterator[Cubie]:
        """Iterate over valud cubies."""
        for i in range(27):
            yield cls(i)

    @classmethod
    def corners(cls) -> Iterator[Cubie]:
        """Iterate over corner cubies."""
        for i in range(8):
            yield cls(i)

    @classmethod
    def edges(cls) -> Iterator[Cubie]:
        """Iterate over edge cubies."""
        for i in range(8, 20):
            yield cls(i)

    @classmethod
    def centers(cls) -> Iterator[Cubie]:
        """Iterate over center cubies."""
        for i in range(20, 26):
            yield cls(i)


class Move(IntEnum):
    """Move enumeration."""
    NONE = NONE  #: No move.
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
        if self == Move.NONE:
            return ""
        str = self.name
        if str[0] in "XYZ":
            str = str[0].lower() + str[1:]
        elif str[1] == "W":
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
        return Layer.from_char(self.name[0]).axis

    @property
    def inverse(self) -> Move:
        """Inverse move."""
        if self == Move.NONE:
            return Move.NONE
        return Move[self.name[:-1] + str(-int(self.name[-1]) % 4)]

    @property
    def layers(self) -> List[Layer]:
        """Move layers."""
        if self == Move.NONE:
            return [Layer.NONE]
        if self.is_rotation:
            return [layer for layer in layer_order if layer.axis == self.axis]
        layers = [Layer.from_char(self.name[0])]
        if self.is_wide:
            layers += [layer for layer in Layer.inners() if layer.axis == self.axis]
        return layers

    @property
    def shifts(self) -> List[int]:
        """Permutation shift for each layer."""
        if self == Move.NONE:
            return [0]
        shift = int(self.name[-1]) if self.name[-1] != "3" else -1
        if self.is_rotation:
            return [shift, -shift, shift if self.axis == Axis.Z else -shift]
        shifts = [shift]
        if self.is_wide:
            if self.axis == Axis.Z:
                shift = -shift
            shifts += [-shift if Move[self.axis.name + "1"].layers[0].char == self.name[0] else shift]
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
    def from_string(cls, string: str) -> Move:
        """
        Return the corresponding :class:`Move` enum.

        Parameters
        ----------
        string : str
            String representation of the move.

        Returns
        -------
        move : Move
            :class:`Move` enum.
        """
        if not isinstance(string, str):
            raise TypeError(f"string must be str, not {type(string).__name__}")
        try:
            return str_move[string]
        except KeyError:
            raise ValueError(f"invalid move string (got '{string}')")

    @classmethod
    def moves(cls) -> Iterator[Move]:
        """Iterate over valid moves."""
        for i in range(54):
            yield cls(i)

    @classmethod
    def face_moves(cls) -> Iterator[Move]:
        """Iterate over face moves."""
        for i in range(18):
            yield cls(i)

    @classmethod
    def slice_moves(cls) -> Iterator[Move]:
        """Iterate over slice moves."""
        for i in range(18, 27):
            yield cls(i)

    @classmethod
    def wide_moves(cls) -> Iterator[Move]:
        """Iterate over wide moves."""
        for i in range(27, 45):
            yield cls(i)

    @classmethod
    def rotations(cls) -> Iterator[Move]:
        """Iterate over rotations."""
        for i in range(45, 54):
            yield cls(i)


layer_axis = {
    Layer.NONE: Axis.NONE,
    Layer.UP: Axis.Y,
    Layer.FRONT: Axis.Z,
    Layer.RIGHT: Axis.X,
    Layer.DOWN: Axis.Y,
    Layer.BACK: Axis.Z,
    Layer.LEFT: Axis.X,
    Layer.MIDDLE: Axis.X,
    Layer.EQUATOR: Axis.Y,
    Layer.STANDING: Axis.Z
}

layer_perm = {
    Layer.NONE: [[Cubie.NONE]],
    Layer.UP: [[Cubie.UBL, Cubie.UBR, Cubie.UFR, Cubie.UFL], [Cubie.UB, Cubie.UR, Cubie.UF, Cubie.UL]],
    Layer.FRONT: [[Cubie.UFR, Cubie.DFR, Cubie.DFL, Cubie.UFL], [Cubie.UF, Cubie.FR, Cubie.DF, Cubie.FL]],
    Layer.RIGHT: [[Cubie.UFR, Cubie.UBR, Cubie.DBR, Cubie.DFR], [Cubie.UR, Cubie.BR, Cubie.DR, Cubie.FR]],
    Layer.DOWN: [[Cubie.DBR, Cubie.DBL, Cubie.DFL, Cubie.DFR], [Cubie.DB, Cubie.DL, Cubie.DF, Cubie.DR]],
    Layer.BACK: [[Cubie.UBL, Cubie.DBL, Cubie.DBR, Cubie.UBR], [Cubie.UB, Cubie.BL, Cubie.DB, Cubie.BR]],
    Layer.LEFT: [[Cubie.UBL, Cubie.UFL, Cubie.DFL, Cubie.DBL], [Cubie.UL, Cubie.FL, Cubie.DL, Cubie.BL]],
    Layer.MIDDLE: [[Cubie.U, Cubie.F, Cubie.D, Cubie.B], [Cubie.UB, Cubie.UF, Cubie.DF, Cubie.DB]],
    Layer.EQUATOR: [[Cubie.F, Cubie.R, Cubie.B, Cubie.L], [Cubie.BL, Cubie.FL, Cubie.FR, Cubie.BR]],
    Layer.STANDING: [[Cubie.U, Cubie.R, Cubie.D, Cubie.L], [Cubie.UL, Cubie.UR, Cubie.DR, Cubie.DL]]
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

face_cubie_index = {
    Face.NONE: (slice(None), slice(None), slice(None)),
    Face.UP: (0, slice(None), slice(None)),
    Face.FRONT: (slice(None), 2, slice(None)),
    Face.RIGHT: (slice(None), slice(None, None, -1), 2),
    Face.DOWN: (2, slice(None, None, -1), slice(None)),
    Face.BACK: (slice(None), 0, slice(None, None, -1)),
    Face.LEFT: (slice(None), slice(None), 0)
}

cubie_axis = {
    Cubie.NONE: Axis.NONE,
    # corners
    Cubie.UBL: Axis.DIAG_1M1,
    Cubie.UFR: Axis.DIAG_111,
    Cubie.DBR: Axis.DIAG_M11,
    Cubie.DFL: Axis.DIAG_11M,
    Cubie.UBR: Axis.DIAG_11M,
    Cubie.UFL: Axis.DIAG_M11,
    Cubie.DBL: Axis.DIAG_111,
    Cubie.DFR: Axis.DIAG_1M1,
    # edges
    Cubie.UB: Axis.EDGE_0M1,
    Cubie.UF: Axis.EDGE_011,
    Cubie.DB: Axis.EDGE_011,
    Cubie.DF: Axis.EDGE_0M1,
    Cubie.UL: Axis.EDGE_M10,
    Cubie.UR: Axis.EDGE_110,
    Cubie.DL: Axis.EDGE_110,
    Cubie.DR: Axis.EDGE_M10,
    Cubie.BL: Axis.EDGE_101,
    Cubie.BR: Axis.EDGE_M01,
    Cubie.FL: Axis.EDGE_M01,
    Cubie.FR: Axis.EDGE_101,
    # core
    Cubie.CORE: Axis.NONE
}

cubie_orbit = {
    Cubie.NONE: Orbit.NONE,
    # corners
    Cubie.UBL: Orbit.TETRAD_111,
    Cubie.UFR: Orbit.TETRAD_111,
    Cubie.DBR: Orbit.TETRAD_111,
    Cubie.DFL: Orbit.TETRAD_111,
    Cubie.UBR: Orbit.TETRAD_M11,
    Cubie.UFL: Orbit.TETRAD_M11,
    Cubie.DBL: Orbit.TETRAD_M11,
    Cubie.DFR: Orbit.TETRAD_M11,
    # edges
    Cubie.UB: Orbit.SLICE_MIDDLE,
    Cubie.UF: Orbit.SLICE_MIDDLE,
    Cubie.DB: Orbit.SLICE_MIDDLE,
    Cubie.DF: Orbit.SLICE_MIDDLE,
    Cubie.UL: Orbit.SLICE_STANDING,
    Cubie.UR: Orbit.SLICE_STANDING,
    Cubie.DL: Orbit.SLICE_STANDING,
    Cubie.DR: Orbit.SLICE_STANDING,
    Cubie.BL: Orbit.SLICE_EQUATOR,
    Cubie.BR: Orbit.SLICE_EQUATOR,
    Cubie.FL: Orbit.SLICE_EQUATOR,
    Cubie.FR: Orbit.SLICE_EQUATOR,
    # core
    Cubie.CORE: Orbit.NONE
}

cubie_index = {
    Cubie.NONE: (1, 1, 1),
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

char_layer = {layer.char: layer for layer in Layer}
char_color = {color.char: color for color in Color}
char_face = {face.char: face for face in Face}
faces_cubie = {tuple(min([cb.faces[i:] + cb.faces[:i] for i in range(len(cb.faces))])if cb.faces else []): cb for cb in Cubie}
layer_order = [Layer.UP, Layer.FRONT, Layer.RIGHT, Layer.DOWN, Layer.BACK, Layer.LEFT] + [*Layer.inners()]
str_move = {move.string: move for move in Move}
str_move.update({move.string[0].lower() + move.string[2:]: move for move in Move.wide_moves()})
