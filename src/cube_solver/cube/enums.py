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
    X = auto()  #: `X` axis along :attr:`Face.RIGHT` and :attr:`Face.LEFT` centers.
    Y = auto()  #: `Y` axis along :attr:`Face.UP` and :attr:`Face.DOWN` centers.
    Z = auto()  #: `Z` axis along :attr:`Face.FRONT` and :attr:`Face.BACK` centers.

    @classmethod
    def axes(cls) -> Iterator[Axis]:
        """Iterate over valid axes."""
        for i in range(3):
            yield cls(i)


class Orbit(IntEnum):
    """Orbit enumeration."""
    NONE = NONE  #: No orbit.
    TETRAD_111 = auto()
    """`Tetrad` orbit containing the :attr:`Cubie.UBL`, :attr:`Cubie.UFR`, :attr:`Cubie.DFL`, and :attr:`Cubie.DBR` corners."""
    TETRAD_M11 = auto()
    """`Tetrad` orbit containing the :attr:`Cubie.UBR`, :attr:`Cubie.UFL`, :attr:`Cubie.DFR`, and :attr:`Cubie.DBL` corners."""

    @classmethod
    def orbits(cls) -> Iterator[Orbit]:
        """Iterate over valid orbits."""
        for i in range(2):
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

    @property
    def char(self) -> str:
        """Character representation of the layer."""
        return self.name[0]

    @property
    def axis(self) -> Axis:
        """Layer axis."""
        return layer_axis[self]

    @property
    def perm(self) -> List[Cubie]:
        """Layer permutation."""
        return layer_perm[self]

    @classmethod
    def from_char(cls, char: str) -> Layer:
        """
        Return the corresponding :class:`Layer` enum.

        Parameters
        ----------
        char : {'N', 'U', 'F', 'R', 'D', 'B', 'L'}
            Character representing the layer.

            * `'N'` means :attr:`Layer.NONE`.
            * `'U'` means :attr:`Layer.UP`.
            * `'F'` means :attr:`Layer.FRONT`.
            * `'R'` means :attr:`Layer.RIGHT`.
            * `'D'` means :attr:`Layer.DOWN`.
            * `'B'` means :attr:`Layer.BACK`.
            * `'L'` means :attr:`Layer.LEFT`.

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
        for i in range(6):
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
    DFL = auto()  #: `Down-Front-Left` corner.
    DBR = auto()  #: `Down-Back-Right` corner.
    UBR = auto()  #: `Up-Back-Right` corner.
    UFL = auto()  #: `Up-Front-Left` corner.
    DFR = auto()  #: `Down-Front-Right` corner.
    DBL = auto()  #: `Down-Back-Left` corner.

    @property
    def orbit(self) -> Orbit:
        """
        Cubie orbit.

        * :attr:`Orbit.NONE`: :attr:`Cubie.NONE`
        * :attr:`Orbit.TETRAD_111`: :attr:`Cubie.UBL`, :attr:`Cubie.UFR`, :attr:`Cubie.DFL`, :attr:`Cubie.DBR`
        * :attr:`Orbit.TETRAD_M11`: :attr:`Cubie.UBR`, :attr:`Cubie.UFL`, :attr:`Cubie.DFR`, :attr:`Cubie.DBL`
        """
        return cubie_orbit[self]

    @property
    def faces(self) -> List[Face]:
        """Cubie visible faces."""
        if self == Cubie.NONE:
            return [Face.NONE]
        faces = [Face.from_char(char) for char in self.name]
        if self.orbit == Orbit.TETRAD_111:
            faces[1:] = faces[2], faces[1]
        return faces

    @property
    def _index(self) -> Tuple[int, ...]:
        """Cubie index of the color representation array."""
        return cubie_index[self]

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
        for i in range(8):
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
        return shifts

    @property
    def is_face(self) -> bool:
        """Whether this is a face move."""
        return 0 <= self < 18

    @property
    def is_rotation(self) -> bool:
        """Whether the move is a rotation."""
        return 18 <= self < 27

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
        for i in range(27):
            yield cls(i)

    @classmethod
    def face_moves(cls) -> Iterator[Move]:
        """Iterate over face moves."""
        for i in range(18):
            yield cls(i)

    @classmethod
    def rotations(cls) -> Iterator[Move]:
        """Iterate over rotations."""
        for i in range(18, 27):
            yield cls(i)


layer_axis = {
    Layer.NONE: Axis.NONE,
    Layer.UP: Axis.Y,
    Layer.FRONT: Axis.Z,
    Layer.RIGHT: Axis.X,
    Layer.DOWN: Axis.Y,
    Layer.BACK: Axis.Z,
    Layer.LEFT: Axis.X
}

layer_perm = {
    Layer.NONE: [Cubie.NONE],
    Layer.UP: [Cubie.UBL, Cubie.UBR, Cubie.UFR, Cubie.UFL],
    Layer.FRONT: [Cubie.UFR, Cubie.DFR, Cubie.DFL, Cubie.UFL],
    Layer.RIGHT: [Cubie.UFR, Cubie.UBR, Cubie.DBR, Cubie.DFR],
    Layer.DOWN: [Cubie.DBR, Cubie.DBL, Cubie.DFL, Cubie.DFR],
    Layer.BACK: [Cubie.UBL, Cubie.DBL, Cubie.DBR, Cubie.UBR],
    Layer.LEFT: [Cubie.UBL, Cubie.UFL, Cubie.DFL, Cubie.DBL]
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
    Face.FRONT: (slice(None), 1, slice(None)),
    Face.RIGHT: (slice(None), slice(None, None, -1), 1),
    Face.DOWN: (1, slice(None, None, -1), slice(None)),
    Face.BACK: (slice(None), 0, slice(None, None, -1)),
    Face.LEFT: (slice(None), slice(None), 0)
}

cubie_orbit = {
    Cubie.NONE: Orbit.NONE,
    # corners
    Cubie.UBL: Orbit.TETRAD_111,
    Cubie.UFR: Orbit.TETRAD_111,
    Cubie.DFL: Orbit.TETRAD_111,
    Cubie.DBR: Orbit.TETRAD_111,
    Cubie.UBR: Orbit.TETRAD_M11,
    Cubie.UFL: Orbit.TETRAD_M11,
    Cubie.DFR: Orbit.TETRAD_M11,
    Cubie.DBL: Orbit.TETRAD_M11,
}

cubie_index = {
    Cubie.NONE: (None, None, None),
    # corners
    Cubie.UBL: (0, 0, 0),
    Cubie.UFR: (0, 1, 1),
    Cubie.DFL: (1, 1, 0),
    Cubie.DBR: (1, 0, 1),
    Cubie.UBR: (0, 0, 1),
    Cubie.UFL: (0, 1, 0),
    Cubie.DFR: (1, 1, 1),
    Cubie.DBL: (1, 0, 0)
}

char_layer = {layer.char: layer for layer in Layer}
char_color = {color.char: color for color in Color}
char_face = {face.char: face for face in Face}
faces_cubie = {tuple(min([cb.faces[i:] + cb.faces[:i] for i in range(len(cb.faces))])if cb.faces else []): cb for cb in Cubie}
layer_order = [Layer.UP, Layer.FRONT, Layer.RIGHT, Layer.DOWN, Layer.BACK, Layer.LEFT]
str_move = {move.string: move for move in Move}
