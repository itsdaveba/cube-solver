"""Cube module."""
from __future__ import annotations

import warnings
import numpy as np
from copy import deepcopy
from typing import Union, Tuple, Dict

from .defs import NONE, SIZE, NUM_DIMS, NUM_CORNERS
from .defs import CORNER_ORIENTATION_SIZE, CORNER_PERMUTATION_SIZE
from .enums import Axis, Orbit, Layer, Color, Face, Cubie, Move
from . import utils

ORIENTATION_AXES = [Axis.Y, Axis.Z, Axis.X]
REPR_ORDER = [Face.UP, Face.LEFT, Face.FRONT, Face.RIGHT, Face.BACK, Face.DOWN]

DEFAULT_COLOR_SCHEME = {
    Face.UP: Color.WHITE,
    Face.FRONT: Color.GREEN,
    Face.RIGHT: Color.RED,
    Face.DOWN: Color.YELLOW,
    Face.BACK: Color.BLUE,
    Face.LEFT: Color.ORANGE
}

INDEX_TO_CUBIE = np.array([
    Cubie.UBL, Cubie.UFR, Cubie.DFL, Cubie.DBR,  # TETRAD_111 orbit
    Cubie.UBR, Cubie.UFL, Cubie.DFR, Cubie.DBL,  # TETRAD_M11 orbit
    Cubie.NONE], dtype=int)

CUBIE_TO_INDEX = np.zeros(max(len(Cubie), max(Cubie) + 1), dtype=int)
for cubie in INDEX_TO_CUBIE:
    CUBIE_TO_INDEX[cubie] = np.where(INDEX_TO_CUBIE == cubie)[0][0]
CUBIE_TO_INDEX[Cubie.NONE] = NONE

CARTESIAN_AXES = [*Axis.axes()]
min_orientation_axes = min([ORIENTATION_AXES[i:] + ORIENTATION_AXES[:i] for i in range(len(ORIENTATION_AXES))])
min_cartesian_axes = min([CARTESIAN_AXES[i:] + CARTESIAN_AXES[:i] for i in range(len(CARTESIAN_AXES))])
SHIFT_MULT = min_orientation_axes == min_cartesian_axes

SWAP_COLORS = {}  # swap colors along axis
for axis in CARTESIAN_AXES:
    shift = CARTESIAN_AXES.index(axis) - 1
    axes = np.roll(np.flip(np.roll(CARTESIAN_AXES, -shift)), shift)
    SWAP_COLORS[axis] = [Axis(ax).name for ax in axes]
COLORS_TYPE = [(axis.name, int) for axis in CARTESIAN_AXES]

warnings.simplefilter("always")


class Cube:
    def __init__(self, scramble: Union[str, None] = None, repr: Union[str, None] = None, random_state: bool = False):
        """
        Create :class:`Cube` object.

        Parameters
        ----------
        scramble : str or None, optional
            Initial scramble. If ``None``, no scramble is applied. Default is ``None``.
        repr : str or None, optional
            Cube string representation. If not ``None``, the ``scramble`` parameter is ignored and
            creates a cube with the given string representation. Default is ``None``.
            See `Notes` for the string representation format.
        random_state : bool, optional
            If ``True``, the ``scramble`` and ``repr`` parameters are ignored and
            creates a cube with a uniform random state. Default is ``False``.

        Notes
        -----
        The ``repr`` parameter must contain characters from `{'W', 'G', 'R', 'Y', 'B', 'O'}`,
        representing the colors :attr:`.Color.WHITE`, :attr:`.Color.GREEN`, :attr:`.Color.RED`,
        :attr:`.Color.YELLOW`, :attr:`.Color.BLUE`, and :attr:`.Color.ORANGE`, respectively.
        The order of the string representation is::

                    =========
                    | 01 02 |
                    ---------
                    | 03 04 |
            =================================
            | 05 06 | 09 10 | 13 14 | 17 18 |
            ---------------------------------
            | 07 08 | 11 12 | 15 16 | 19 20 |
            =================================
                    | 21 22 |
                    ---------
                    | 24 23 |
                    =========

        If the :attr:`orientation`, :attr:`permutation`, or :attr:`permutation_parity` values
        cannot be determined correctly from the string representation,
        the :attr:`orientation` and :attr:`permutation` arrays will contain ``-1`` at those positions,
        and :attr:`permutation_parity` will be set to ``None``.

        The default color scheme used for the cube is as follows (note: this may differ when using the ``repr`` parameter):

        * :attr:`.Face.UP`: :attr:`.Color.WHITE`
        * :attr:`.Face.FRONT`: :attr:`.Color.GREEN`
        * :attr:`.Face.RIGHT`: :attr:`.Color.RED`
        * :attr:`.Face.DOWN`: :attr:`.Color.YELLOW`
        * :attr:`.Face.BACK`: :attr:`.Color.BLUE`
        * :attr:`.Face.LEFT`: :attr:`.Color.ORANGE`

        Examples
        --------
        >>> from cube_solver import Cube

        Initial scramble.

        >>> cube = Cube("U F2 R'")
        >>> cube  # string representation of the cube state
        WBYOGROBGWRYBROGYOWBWGYR

        Initial string representation.

        >>> cube = Cube(repr="WBYOGROBGWRYBROGYOWBWGYR")
        >>> print(cube)  # print a visual layout of the cube state
              -------
              | W B |
              | Y O |
        -------------------------
        | G R | G W | B R | Y O |
        | O B | R Y | O G | W B |
        -------------------------
              | W G |
              | Y R |
              -------

        Initial random state.

        >>> cube = Cube(random_state=True)
        >>> cube.coords  # coordinates of the cube state (result might differ) # doctest: +SKIP
        (632, 3766)
        """
        if scramble is not None and not isinstance(scramble, str):
            raise TypeError(f"scramble must be str or None, not {type(scramble).__name__}")
        if repr is not None and not isinstance(repr, str):
            raise TypeError(f"repr must be str or None, not {type(repr).__name__}")
        if not isinstance(random_state, bool):
            raise TypeError(f"random_state must be bool, not {type(random_state).__name__}")

        self._color_scheme: Dict[Face, Color]
        """
        Color shceme of the cube.
        Used to generate and parse the string representation.
        """
        self._colors: np.ndarray
        """
        Color representation array of the cube.
        Used to generate and parse the string representation.
        """
        self.orientation: np.ndarray
        """
        Orientation array.

        The ``orientation`` array contains the orientation values of the ``8`` corners.

        A corner is correctly oriented when the `top` or `bottom` facelet of the corner piece matches either the `top` or
        `bottom` color of the cube. The cube's color scheme is determined with respect to the :attr:`.Cubie.DBL` corner.
        Corner orientation values are:

        * ``0`` if the corner is `correctly` oriented.
        * ``1`` if the corner is `twisted clockwise` relative to the correct orientation.
        * ``2`` if the corner is `twisted counter-clockwise` relative to the correct orientation.
        """
        self.permutation: np.ndarray
        """
        Permutation array.

        The ``permutation`` array contains the permutation values of the ``8`` corners.

        The `solved state` permutation goes from ``0`` to ``7`` for the corners. The piece ordering in the solved state is:

        * Corners: [``UBL``, ``UFR``, ``DFL``, ``DBR``, ``UBR``, ``UFL``, ``DFR``, ``DBL``]
        """
        self.permutation_parity: Union[bool, None]
        """
        Permutation parity.

        The ``permutation_parity`` indicates the parity of the `corner` permutation,
        ``True`` for ``odd`` parity, ``False`` for ``even`` parity,
        and ``None`` if the parity cannot be determined.

        The `solved state` permutation parity starts with ``even`` corner parity.
        """
        self.reset()
        if random_state:
            self.set_random_state()
        elif repr is not None:
            self._parse_repr(repr)
        elif scramble is not None:
            self.apply_maneuver(scramble)

    @property
    def coords(self) -> Tuple[int, int]:
        """
        Cube coordinates.

        `Corner orientation` and `corner permutation`.

        See Also
        --------
        get_coords
        set_coords

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube("U F2 R'")

        Get cube coordinates.

        >>> cube.coords
        (183, 3675)

        Set cube coordinates.

        >>> cube.coords = (0, 0)  # solved state
        >>> cube
        WWWWOOOOGGGGRRRRBBBBYYYY
        """
        return self.get_coords()

    @coords.setter
    def coords(self, coords: Tuple[int, int]):
        self.set_coords(coords)

    @property
    def is_solved(self) -> bool:
        """
        Whether the cube is solved.

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube()
        >>> cube.is_solved
        True
        >>> cube.apply_maneuver("U F2 R'")
        >>> cube.is_solved
        False
        """
        return self.coords == (0, 0)

    def __ne__(self, other: object) -> bool:
        """Negation of equality comparison."""
        return not self.__eq__(other)

    def __eq__(self, other: object) -> bool:
        """Equality comparison."""
        if not isinstance(other, Cube):
            return False
        return repr(self) == repr(other)

    def __repr__(self) -> str:
        """String representation of the :class:`Cube` object."""
        solved_colors = np.full_like(self._colors, Color.NONE)
        for face in Face.faces():
            solved_colors[face._index][face.axis.name] = self._color_scheme[face]

        for cubie in Cubie.cubies():
            cubie_orientation = self.orientation[CUBIE_TO_INDEX[cubie]]
            cubie_permutation = Cubie(INDEX_TO_CUBIE[self.permutation[CUBIE_TO_INDEX[cubie]]])
            cubie_colors = np.array(solved_colors[cubie_permutation._index], dtype=COLORS_TYPE)
            if cubie_permutation != Cubie.NONE:
                if cubie.orbit != cubie_permutation.orbit:
                    cubie_colors = np.array(cubie_colors[SWAP_COLORS[Axis.Y]], dtype=COLORS_TYPE)
                if cubie_orientation:
                    shift = (cubie_orientation if cubie.orbit == Orbit.TETRAD_M11 else -cubie_orientation) * SHIFT_MULT
                    cubie_colors = np.array(tuple(np.roll(cubie_colors.tolist(), shift)), dtype=COLORS_TYPE)
            else:
                cubie_colors = Color.NONE
            self._colors[cubie._index] = cubie_colors

        return "".join([Color(c).char for c in np.ravel([self._colors[face._index][face.axis.name] for face in REPR_ORDER])])

    def __str__(self) -> str:
        """Print representation of the `Cube` object."""
        repr = self.__repr__()

        # up face
        str = "  " * SIZE + "  " + "--" * SIZE + "---\n"
        for i in range(SIZE):
            j = REPR_ORDER.index(Face.UP) * SIZE * SIZE + i * SIZE
            str += "  " * SIZE + "  | " + " ".join(repr[j:j+SIZE]) + " |\n"

        # lateral faces
        str += "--------" * SIZE + "---------\n"
        for i in range(SIZE):
            js = [REPR_ORDER.index(face) * SIZE * SIZE + i * SIZE for face in [Face.LEFT, Face.FRONT, Face.RIGHT, Face.BACK]]
            str += "| " + " | ".join(" ".join(repr[j:j+SIZE]) for j in js) + " |\n"
        str += "--------" * SIZE + "---------\n"

        # down face
        for i in range(SIZE):
            j = REPR_ORDER.index(Face.DOWN) * SIZE * SIZE + i * SIZE
            str += "  " * SIZE + "  | " + " ".join(repr[j:j+SIZE]) + " |\n"
        str += "  " * SIZE + "  " + "--" * SIZE + "---"

        return str

    def _parse_repr(self, repr: str):
        """
        Parse a string representation into a cube state.

        If the :attr:`orientation`, :attr:`permutation`, or :attr:`permutation_parity` values
        cannot be determined correctly from the string representation,
        the :attr:`orientation` and :attr:`permutation` arrays will contain ``-1`` at those positions,
        and :attr:`permutation_parity` will be set to ``None``.

        Parameters
        ----------
        repr : str
            String representation of the cube.
        """
        if len(repr) != len(REPR_ORDER)*SIZE*SIZE:
            raise ValueError(f"repr length must be {len(REPR_ORDER)*SIZE*SIZE} (got {len(repr)})")

        face_repr = [repr[i:i+SIZE*SIZE] for i in range(0, len(REPR_ORDER)*SIZE*SIZE, SIZE*SIZE)]
        for face, _repr in zip(REPR_ORDER, face_repr):
            self._colors[face._index][face.axis.name] = np.reshape([*map(Color.from_char, _repr)], (SIZE, SIZE))

        # color scheme
        inv_color_scheme = {color: Face.NONE for color in Color.colors()}
        for face in Face.faces():
            color = Color(self._colors[Cubie.DBL._index][face.axis.name])
            if face in (Face.UP, Face.FRONT, Face.RIGHT):
                try:
                    default_face = [*DEFAULT_COLOR_SCHEME.keys()][[*DEFAULT_COLOR_SCHEME.values()].index(color)]
                    color = DEFAULT_COLOR_SCHEME[default_face.opposite]
                except ValueError:
                    color = Color.NONE
            self._color_scheme[face] = color
        inv_color_scheme.update({color: face for face, color in self._color_scheme.items()})

        for cubie in Cubie.cubies():
            cubie_colors = [self._colors[cubie._index][axis.name] for axis in ORIENTATION_AXES]
            cubie_faces = np.array([inv_color_scheme[color] for color in cubie_colors if color != Color.NONE], dtype=Face)
            try:
                if cubie.orbit == Orbit.TETRAD_111:
                    x, z = [ORIENTATION_AXES.index(axis) for axis in (Axis.X, Axis.Z)]
                    cubie_faces[x], cubie_faces[z] = cubie_faces[z], cubie_faces[x]  # swap along `Y` axis
                self.orientation[CUBIE_TO_INDEX[cubie]] = [face.axis for face in cubie_faces].index(Axis.Y)
                self.permutation[CUBIE_TO_INDEX[cubie]] = CUBIE_TO_INDEX[Cubie.from_faces(cubie_faces.tolist())]
            except Exception:
                self.orientation[CUBIE_TO_INDEX[cubie]] = NONE
                self.permutation[CUBIE_TO_INDEX[cubie]] = CUBIE_TO_INDEX[Cubie.NONE]

        if np.any(self.permutation == CUBIE_TO_INDEX[Cubie.NONE]):
            warnings.warn("invalid string representation, setting undefined orientation and permutation values with -1")
            self.permutation_parity = None
        else:
            if np.sum(self.orientation) % 3 != 0:
                warnings.warn("invalid corner orientation")
            if len(set(self.permutation)) == NUM_CORNERS:
                self.permutation_parity = utils.get_permutation_parity(self.permutation)
            else:
                warnings.warn("invalid corner permutation")
                self.permutation_parity = None

    def reset(self):
        """
        Reset the cube to the solved state using the default color scheme.

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube(random_state=True)
        >>> cube.reset()
        >>> cube
        WWWWOOOOGGGGRRRRBBBBYYYY
        """
        self._color_scheme = DEFAULT_COLOR_SCHEME.copy()
        self._colors = np.full((SIZE,) * NUM_DIMS, Color.NONE, dtype=COLORS_TYPE)
        self.orientation = np.zeros(NUM_CORNERS, dtype=int)
        self.permutation = np.arange(NUM_CORNERS, dtype=int)
        self.permutation_parity = False

    def set_random_state(self):
        """
        Set a uniform random state.

        Sets random :attr:`coords` (`corner orientation` and `corner permutation`).

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube()
        >>> cube.set_random_state()
        >>> cube.coords  # result might differ # doctest: +SKIP
        (632, 3766)
        """
        self.set_coord("co", np.random.randint(CORNER_ORIENTATION_SIZE))
        self.set_coord("cp", np.random.randint(CORNER_PERMUTATION_SIZE))

    def apply_move(self, move: Move):
        """
        Apply a move to the cube.

        Parameters
        ----------
        move : Move
            Move to apply.

        Examples
        --------
        >>> from cube_solver import Cube, Move
        >>> cube = Cube()
        >>> cube.apply_move(Move.U1)   # U face move
        >>> cube.apply_move(Move.X1)   # x rotation
        >>> cube
        RRGGGOGOYYYYRBRBWWWWBBOO
        """
        if not isinstance(move, Move):
            raise TypeError(f"move must be Move, not {type(move).__name__}")

        if move.is_face:
            layer = move.layers[0]
            shift = move.shifts[0]
            cubies = np.roll(layer.perm, shift, axis=-1)
            orientation = self.orientation[CUBIE_TO_INDEX[cubies]]
            if shift % 2 == 1:
                if move.axis == Axis.Z:
                    orientation = np.where(orientation != NONE, (orientation + [1, 2, 1, 2]) % 3, NONE)
                elif move.axis == Axis.X:
                    orientation = np.where(orientation != NONE, (orientation + [2, 1, 2, 1]) % 3, NONE)
                if self.permutation_parity is not None:
                    self.permutation_parity = not self.permutation_parity
            self.orientation[CUBIE_TO_INDEX[layer.perm]] = orientation
            self.permutation[CUBIE_TO_INDEX[layer.perm]] = self.permutation[CUBIE_TO_INDEX[cubies]]
            if layer in (Layer.DOWN, Layer.BACK, Layer.LEFT):
                self._parse_repr(repr(self))

        elif move.is_rotation:
            for layer, shift in zip(move.layers, move.shifts):
                self.apply_move(Move[layer.char + str(shift % 4)])

    def apply_maneuver(self, maneuver: str):
        """
        Apply a sequence of moves to the cube.

        Accepts the following move types:

        * Face moves (e.g. `U`, `F2`, `R'`).
        * Rotations (e.g. `x`, `y2`, `z'`).

        Parameters
        ----------
        maneuver : str
            The sequence of moves to apply.

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube()
        >>> cube.apply_maneuver("U x")
        >>> cube
        RRGGGOGOYYYYRBRBWWWWBBOO
        """
        if not isinstance(maneuver, str):
            raise TypeError(f"maneuver must be str, not {type(maneuver).__name__}")

        # get moves from attr `moves` if maneuver is an instance of the `Maneuver` str subclass
        moves = getattr(maneuver, "moves", [Move.from_string(move_str) for move_str in maneuver.split()])
        for move in moves:
            self.apply_move(move)

    def get_coord(self, coord_name: str) -> int:
        """
        Get cube coordinate value.

        Parameters
        ----------
        coord_name : {'co', 'cp'}
            Get the specified cube coordinate.

            * 'co' means `corner orientation`.
            * 'cp' means `corner permutation`.

        Returns
        -------
        coord : int
            Cube coordinate value.

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube("U F2 R'")

        Get corner coordinates.

        >>> cube.get_coord('co')   # corner orientation
        183
        >>> cube.get_coord('cp')   # corner permutation
        3675
        """
        if not isinstance(coord_name, str):
            raise TypeError(f"coord_name must be str, not {type(coord_name).__name__}")

        if coord_name == "co":
            return utils.get_orientation_coord(self.orientation[:-1], 3, is_modulo=True)
        if coord_name == "cp":
            return utils.get_permutation_coord(self.permutation[:-1])

        raise ValueError(f"coord_name must be one of 'co', 'cp' (got '{coord_name}')")

    def set_coord(self, coord_name: str, coord: int):
        """
        Set cube coordinate value.

        Parameters
        ----------
        coord_name : {'co', 'cp'}
            Set the specified cube coordinate.

            * 'co' means `corner orientation`.
            * 'cp' means `corner permutation`.

        coord : int
            Cube coordinate value.

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube()

        Set corner coordinates.

        >>> cube.set_coord('co', 183)   # corner orientation
        >>> cube.orientation
        array([0, 2, 0, 2, 1, 0, 1, 0])
        >>> cube.set_coord('cp', 3675)  # corner permutation
        >>> cube.permutation
        array([5, 0, 4, 1, 3, 6, 2, 7])
        """
        if not isinstance(coord_name, str):
            raise TypeError(f"coord_name must be str, not {type(coord_name).__name__}")
        if not isinstance(coord, int):
            raise TypeError(f"coord must be int, not {type(coord).__name__}")

        if coord_name == "co":
            orientation = utils.get_orientation_array(coord, 3, len(self.orientation) - 1, force_modulo=True)
            self.orientation = np.concatenate([orientation, [0]])
        elif coord_name == "cp":
            permutation, self.permutation_parity = utils.get_permutation_array(coord, len(self.permutation) - 1)
            self.permutation = np.concatenate([permutation, [NUM_CORNERS - 1]])
        else:
            raise ValueError(f"coord_name must be one of 'co', 'cp' (got '{coord_name}')")

    def get_coords(self) -> Tuple[int, int]:
        """
        Get cube coordinates.

        Get the `corner orientation` and `corner permutation` coordinates.

        Returns
        -------
        coords : tuple of (int, int)
            Cube coordinates in the following order:
            `corner orientation`, corner permutation`.

        See Also
        --------
        get_coord

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube("U F2 R'")

        Get cube coordinates.

        >>> cube.get_coords()
        (183, 3675)
        """
        return (self.get_coord("co"), self.get_coord("cp"))

    def set_coords(self, coords: Tuple[int, int]):
        """
        Set cube coordinates.

        Set the `corner orientation` and `corner permutation`.

        Parameters
        ----------
        coords : tuple of (int, int)
            Cube coordinates in the following order:
            `corner orientation`, `corner permutation`.

        See Also
        --------
        set_coord

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube()

        Set cube coordinates.

        >>> coords = (183, 3675)
        >>> cube.set_coords(coords)
        >>> cube
        WBYOGROBGWRYBROGYOWBWGYR
        """
        if not isinstance(coords, tuple):
            raise TypeError(f"coords must be tuple, not {type(coords).__name__}")
        if len(coords) != 2:
            raise ValueError(f"coords tuple length must be 2 (got {len(coords)})")

        self.set_coord("co", coords[0])
        self.set_coord("cp", coords[1])

    def copy(self) -> Cube:
        """
        Return a copy of the cube.

        Returns
        -------
        cube : Cube
            Copy of cube object.

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube("U F2 R'")
        >>> cube_copy = cube.copy()
        >>> cube_copy
        WBYOGROBGWRYBROGYOWBWGYR
        >>> cube_copy == cube
        True
        """
        return deepcopy(self)


def apply_move(cube: Cube, move: Move) -> Cube:
    """
    Return a copy of the the cube with the move applyed.

    Parameters
    ----------
    cube : Cube
        Cube object.
    move : Move
        Move to apply.

    Returns
    -------
    cube : Cube
        Copy of the cube with the move applied.

    Examples
    --------
    >>> from cube_solver import Cube, Move, apply_move
    >>> cube = Cube()
    >>> apply_move(cube, Move.U1)   # U face move
    WWWWGGOORRGGBBRROOBBYYYY
    >>> apply_move(cube, Move.X1)   # x rotation
    GGGGOOOOYYYYRRRRWWWWBBBB
    """
    if not isinstance(cube, Cube):
        raise TypeError(f"cube must be Cube, not {type(cube).__name__}")

    cube = cube.copy()
    cube.apply_move(move)
    return cube


def apply_maneuver(cube: Cube, maneuver: str) -> Cube:
    """
    Return a copy of the cube with the sequence of moves applied.

    Accepts the following move types:

    * Face moves (e.g. `U`, `F2`, `R'`).
    * Rotations (e.g. `x`, `y2`, `z'`).

    Parameters
    ----------
    cube : Cube
        Cube object.
    maneuver : str
        The sequence of moves to apply.

    Returns
    -------
    cube : Cube
        Copy of the cube with the sequence of moves applied.

    Examples
    --------
    >>> from cube_solver import Cube, apply_maneuver
    >>> cube = Cube()
    >>> apply_maneuver(cube, "U x")
    RRGGGOGOYYYYRBRBWWWWBBOO
    """
    if not isinstance(cube, Cube):
        raise TypeError(f"cube must be Cube, not {type(cube).__name__}")

    cube = cube.copy()
    cube.apply_maneuver(maneuver)
    return cube
