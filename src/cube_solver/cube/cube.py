"""Cube module."""
from __future__ import annotations

import math
import warnings
import numpy as np
from copy import deepcopy
from itertools import chain
from typing import Union, Tuple, Dict

from ..defs import CoordType, CoordsType
from .defs import NONE, SIZE, NUM_DIMS, NUM_CORNERS, NUM_EDGES, NUM_ORBIT_ELEMS
from .defs import CORNER_ORIENTATION_SIZE, EDGE_ORIENTATION_SIZE, CORNER_PERMUTATION_SIZE, EDGE_PERMUTATION_SIZE
from .enums import Axis, Orbit, Layer, Color, Face, Cubie, Move
from . import utils

ORIENTATION_AXES = [Axis.Y, Axis.Z, Axis.X]
REPR_ORDER = [Face.UP, Face.LEFT, Face.FRONT, Face.RIGHT, Face.BACK, Face.DOWN]

CORNER_ORBITS = [Orbit.TETRAD_111, Orbit.TETRAD_M11]
EDGE_ORBITS = [Orbit.SLICE_MIDDLE, Orbit.SLICE_EQUATOR, Orbit.SLICE_STANDING]

DEFAULT_COLOR_SCHEME = {
    Face.UP: Color.WHITE,
    Face.FRONT: Color.GREEN,
    Face.RIGHT: Color.RED,
    Face.DOWN: Color.YELLOW,
    Face.BACK: Color.BLUE,
    Face.LEFT: Color.ORANGE
}

INDEX_TO_CUBIE = np.array([
    Cubie.UBL, Cubie.UFR, Cubie.DBR, Cubie.DFL,  # TETRAD_111 orbit
    Cubie.UBR, Cubie.UFL, Cubie.DBL, Cubie.DFR,  # TETRAD_M11 orbit
    Cubie.UB, Cubie.UF, Cubie.DB, Cubie.DF,      # SLICE_MIDDLE orbit
    Cubie.UL, Cubie.UR, Cubie.DL, Cubie.DR,      # SLICE_STANDING orbit
    Cubie.BL, Cubie.BR, Cubie.FL, Cubie.FR,      # SLICE_EQUATOR orbit
    Cubie.NONE], dtype=int)

CUBIE_TO_INDEX = np.zeros(max(len(Cubie), max(Cubie) + 1), dtype=int)
for cubie in INDEX_TO_CUBIE:
    CUBIE_TO_INDEX[cubie] = np.where(INDEX_TO_CUBIE == cubie)[0][0]
CUBIE_TO_INDEX[Cubie.NONE] = NONE

CARTESIAN_AXES = [*Axis.cartesian_axes()]
min_orientation_axes = min([ORIENTATION_AXES[i:] + ORIENTATION_AXES[:i] for i in range(len(ORIENTATION_AXES))])
min_cartesian_axes = min([CARTESIAN_AXES[i:] + CARTESIAN_AXES[:i] for i in range(len(CARTESIAN_AXES))])
SHIFT_MULT = min_orientation_axes == min_cartesian_axes

SWAP_COLORS = {}  # swap colors along axis
for axis in CARTESIAN_AXES:
    shift = CARTESIAN_AXES.index(axis) - 1
    axes = np.roll(np.flip(np.roll(CARTESIAN_AXES, -shift)), shift)
    SWAP_COLORS[axis] = [Axis(ax).name for ax in axes]
COLORS_TYPE = [(axis.name, int) for axis in CARTESIAN_AXES]

ORBIT_OFFSET = {
    Orbit.TETRAD_111: CUBIE_TO_INDEX[Cubie.UBL],
    Orbit.TETRAD_M11: CUBIE_TO_INDEX[Cubie.UBR],
    Orbit.SLICE_MIDDLE: CUBIE_TO_INDEX[Cubie.UB],
    Orbit.SLICE_STANDING: CUBIE_TO_INDEX[Cubie.UL],
    Orbit.SLICE_EQUATOR: CUBIE_TO_INDEX[Cubie.BL]
}

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

                       ------------
                       | 01 02 03 |
                       | 04 05 06 |
                       | 07 08 09 |
            ---------------------------------------------
            | 10 11 12 | 19 20 21 | 28 29 30 | 37 38 39 |
            | 13 14 15 | 22 23 24 | 31 32 33 | 40 41 42 |
            | 16 17 18 | 25 26 27 | 34 35 36 | 43 44 45 |
            ---------------------------------------------
                       | 46 47 48 |
                       | 49 50 51 |
                       | 52 53 54 |
                       ------------

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
        WWBWWBYYOGGROOROOBGGWGGWRRYBRRBRROOGYOOYBBWBBWWGYYGYYR

        Initial string representation.

        >>> cube = Cube(repr="WWBWWBYYOGGROOROOBGGWGGWRRYBRRBRROOGYOOYBBWBBWWGYYGYYR")
        >>> print(cube)  # print a visual layout of the cube state
                ---------
                | W W B |
                | W W B |
                | Y Y O |
        ---------------------------------
        | G G R | G G W | B R R | Y O O |
        | O O R | G G W | B R R | Y B B |
        | O O B | R R Y | O O G | W B B |
        ---------------------------------
                | W W G |
                | Y Y G |
                | Y Y R |
                ---------

        Initial random state.

        >>> cube = Cube(random_state=True)
        >>> cube.coords  # coordinates of the cube state (result might differ) # doctest: +SKIP
        (167, 48, 22530, 203841327)
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

        The ``orientation`` array contains the orientation values of the ``8`` corners and ``12`` edges of the cube.
        The first ``8`` elements represent the `corner` orientation values, and the remaining ``12`` elements represent the
        `edge` orientation values.

        A corner is correctly oriented when the `top` or `bottom` facelet of the corner piece matches eihter the `top` or
        `bottom` color of the cube. Corner orientation values are:

        * ``0`` if the corner is `correctly` oriented.
        * ``1`` if the corner is `twisted clockwise` relative to the correct orientation.
        * ``2`` if the corner is `twisted counter-clockwise` relative to the correct orientation.

        An edge is correctly oriented if, when placed in its correct position using only
        :attr:`.Face.UP`, :attr:`.Face.DOWN`, :attr:`.Face.RIGHT` and :attr:`.Face.LEFT` face turns,
        it does not appear `flipped`. Edge orientation values are:

        * ``0`` if the edge is `correctly` oriented.
        * ``1`` if the edge is incorrectly oriented (i.e. `flipped`).
        """
        self.permutation: np.ndarray
        """
        Permutation array.

        The ``permutation`` array contains the permutation values of the ``8`` corners and ``12`` edges of the cube.
        The first ``8`` elements represent the `corner` permutation values, and the remaining ``12`` elements represent the
        `edge` permutation values.

        The `solved state` permutation goes from ``0`` to ``7`` for the corners,
        and from ``8`` to ``19`` for the edges. The piece ordering in the solved state is:

        * Corners: [``UBL``, ``UFR``, ``DBR``, ``DFL``, ``UBR``, ``UFL``, ``DBL``, ``DFR``]
        * Edges: [``UB``, ``UF``, ``DB``, ``DF``, ``UL``, ``UR``, ``DL``, ``DR``, ``BL``, ``BR``, ``FL``, ``FR``]
        """
        self.permutation_parity: Union[bool, None]
        """
        Permutation parity.

        The ``permutation_parity`` indicates the parity of both `corner` and `edge` permutations
        (i.e. both parities are always the same), ``True`` for ``odd`` parity, ``False`` for ``even`` parity,
        and ``None`` if the parity cannot be determined.

        The `solved state` permutation parity starts with ``even`` corner and endge parity.
        """
        self.reset()
        if random_state:
            self.set_random_state()
        elif repr is not None:
            self._parse_repr(repr)
        elif scramble is not None:
            self.apply_maneuver(scramble)

    @property
    def coords(self) -> Tuple[int, ...]:
        """
        Cube coordinates.

        `Corner orientation`, `edge orientation`,
        `corner permutation`, and `edge permutation` coordinates.

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
        (657, 0, 25253, 85684063)

        Set cube coordinates.

        >>> cube.coords = (0, 0, 0, 0)  # solved state
        >>> cube
        WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY
        """
        coords = self.get_coords()
        return tuple(coord if isinstance(coord, int) else coord[0] for coord in coords)

    @coords.setter
    def coords(self, coords: Tuple[int, ...]):
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
        return self.coords == (0, 0, 0, 0)

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
        solved_colors = np.full_like(self._colors, (Color.NONE,) * NUM_DIMS)
        for face in Face.faces():
            solved_colors[face._index][face.axis.name] = self._color_scheme[face]

        # centers
        for center in Cubie.centers():
            self._colors[center._index] = solved_colors[center._index]

        # corners and edges
        for cubie in chain(Cubie.corners(), Cubie.edges()):
            cubie_orientation = self.orientation[CUBIE_TO_INDEX[cubie]]
            cubie_permutation = Cubie(INDEX_TO_CUBIE[self.permutation[CUBIE_TO_INDEX[cubie]]])
            cubie_colors = np.array(solved_colors[cubie_permutation._index], dtype=COLORS_TYPE)
            if cubie_permutation != Cubie.NONE:
                if cubie.is_corner:
                    if cubie.orbit != cubie_permutation.orbit:
                        cubie_colors = np.array(cubie_colors[SWAP_COLORS[Axis.Y]], dtype=COLORS_TYPE)
                    if cubie_orientation:
                        shift = (cubie_orientation if cubie.orbit == Orbit.TETRAD_M11 else -cubie_orientation) * SHIFT_MULT
                        cubie_colors = np.array(tuple(np.roll(cubie_colors.tolist(), shift)), dtype=COLORS_TYPE)
                else:
                    orbits = {cubie.orbit, cubie_permutation.orbit}
                    if orbits == {Orbit.SLICE_MIDDLE, Orbit.SLICE_EQUATOR}:
                        shift = (1 if cubie.orbit == Orbit.SLICE_EQUATOR else -1) * SHIFT_MULT
                        cubie_colors = np.array(tuple(np.roll(cubie_colors.tolist(), shift)), dtype=COLORS_TYPE)
                    elif orbits == {Orbit.SLICE_MIDDLE, Orbit.SLICE_STANDING}:
                        cubie_colors = np.array(cubie_colors[SWAP_COLORS[Axis.Y]], dtype=COLORS_TYPE)
                    elif orbits == {Orbit.SLICE_EQUATOR, Orbit.SLICE_STANDING}:
                        cubie_colors = np.array(cubie_colors[SWAP_COLORS[Axis.X]], dtype=COLORS_TYPE)
                    if cubie_orientation:
                        axis = Layer[cubie.orbit.name.split("_")[1]].axis
                        cubie_colors = np.array(cubie_colors[SWAP_COLORS[axis]], dtype=COLORS_TYPE)
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

        # centers
        inv_color_scheme = {color: Face.NONE for color in Color.colors()}
        for center in Cubie.centers():
            self._color_scheme[Face.from_char(center.name)] = Color(self._colors[center._index][center.axis.name])
        inv_color_scheme.update({color: face for face, color in self._color_scheme.items()})

        # corners and edges
        for cubie in chain(Cubie.corners(), Cubie.edges()):
            cubie_colors = [self._colors[cubie._index][axis.name] for axis in ORIENTATION_AXES]
            cubie_faces = np.array([inv_color_scheme[color] for color in cubie_colors if color != Color.NONE], dtype=Face)
            try:
                if cubie.is_corner:
                    if cubie.orbit == Orbit.TETRAD_111:
                        x, z = [ORIENTATION_AXES.index(axis) for axis in (Axis.X, Axis.Z)]
                        cubie_faces[x], cubie_faces[z] = cubie_faces[z], cubie_faces[x]  # swap along `Y` axis
                    self.orientation[CUBIE_TO_INDEX[cubie]] = [face.axis for face in cubie_faces].index(Axis.Y)
                else:
                    axis_importance = [ORIENTATION_AXES.index(face.axis) for face in cubie_faces]
                    self.orientation[CUBIE_TO_INDEX[cubie]] = np.diff(axis_importance)[0] < 0
                self.permutation[CUBIE_TO_INDEX[cubie]] = CUBIE_TO_INDEX[Cubie.from_faces(cubie_faces.tolist())]
            except Exception:
                self.orientation[CUBIE_TO_INDEX[cubie]] = NONE
                self.permutation[CUBIE_TO_INDEX[cubie]] = CUBIE_TO_INDEX[Cubie.NONE]

        if np.any(self.permutation == CUBIE_TO_INDEX[Cubie.NONE]):
            warnings.warn("invalid string representation, setting undefined orientation and permutation values with -1")
            self.permutation_parity = None
        else:
            if np.sum(self.orientation[:NUM_CORNERS]) % 3 != 0:
                warnings.warn("invalid corner orientation")
            if np.sum(self.orientation[NUM_CORNERS:]) % 2 != 0:
                warnings.warn("invalid edge orientation")
            is_valid_corner_perm = len(set(self.permutation[:NUM_CORNERS])) == NUM_CORNERS
            is_valid_edge_perm = len(set(self.permutation[NUM_CORNERS:])) == NUM_EDGES
            if is_valid_corner_perm and is_valid_edge_perm:
                corner_parity = utils.get_permutation_parity(self.permutation[:NUM_CORNERS])
                edge_parity = utils.get_permutation_parity(self.permutation[NUM_CORNERS:])
                if corner_parity == edge_parity:
                    self.permutation_parity = corner_parity
                else:
                    warnings.warn("invalid cube parity")
                    self.permutation_parity = None
            else:
                if not is_valid_corner_perm:
                    warnings.warn("invalid corner permutation")
                if not is_valid_edge_perm:
                    warnings.warn("invalid edge permutation")
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
        WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY
        """
        self._color_scheme = DEFAULT_COLOR_SCHEME.copy()
        self._colors = np.full((SIZE,) * NUM_DIMS, (Color.NONE,) * NUM_DIMS, dtype=COLORS_TYPE)
        self.orientation = np.zeros(NUM_CORNERS + NUM_EDGES, dtype=int)
        self.permutation = np.arange(NUM_CORNERS + NUM_EDGES, dtype=int)
        self.permutation_parity = False

    def set_random_state(self):
        """
        Set a uniform random state.

        Sets random :attr:`coords` (`corner orientation`, `edge orientation`,
        `corner permutation`, and `edge permutation` coordinates).

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube()
        >>> cube.set_random_state()
        >>> cube.coords  # result might differ # doctest: +SKIP
        (167, 48, 22530, 203841327)
        """
        self.set_coord("co", np.random.randint(CORNER_ORIENTATION_SIZE))
        self.set_coord("eo", np.random.randint(EDGE_ORIENTATION_SIZE))
        self.set_coord("cp", np.random.randint(CORNER_PERMUTATION_SIZE))
        self.set_coord("ep", np.random.randint(EDGE_PERMUTATION_SIZE))

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
        >>> cube.apply_move(Move.M2)   # M2 slice move
        >>> cube.apply_move(Move.FW3)  # Fw' wide move
        >>> cube.apply_move(Move.X1)   # x rotation
        >>> cube
        RGGBBORGGWYWWYWGOOGOOGOOYWYYWYYWYRRBRRBRRBWYWBRBBGBOGO
        """
        if not isinstance(move, Move):
            raise TypeError(f"move must be Move, not {type(move).__name__}")

        if move.is_face:
            layer = move.layers[0]
            shift = move.shifts[0]
            cubies = np.roll(layer.perm, shift, axis=1)
            orientation = self.orientation[CUBIE_TO_INDEX[cubies]]
            if shift % 2 == 1:
                if move.axis == Axis.Z:
                    orientation = np.where(orientation != NONE, (orientation + [[1, 2, 1, 2], [1] * 4]) % ([3], [2]), NONE)
                elif move.axis == Axis.X:
                    orientation[0] = np.where(orientation[0] != NONE, (orientation[0] + [2, 1, 2, 1]) % 3, NONE)
            self.orientation[CUBIE_TO_INDEX[layer.perm]] = orientation
            self.permutation[CUBIE_TO_INDEX[layer.perm]] = self.permutation[CUBIE_TO_INDEX[cubies]]
            if self.permutation_parity is not None:
                self.permutation_parity = not self.permutation_parity

        elif move.is_slice:
            shift = move.shifts[0]
            base_move = Move[move.axis.name + "1"].layers[move.axis == Axis.Z].char
            self.apply_move(Move[base_move + "W" + str(-shift % 4)])  # wide move
            self.apply_move(Move[base_move + str(shift % 4)])  # face move

        elif move.is_wide:
            shift = move.shifts[0]
            mult = 1 if Move[move.axis.name + "1"].layers[0].char == move.name[0] else -1
            self.apply_move(Move[move.axis.name + str(mult * shift % 4)])  # rotation
            self.apply_move(Move[Face.from_char(move.name[0]).opposite.char + str(shift % 4)])  # face move

        elif move.is_rotation:
            layers_perm = [layer.perm for layer in move.layers]
            layers_shifted = [np.roll(layer, shift, axis=1) for layer, shift in zip(layers_perm, move.shifts)]
            rotation = {center: center for center in Cubie.centers()}
            rotation.update({Cubie(key): Cubie(val) for key, val in zip(np.ravel(layers_shifted), np.ravel(layers_perm))})
            rotation[Cubie.NONE] = Cubie.NONE

            # centers
            color_scheme = self._color_scheme.copy()
            for center in Cubie.centers():
                self._color_scheme[Face.from_char(rotation[center].name)] = color_scheme[Face.from_char(center.name)]

            # corners and edges
            cubies = [*Cubie.corners()] + [*Cubie.edges()]
            rotations = [rotation[cubie] for cubie in cubies]
            orientation = self.orientation[CUBIE_TO_INDEX[cubies]]
            cubie_permutation = [Cubie(INDEX_TO_CUBIE[perm]) for perm in self.permutation[CUBIE_TO_INDEX[cubies]]]
            if move.shifts[0] % 2 == 1:
                is_corner = np.array([cubie.is_corner for cubie in cubies])
                cubie_orbits = np.array([cubie.orbit for cubie in cubies])
                perm_orbits = np.array([perm.orbit for perm in cubie_permutation])
                corner_comp = edge_comp = [perm_orbits, cubie_orbits]
                if move.axis == Axis.X:
                    corner_comp = [cubie_orbits, Orbit.TETRAD_M11]
                    edge_comp = [Orbit.SLICE_MIDDLE, Orbit.SLICE_MIDDLE]
                elif move.axis == Axis.Y:
                    edge_comp = [Orbit.SLICE_EQUATOR, Orbit.SLICE_EQUATOR]
                else:
                    corner_comp = [cubie_orbits, Orbit.TETRAD_111]
                axis_comp = perm_orbits != np.where(is_corner, corner_comp[0], edge_comp[0])
                condition = cubie_orbits == np.where(is_corner, corner_comp[1], edge_comp[1])
                incr = np.where(condition, axis_comp, np.where(is_corner, -axis_comp.astype(int), ~axis_comp))
                orientation = np.where(orientation != NONE, (orientation + incr) % np.where(is_corner, 3, 2), NONE)
            self.orientation[CUBIE_TO_INDEX[rotations]] = orientation
            self.permutation[CUBIE_TO_INDEX[rotations]] = [CUBIE_TO_INDEX[rotation[perm]] for perm in cubie_permutation]

    def apply_maneuver(self, maneuver: str):
        """
        Apply a sequence of moves to the cube.

        Accepts the following move types:

        * Face moves (e.g. `U`, `F2`, `R'`).
        * Slice moves (e.g. `M`, `E2`, `S'`).
        * Wide moves (e.g. `Uw`, `Fw2`, `Rw'` or `u`, `f2`, `r'`).
        * Rotations (e.g. `x`, `y2`, `z'`).

        Parameters
        ----------
        maneuver : str
            The sequence of moves to apply.

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube()
        >>> cube.apply_maneuver("U M2 Fw' x")
        >>> cube
        RGGBBORGGWYWWYWGOOGOOGOOYWYYWYYWYRRBRRBRRBWYWBRBBGBOGO
        """
        if not isinstance(maneuver, str):
            raise TypeError(f"maneuver must be str, not {type(maneuver).__name__}")

        # get moves from attr `moves` if maneuver is an instance of the `Maneuver` str subclass
        moves = getattr(maneuver, "moves", [Move.from_string(move_str) for move_str in maneuver.split()])
        for move in moves:
            self.apply_move(move)

    def get_coord(self, coord_name: str) -> CoordType:
        """
        Get cube coordinate value.

        Parameters
        ----------
        coord_name : {'co', 'eo', 'cp', 'ep', 'pcp', 'pep'}
            Get the specified cube coordinate.

            * 'co' means `corner orientation`.
            * 'eo' means `edge orientation`.
            * 'cp' means `corner permutation`.
            * 'ep' means `edge permutation`.
            * 'pcp' means `partial corner permutation` (a value for each corner `orbit`).
            * 'pep' means `partial edge permutation` (a value for each edge `orbit`).

        Returns
        -------
        coord : int or tuple of int
            Cube coordinate value. For partial coordinate values, a value of ``-1``
            indicates an `orbit` with no permutation values (e.g., ``(-1, -1, 0)``),
            meaning the permutation (and orientation) values for that `orbit` are set to ``-1``.
            If only the value of the first `orbit` is available (i.e., ``(coord, -1, -1)``),
            returns an `int` representing that partial coordinate value.

            The `orbit` order for partial coordinate values is:

            * Corner orbits: :attr:`.Orbit.TETRAD_111`, :attr:`.Orbit.TETRAD_M11`
            * Edge orbits: :attr:`.Orbit.SLICE_MIDDLE`, :attr:`.Orbit.SLICE_EQUATOR`, :attr:`.Orbit.SLICE_STANDING`

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube("U F2 R'")

        Get corner coordinates.

        >>> cube.get_coord('co')   # corner orientation
        657
        >>> cube.get_coord('cp')   # corner permutation
        25253
        >>> cube.get_coord('pcp')  # partial corner permutation
        (984, 679)

        Get edge coordinates.

        >>> cube.get_coord('eo')   # edge orientation
        0
        >>> cube.get_coord('ep')   # edge permutation
        85684063
        >>> cube.get_coord('pep')  # partial edge permutation
        (8087, 7016, 3576)
        """
        if not isinstance(coord_name, str):
            raise TypeError(f"coord_name must be str, not {type(coord_name).__name__}")

        if coord_name in ("co", "eo"):
            orientation = self.orientation[:NUM_CORNERS] if coord_name == "co" else self.orientation[NUM_CORNERS:]
            return utils.get_orientation_coord(orientation, 3 if coord_name == "co" else 2, is_modulo=True)

        if coord_name in ("cp", "ep", "pcp", "pep"):
            permutation = self.permutation[:NUM_CORNERS] if coord_name in ("cp", "pcp") else self.permutation[NUM_CORNERS:]
            if coord_name in ("cp", "ep"):
                coord = utils.get_permutation_coord(permutation)
                if coord_name == "ep":
                    return coord // 2
                return coord
            orbits = CORNER_ORBITS if coord_name == "pcp" else EDGE_ORBITS
            combs = [np.where(np.array([Cubie(INDEX_TO_CUBIE[p]).orbit for p in permutation]) == orbit)[0] for orbit in orbits]
            coord = [utils.get_partial_permutation_coord(permutation[comb], comb) if len(comb) else NONE for comb in combs]
            if any(c != NONE for c in coord[1:]):
                return tuple(coord)
            return coord[0]

        raise ValueError(f"coord_name must be one of 'co', 'eo', 'cp', 'ep', 'pcp', 'pep' (got '{coord_name}')")

    def set_coord(self, coord_name: str, coord: CoordType):
        """
        Set cube coordinate value.

        Parameters
        ----------
        coord_name : {'co', 'eo', 'cp', 'ep', 'pcp', 'pep'}
            Set the specified cube coordinate.

            * 'co' means `corner orientation`.
            * 'eo' means `edge orientation`.
            * 'cp' means `corner permutation`.
            * 'ep' means `edge permutation`.
            * 'pcp' means `partial corner permutation` (a value for each corner `orbit`).
            * 'pep' means `partial edge permutation` (a value for each edge `orbit`).

        coord : int or tuple of int
            Cube coordinate value. For partial coordinate values, a value of ``-1``
            indicates an `orbit` with no permutation values (e.g., ``(-1, -1, 0)``),
            meaning the permutation (and orientation) values for that `orbit` are set to ``-1``.
            If an `int` is passed as a partial coordinate value, the value
            will be applied only to the first `orbit` (i.e., ``(coord, -1, -1)``).

            The `orbit` order for partial coordinate values is:

            * Corner orbits: :attr:`.Orbit.TETRAD_111`, :attr:`.Orbit.TETRAD_M11`
            * Edge orbits: :attr:`.Orbit.SLICE_MIDDLE`, :attr:`.Orbit.SLICE_EQUATOR`, :attr:`.Orbit.SLICE_STANDING`

        Notes
        -----
        Corner and edge permutation parities are always either both ``even`` or both ``odd``.
        This constraint is enforced when setting the `(partial) corner permutation`
        and the normal `edge permutation` coordinates by adjusting the permutation of the edges
        in the :attr:`permutation` array while preserving the same normal `edge permutation` coordinate
        (i.e. the number of valid edge permutations is effectively halved).
        When setting the `partial edge permutation` coordinate, if the edge permutation parity
        does not match the corner permutation parity, :attr:`permutation_parity` will be set to ``None``.

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube()

        Set corner coordinates.

        >>> cube.set_coord('co', 456)           # corner orientation
        >>> cube.orientation[:8]
        array([0, 1, 2, 1, 2, 2, 0, 1])
        >>> cube.set_coord('cp', 28179)         # corner permutation
        >>> cube.permutation[:8]
        array([5, 4, 0, 7, 1, 3, 6, 2])
        >>> cube.set_coord('pcp', (1273, 391))  # partial corner permutation
        >>> cube.permutation[:8]
        array([5, 4, 0, 7, 1, 3, 6, 2])

        Set edge coordinates.

        >>> cube.set_coord('eo', 673)                  # edge orientation
        >>> cube.orientation[8:]
        array([0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0])
        >>> cube.set_coord('ep', 96690777)             # edge permutation
        >>> cube.permutation[8:]
        array([12, 18, 10, 19,  9, 13, 14, 17, 16,  8, 11, 15])
        >>> cube.set_coord('pep', (7262, 2633, 8640))  # partial edge permutation
        >>> cube.permutation[8:]
        array([12, 18, 10, 19,  9, 13, 14, 17, 16,  8, 11, 15])

        Set some partial coordinates with ``-1``.

        >>> cube.set_coord('pcp', (1273, -1))      # same as cube.set_coord('pcp', 1273)
        >>> cube.permutation[:8]
        array([-1, -1,  0, -1,  1,  3, -1,  2])
        >>> cube.set_coord('pep', (7262, -1, -1))  # same as cube.set_coord('pep', 7262)
        >>> cube.permutation[8:]
        array([-1, -1, 10, -1,  9, -1, -1, -1, -1,  8, 11, -1])
        """
        if not isinstance(coord_name, str):
            raise TypeError(f"coord_name must be str, not {type(coord_name).__name__}")
        if not isinstance(coord, (int, tuple)):
            raise TypeError(f"coord must be int or tuple, not {type(coord).__name__}")

        if coord_name in ("co", "eo"):
            if not isinstance(coord, int):
                raise TypeError(f"coord must be int for coord_name '{coord_name}', not {type(coord).__name__}")
            orientation = self.orientation[:NUM_CORNERS] if coord_name == "co" else self.orientation[NUM_CORNERS:]
            v = 3 if coord_name == "co" else 2
            orientation[:] = utils.get_orientation_array(coord, v, len(orientation), force_modulo=True)
        elif coord_name in ("cp", "ep", "pcp", "pep"):
            permutation = self.permutation[:NUM_CORNERS] if coord_name in ("cp", "pcp") else self.permutation[NUM_CORNERS:]
            if coord_name in ("cp", "ep"):
                if not isinstance(coord, int):
                    raise TypeError(f"coord must be int for coord_name '{coord_name}', not {type(coord).__name__}")
                permutation[:], permutation_parity = utils.get_permutation_array(coord, len(permutation), coord_name == "ep")
                if coord_name == "ep":
                    permutation += NUM_CORNERS
                if self.permutation_parity is None:
                    other_perm = self.permutation[NUM_CORNERS:] if coord_name == "cp" else self.permutation[:NUM_CORNERS]
                    if not np.any(other_perm == CUBIE_TO_INDEX[Cubie.NONE]):
                        other_parity = utils.get_permutation_parity(other_perm)
                        self.permutation_parity = permutation_parity if permutation_parity == other_parity else other_parity
            else:
                orbits = CORNER_ORBITS if coord_name == "pcp" else EDGE_ORBITS
                if isinstance(coord, int):
                    coord_tuple = (coord,) + (NONE,) * (len(orbits) - 1)
                else:
                    coord_tuple = coord
                size = len(coord_tuple)
                if size != len(orbits):
                    raise ValueError(f"coord tuple length must be {len(orbits)} for coord_name '{coord_name}' (got {size})")
                perm = np.full_like(permutation, CUBIE_TO_INDEX[Cubie.NONE])
                for coord, orbit in zip(coord_tuple, orbits):
                    if not isinstance(coord, int):
                        raise TypeError(f"coord tuple elements must be int, not {type(coord).__name__}")
                    if coord != NONE:
                        if coord < 0 or coord >= math.perm(len(perm), NUM_ORBIT_ELEMS):
                            raise ValueError(f"coord must be >= 0 and < {math.perm(len(perm), NUM_ORBIT_ELEMS)} (got {coord})")
                        partial_permuttion, combination = utils.get_partial_permutation_array(coord, NUM_ORBIT_ELEMS)
                        if np.any(perm[combination] != CUBIE_TO_INDEX[Cubie.NONE]):
                            raise ValueError(f"invalid partial coordinates, overlapping detected (got {coord_tuple})")
                        perm[combination] = partial_permuttion + ORBIT_OFFSET[orbit]
                permutation[:] = perm
                if np.any(self.permutation == CUBIE_TO_INDEX[Cubie.NONE]):
                    self.orientation = np.where(self.permutation == CUBIE_TO_INDEX[Cubie.NONE], NONE, self.orientation)
                    self.permutation_parity = None
                    permutation_parity = None
                else:
                    corner_parity = utils.get_permutation_parity(self.permutation[:NUM_CORNERS])
                    edge_parity = utils.get_permutation_parity(self.permutation[NUM_CORNERS:])
                    permutation_parity = corner_parity
                    if corner_parity == edge_parity:
                        self.permutation_parity = corner_parity
                    else:
                        if coord_name == "pcp":
                            self.permutation_parity = edge_parity
                        else:
                            warnings.warn("invalid cube parity")
                            self.permutation_parity = None
            if self.permutation_parity is not None and self.permutation_parity != permutation_parity:
                self.permutation[-2:] = self.permutation[[-1, -2]]
                if coord_name in ("cp", "pcp"):
                    self.permutation_parity = permutation_parity
            condition = (self.permutation != CUBIE_TO_INDEX[Cubie.NONE]) & (self.orientation == NONE)
            self.orientation = np.where(condition, 0, self.orientation)
        else:
            raise ValueError(f"coord_name must be one of 'co', 'eo', 'cp', 'ep', 'pcp', 'pep' (got '{coord_name}')")

    def get_coords(
            self,
            partial_corner_perm: bool = False,
            partial_edge_perm: bool = False) -> CoordsType:
        """
        Get cube coordinates.

        Get the `corner orientation`, `edge orientation`,
        `(partial) corner permutation` and `(partial) edge permutation` coordinates.

        Parameters
        ----------
        partial_corner_perm : bool, optional
            If ``True``, returns the `partial corner permutation` coordinate,
            otherwise returns the normal `corner permutation` coordinate. Default is ``False``.
        partial_edge_perm : bool, optional
            If ``True``, returns the `partial edge permutation` coordinate,
            otherwise returns the normal `edge permutation` coordinate. Default is ``False``.

        Returns
        -------
        coords : tuple of (int or tuple of int)
            Cube coordinates in the following order:
            `corner orientation`, `edge orientation`, `(partial) corner permutation`, `(partial) edge permutation`.

        See Also
        --------
        get_coord

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube("U F2 R'")

        Get cube coordinates.

        >>> cube.get_coords()
        (657, 0, 25253, 85684063)

        Get cube coordinates with `partial corner permutation` and `partial edge permutation`.

        >>> cube.get_coords(partial_corner_perm=True, partial_edge_perm=True)
        (657, 0, (984, 679), (8087, 7016, 3576))
        """
        if not isinstance(partial_corner_perm, bool):
            raise TypeError(f"partial_corner_perm must be bool, not {type(partial_corner_perm).__name__}")
        if not isinstance(partial_edge_perm, bool):
            raise TypeError(f"partial_edge_perm must be bool, not {type(partial_edge_perm).__name__}")

        return (self.get_coord("co"), self.get_coord("eo"),
                self.get_coord("pcp" if partial_corner_perm else "cp"),
                self.get_coord("pep" if partial_edge_perm else "ep"))

    def set_coords(
            self,
            coords: CoordsType,
            partial_corner_perm: bool = False,
            partial_edge_perm: bool = False):
        """
        Set cube coordinates.

        Set the `corner orientation`, `edge orientation`,
        `(partial) corner permutation` and `(partial) edge permutation` coordinates.

        Parameters
        ----------
        coords : tuple of (int or tuple of int)
            Cube coordinates in the following order:
            `corner orientation`, `edge orientation`, `(partial) corner permutation`, `(partial) edge permutation`.
        partial_corner_perm : bool, optional
            If ``True``, sets the `partial corner permutation` coordinate,
            otherwise sets the normal `corner permutation` coordinate. Default is ``False``.
        partial_edge_perm : bool, optional
            If ``True``, sets the `partial edge permutation` coordinate,
            otherwise sets the normal `edge permutation` coordinate. Default is ``False``.

        See Also
        --------
        set_coord

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube()

        Set cube coordinates.

        >>> coords = (657, 0, 25253, 85684063)
        >>> cube.set_coords(coords)
        >>> cube
        WWBWWBYYOGGROOROOBGGWGGWRRYBRRBRROOGYOOYBBWBBWWGYYGYYR

        Set cube coordinates with `partial corner permutation` and `partial edge permutation`.

        >>> coords = (657, 0, (984, 679), (8087, 7016, 3576))
        >>> cube.set_coords(coords, partial_corner_perm=True, partial_edge_perm=True)
        >>> cube
        WWBWWBYYOGGROOROOBGGWGGWRRYBRRBRROOGYOOYBBWBBWWGYYGYYR
        """
        if not isinstance(coords, tuple):
            raise TypeError(f"coords must be tuple, not {type(coords).__name__}")
        if not isinstance(partial_corner_perm, bool):
            raise TypeError(f"partial_corner_perm must be bool, not {type(partial_corner_perm).__name__}")
        if not isinstance(partial_edge_perm, bool):
            raise TypeError(f"partial_edge_perm must be bool, not {type(partial_edge_perm).__name__}")
        if len(coords) < 4:
            raise ValueError(f"coords tuple length must be 4 (got {len(coords)})")

        self.set_coord("co", coords[0])
        self.set_coord("eo", coords[1])
        self.set_coord("pcp" if partial_corner_perm else "cp", coords[2])
        self.set_coord("pep" if partial_edge_perm else "ep", coords[3])

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
        WWBWWBYYOGGROOROOBGGWGGWRRYBRRBRROOGYOOYBBWBBWWGYYGYYR
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
    WWWWWWWWWGGGOOOOOORRRGGGGGGBBBRRRRRROOOBBBBBBYYYYYYYYY
    >>> apply_move(cube, Move.M2)   # M2 slice move
    WYWWYWWYWOOOOOOOOOGBGGBGGBGRRRRRRRRRBGBBGBBGBYWYYWYYWY
    >>> apply_move(cube, Move.FW3)  # Fw' wide move
    WWWRRRRRROWWOWWOWWGGGGGGGGGYYRYYRYYRBBBBBBBBBOOOOOOYYY
    >>> apply_move(cube, Move.X1)   # x rotation
    GGGGGGGGGOOOOOOOOOYYYYYYYYYRRRRRRRRRWWWWWWWWWBBBBBBBBB
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
    * Slice moves (e.g. `M`, `E2`, `S'`).
    * Wide moves (e.g. `Uw`, `Fw2`, `Rw'` or `u`, `f2`, `r'`).
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
    >>> apply_maneuver(cube, "U M2 Fw' x")
    RGGBBORGGWYWWYWGOOGOOGOOYWYYWYYWYRRBRRBRRBWYWBRBBGBOGO
    """
    if not isinstance(cube, Cube):
        raise TypeError(f"cube must be Cube, not {type(cube).__name__}")

    cube = cube.copy()
    cube.apply_maneuver(maneuver)
    return cube
