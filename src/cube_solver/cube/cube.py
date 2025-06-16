"""Cube module."""
import math
import warnings
import numpy as np
from itertools import chain

from .defs import NONE, SIZE, NUM_DIMS, NUM_CORNERS, NUM_EDGES, NUM_ORBIT_ELEMS
from .defs import CORNER_ORIENTATION_SIZE, EDGE_ORIENTATION_SIZE, CORNER_PERMUTATION_SIZE, EDGE_PERMUTATION_SIZE
from .enums import Axis, Orbit, Layer, Color, Face, Cubie, Move
from . import utils


def fromatwarning(message, category, *args, **kwargs):
    return f"{category.__name__}: {message}\n"


warnings.simplefilter("always")
warnings.formatwarning = fromatwarning

# TODO make all public just to check documentation


REPR_ORDER = [Face.UP, Face.LEFT, Face.FRONT, Face.RIGHT, Face.BACK, Face.DOWN]
COLORS_TYPE = [(axis.name, int) for axis in Axis.cartesian_axes()]
AXIS_ORIENTATION_ORDER = (Axis.Y, Axis.Z, Axis.X)
SWAP_COLORS = {  # swap colors along axis
    Axis.X: [Axis.X.name, Axis.Z.name, Axis.Y.name],
    Axis.Y: [Axis.Z.name, Axis.Y.name, Axis.X.name],
    Axis.Z: [Axis.Y.name, Axis.X.name, Axis.Z.name]
}
ORBIT_OFFSET = {   # TODO improve when testing differnt order of cubies
    Orbit.SLICE_MIDDLE: 8,
    Orbit.SLICE_EQUATOR: 16,
    Orbit.SLICE_STANDING: 12,
    Orbit.TETRAD_111: 0,
    Orbit.TETRAD_M11: 4
}


class Cube:
    def __init__(self, scramble: str | None = None, repr: str | None = None, random_state: bool = False):
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
        representing the colors :attr:`Color.WHITE`, :attr:`Color.GREEN`, :attr:`Color.RED`,
        :attr:`Color.YELLOW`, :attr:`Color.BLUE`, and :attr:`Color.ORANGE`, respectively.
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

        If the ``orientation`` or ``permutation`` values could not be determined correctly from the
        string representation, the orientation array and permutation array will have a ``-1`` for those locations.

        The default color scheme used for the cube is as follows (note: this may differ when using the ``repr`` parameter):

        * :attr:`Face.UP`: :attr:`Color.WHITE`
        * :attr:`Face.FRONT`: :attr:`Color.GREEN`
        * :attr:`Face.RIGHT`: :attr:`Color.RED`
        * :attr:`Face.DOWN`: :attr:`Color.YELLOW`
        * :attr:`Face.BACK`: :attr:`Color.BLUE`
        * :attr:`Face.LEFT`: :attr:`Color.ORANGE`

        Examples
        --------
        >>> from cube_solver import Cube

        Initial scramble.

        >>> cube = Cube("U F R")
        >>> cube  # string representation of the cube state
        WWRWWROORGGYOOYOOYGGBGGYGGYWWWRRBRRBGOOWBBWBBRRBYYBYYO

        Initial string representation.

        >>> cube = Cube(repr="WWRWWROORGGYOOYOOYGGBGGYGGYWWWRRBRRBGOOWBBWBBRRBYYBYYO")
        >>> print(cube)  # print a visual layout of the cube state
                ---------
                | W W R |
                | W W R |
                | O O R |
        ---------------------------------
        | G G Y | G G B | W W W | G O O |
        | O O Y | G G Y | R R B | W B B |
        | O O Y | G G Y | R R B | W B B |
        ---------------------------------
                | R R B |
                | Y Y B |
                | Y Y O |
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

        self._color_scheme: dict[Face, Color]
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
        :attr:`Face.UP`, :attr:`Face.DOWN`, :attr:`Face.RIGHT` and :attr:`Face.LEFT` face turns, it does not appear `flipped`.
        Edge orientation values are:

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
        self.permutation_parity: bool | None
        """
        Permutation parity.

        The ``permutation_parity`` indicates the parity of both `corner` and `edge` permutations
        (i.e. both parities are always the same), ``True`` for ``odd`` parity, ``False`` for ``even`` parity.,
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
    def coords(self) -> tuple[int, int, int, int]:
        """
        Cube coordinates.

        Corner orientation, edge orientation,
        corner permutation, and edge permutation coordinates.

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube("U F R")

        Get cube coordinates.

        >>> cube.coords
        (456, 673, 28179, 193381554)

        Set cube coordinates.

        >>> cube.coords = (0, 0, 0, 0)  # solved state
        >>> cube
        WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY
        """
        return self._get_coords()

    @coords.setter
    def coords(self, coords: tuple[int | tuple[int, ...], ...]):
        self._set_coords(coords)

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
            cubie_orientation = self.orientation[cubie]
            cubie_permutation = Cubie(self.permutation[cubie])
            cubie_colors = np.array(solved_colors[cubie_permutation._index], dtype=COLORS_TYPE)
            if cubie.is_corner:
                if cubie_permutation != Cubie.NONE and cubie.orbit != cubie_permutation.orbit:
                    cubie_colors = np.array(cubie_colors[SWAP_COLORS[Axis.Y]], dtype=COLORS_TYPE)
                if cubie_orientation != NONE and cubie_orientation:
                    shift = cubie_orientation if cubie.orbit == Orbit.TETRAD_M11 else -cubie_orientation
                    cubie_colors = np.array(tuple(np.roll(cubie_colors.tolist(), shift)), dtype=COLORS_TYPE)
            else:
                orbits = {cubie.orbit, cubie_permutation.orbit}
                if orbits == {Orbit.SLICE_MIDDLE, Orbit.SLICE_EQUATOR}:
                    shift = 1 if cubie.orbit == Orbit.SLICE_EQUATOR else -1
                    cubie_colors = np.array(tuple(np.roll(cubie_colors.tolist(), shift)), dtype=COLORS_TYPE)
                elif orbits == {Orbit.SLICE_MIDDLE, Orbit.SLICE_STANDING}:
                    cubie_colors = np.array(cubie_colors[SWAP_COLORS[Axis.Y]], dtype=COLORS_TYPE)
                elif orbits == {Orbit.SLICE_EQUATOR, Orbit.SLICE_STANDING}:
                    cubie_colors = np.array(cubie_colors[SWAP_COLORS[Axis.X]], dtype=COLORS_TYPE)
                if cubie_orientation != NONE and cubie_orientation:
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

        Parameters
        ----------
        repr : str
            String representation of the cube.
        """
        if len(repr) != len(REPR_ORDER)*SIZE*SIZE:
            raise ValueError(f"repr length must be {len(REPR_ORDER)*SIZE*SIZE} (got {len(repr)})")

        face_repr = [repr[i:i+SIZE*SIZE] for i in range(0, len(REPR_ORDER)*SIZE*SIZE, SIZE*SIZE)]
        for face, _repr in zip(REPR_ORDER, face_repr):
            self._colors[face._index][face.axis.name] = np.reshape([*map(Color.from_char, _repr)], shape=(SIZE, SIZE))

        # centers
        inv_color_scheme = {color: Face.NONE for color in Color.colors()}
        for center in Cubie.centers():
            self._color_scheme[Face.from_char(center.name)] = Color(self._colors[center._index][center.axis.name])
        inv_color_scheme.update({color: face for face, color in self._color_scheme.items()})

        # corners and edges
        for cubie in chain(Cubie.corners(), Cubie.edges()):
            cubie_colors = [self._colors[cubie._index][axis.name] for axis in AXIS_ORIENTATION_ORDER]
            cubie_faces = np.array([inv_color_scheme[color] for color in cubie_colors if color != Color.NONE], dtype=Face)
            try:
                if cubie.is_corner:
                    if cubie.orbit == Orbit.TETRAD_111:
                        x, z = [AXIS_ORIENTATION_ORDER.index(axis) for axis in (Axis.X, Axis.Z)]
                        cubie_faces[x], cubie_faces[z] = cubie_faces[z], cubie_faces[x]  # swap along `Y` axis
                    self.orientation[cubie] = [face.axis for face in cubie_faces].index(Axis.Y)
                else:
                    axis_importance = [AXIS_ORIENTATION_ORDER.index(face.axis) for face in cubie_faces]
                    self.orientation[cubie] = np.diff(axis_importance)[0] < 0
                self.permutation[cubie] = Cubie.from_faces(cubie_faces.tolist())
            except Exception:
                self.orientation[cubie] = NONE
                self.permutation[cubie] = Cubie.NONE

        if np.any(self.permutation == Cubie.NONE):
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
        >>> cube = Cube(random_state=True)  # doctest: +SKIP
        >>> cube.reset()  # doctest: +SKIP
        >>> cube  # doctest: +SKIP
        WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY
        """
        self._color_scheme = {face: color for face, color in zip(Face.faces(), Color.colors())}
        self._colors = np.full((SIZE,) * NUM_DIMS, (Color.NONE,) * NUM_DIMS, dtype=COLORS_TYPE)
        self.orientation = np.zeros(NUM_CORNERS + NUM_EDGES, dtype=int)
        self.permutation = np.arange(NUM_CORNERS + NUM_EDGES, dtype=int)
        self.permutation_parity = False

    def set_random_state(self):
        """
        Set a uniform random state.

        Sets random corner orientation, edge orientation,
        corner permutation, and edge permutation coordinates.
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
            The move to apply.

        Examples
        --------
        >>> from cube_solver import Cube
        >>> from cube_solver.cube import Move
        >>> cube = Cube()
        >>> cube.apply_move(Move.U1)   # U
        >>> cube.apply_move(Move.M2)   # M2
        >>> cube.apply_move(Move.FW3)  # Fw'
        >>> cube.apply_move(Move.X1)   # x
        >>> cube
        RGGBBORGGWYWWYWGOOGOOGOOYWYYWYYWYRRBRRBRRBWYWBRBBGBOGO
        """
        if not isinstance(move, Move):
            raise TypeError(f"move must be Move, not {type(move).__name__}")

        if move.is_face:
            layer = move.layers[0]
            shift = move.shifts[0]
            cubies = np.roll(layer.perm, shift, axis=1)
            orientation = self.orientation[cubies]
            if shift % 2 == 1:
                if move.axis == Axis.Z:
                    orientation = np.where(orientation != NONE, (orientation + [[1, 2, 1, 2], [1] * 4]) % ([3], [2]), NONE)
                elif move.axis == Axis.X:
                    orientation[0] = np.where(orientation[0] != NONE, (orientation[0] + [2, 1, 2, 1]) % 3, NONE)
            self.orientation[layer.perm] = orientation
            self.permutation[layer.perm] = self.permutation[cubies]

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
            orientation = self.orientation[cubies]
            permutation = self.permutation[cubies]
            if move.shifts[0] % 2 == 1:
                is_corner = np.array([cubie.is_corner for cubie in cubies])
                cubie_orbits = np.array([cubie.orbit for cubie in cubies])
                perm_orbits = np.array([Cubie(perm).orbit for perm in permutation])
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
            self.orientation[rotations] = orientation
            self.permutation[rotations] = [rotation[perm] for perm in permutation]

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

        for move_str in maneuver.split():
            self.apply_move(Move.from_string(move_str))

    def get_coord(self, coord_type: str) -> int | tuple[int, ...]:
        """
        Get cube coordinate value.

        Parameters
        ----------
        coord_type : {'co', 'eo', 'cp', 'ep', 'pcp', 'pep'}
            Get the specified cube coordinate.

            * 'co' means corner orientation.
            * 'eo' means edge orientation.
            * 'cp' means corner permutation.
            * 'ep' means edge permutation.
            * 'pcp' means partial corner permutation (a value for each corner `orbit`).
            * 'pep' means partial edge permutation (a value for each edge `orbit`).

        Returns
        -------
        coord : int or tuple of int
            Cube coordinate value. For partial coordinate values, a value of ``-1``
            indicates an `orbit` with no permutation values (e.g., ``(-1, -1, 0)``),
            meaning the permutation values for that `orbit` are set to ``-1``.
            If only the value of the first `orbit` is available (i.e., ``(coord, -1, -1)``),
            returns an `int` representing that partial coordinate value.

            The `orbit` order for partial coordinate values is:

            * Corner orbits: :attr:`Orbit.TETRAD_111`, :attr:`Orbit.TETRAD_M11`
            * Edge orbits: :attr:`Orbit.SLICE_MIDDLE`, :attr:`Orbit.SLICE_EQUATOR`, :attr:`Orbit.SLICE_STANDING`

        See Also
        --------
        set_coord

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube("U F R")

        Get corner coordinates.

        >>> cube.get_coord('co')
        456
        >>> cube.get_coord('cp')
        28179
        >>> cube.get_coord('pcp')
        (1273, 391)

        Get edge coordinates.

        >>> cube.get_coord('eo')
        673
        >>> cube.get_coord('ep')
        193381554
        >>> cube.get_coord('pep')
        (2633, 8640, 7262)
        """
        if not isinstance(coord_type, str):
            raise TypeError(f"coord_type must be str, not {type(coord_type).__name__}")

        if coord_type in ("co", "eo"):
            orientation = self.orientation[:NUM_CORNERS] if coord_type == "co" else self.orientation[NUM_CORNERS:]
            return utils.get_orientation_coord(orientation, 3 if coord_type == "co" else 2, is_modulo=True)

        if coord_type in ("cp", "ep", "pcp", "pep"):
            permutation = self.permutation[:NUM_CORNERS] if coord_type in ("cp", "pcp") else self.permutation[NUM_CORNERS:]
            if coord_type in ("cp", "ep"):
                coord = utils.get_permutation_coord(permutation)
                if coord_type == "ep":
                    return coord // 2
                return coord
            orbits = [*Orbit.tetrads()] if coord_type == "pcp" else [*Orbit.slices()]
            combs = [np.where(np.array([Cubie(perm).orbit for perm in permutation]) == orbit)[0] for orbit in orbits]
            coord = tuple(utils.get_partial_permutation_coord(permutation[comb], comb) if len(comb) else -1 for comb in combs)
            if any(c != NONE for c in coord[1:]):
                return coord
            return coord[0]

        raise ValueError(f"coord_type must be one of 'co', 'eo', 'cp', 'ep', 'pcp', 'pep' (got '{coord_type}')")

    def _get_coords(
            self,
            partial_corner_perm: bool = False,
            partial_edge_perm: bool = False) -> tuple[int | tuple[int, ...], ...]:
        """
        Get cube coordinates.

        Get the corner orientation, edge orientation,
        (partial) corner permutation and (partial) edge permutation coordinates.

        Parameters
        ----------
        partial_corner_perm : bool, optional
            If ``True``, returns the partial corner permutation coordinate,
            otherwise returns the normal corner permutation coordinate. Default is ``False``.
        partial_edge_perm : bool, optional
            If ``True``, returns the partial edge permutation coordinate,
            otherwise returns the normal edge permutation coordinate. Default is ``False``.

        Returns
        -------
        coords : tuple of (int or tuple of int)
            Cube coordinates in the following order:
            corner orientation, edge orientation, (partial) corner permutation, (partial) edge permutation.

        See Also
        --------
        get_coord

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube("U F R")

        Get cube coordinates.

        >>> cube.get_coords()
        (456, 673, 28179, 193381554)

        Get cube coordinates with partial corner permutation and partial edge permutation.

        >>> cube.get_coords(partial_corner_perm=True, partial_edge_perm=True)
        (456, 673, (1273, 391), (2633, 8640, 7262))
        """
        # print("getter")
        return (self.get_coord("co"), self.get_coord("eo"),
                self.get_coord("pcp" if partial_corner_perm else "cp"),
                self.get_coord("pep" if partial_edge_perm else "ep"))

    def set_coord(self, coord_type: str, coord: int | tuple[int, ...]):
        """
        Set cube coordinate value.

        Parameters
        ----------
        coord_type : {'co', 'eo', 'cp', 'ep', 'pcp', 'pep'}
            Set the specified cube coordinate.

            * 'co' means corner orientation.
            * 'eo' means edge orientation.
            * 'cp' means corner permutation.
            * 'ep' means edge permutation.
            * 'pcp' means partial corner permutation (a value for each corner `orbit`).
            * 'pep' means partial edge permutation (a value for each edge `orbit`).

        coord : int or tuple of int
            Cube coordinate value. For partial coordinate values, a value of ``-1``
            indicates an `orbit` with no permutation values (e.g., ``(-1, -1, 0)``),
            meaning the permutation values for that `orbit` are set to ``-1``.
            If an `int` is passed as a partial coordinate value, the value
            will be applied only to the first `orbit` (i.e., ``(coord, -1, -1)``).

            The `orbit` order for partial coordinate values is:

            * Corner orbits: :attr:`Orbit.TETRAD_111`, :attr:`Orbit.TETRAD_M11`
            * Edge orbits: :attr:`Orbit.SLICE_MIDDLE`, :attr:`Orbit.SLICE_EQUATOR`, :attr:`Orbit.SLICE_STANDING`

        # Notes
        # -----
        # Corner and edge permutation parities are always either both `odd` or both `even`.
        # This constraint is enforced when setting the normal corner and edge permutation coordinates
        # by modifying the edge permutation accordingly to ensure both parities match (i.e., the number of valid
        # edge permutations is halved). For this reason, it is recommended to set the corner permutation coordinate
        # before the edge permutation coordinate. The parity constraint is not enforced when setting
        # partial permutation coordinates for either corners or edges.

        See Also
        --------
        get_coord

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube()

        Set corner coordinates.

        >>> cube.set_coord('co', 456)
        >>> cube.orientation[:8]
        array([0, 1, 2, 1, 2, 2, 0, 1])
        >>> cube.set_coord('cp', 28179)
        >>> cube.permutation[:8]
        array([5, 4, 0, 7, 1, 3, 6, 2])
        >>> cube.set_coord('pcp', (1273, 391))
        >>> cube.permutation[:8]
        array([5, 4, 0, 7, 1, 3, 6, 2])

        Set edge coordinates.

        >>> cube.set_coord('eo', 673)
        >>> cube.orientation[8:]
        array([0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0])
        >>> cube.set_coord('ep', 193381554)
        >>> cube.permutation[8:]
        array([12, 18, 10, 19,  9, 13, 14, 17, 16,  8, 11, 15])
        >>> cube.set_coord('pep', (2633, 8640, 7262))
        >>> cube.permutation[8:]
        array([12, 18, 10, 19,  9, 13, 14, 17, 16,  8, 11, 15])

        Set some partial coordinates with ``-1``.

        >>> cube.set_coord('pcp', (1273, -1))      # same as cube.set_coord('pcp', 1273)
        >>> cube.permutation[:8]
        array([-1, -1,  0, -1,  1,  3, -1,  2])
        >>> cube.set_coord('pep', (2633, -1, -1))  # same as cube.set_coord('pep', 2633)
        >>> cube.permutation[8:]
        array([-1, 18, -1, 19, -1, -1, -1, 17, 16, -1, -1, -1])
        """
        if not isinstance(coord_type, str):  # TODO add that orientation is also not enforced for partial permutation
            raise TypeError(f"coord_type must be str, not {type(coord_type).__name__}")
        if not isinstance(coord, (int, tuple)):
            raise TypeError(f"coord must be int or tuple, not {type(coord).__name__}")

        if coord_type in ("co", "eo"):
            if not isinstance(coord, int):
                raise TypeError(f"coord must be int for coord_type '{coord_type}', not {type(coord).__name__}")
            orientation = self.orientation[:NUM_CORNERS] if coord_type == "co" else self.orientation[NUM_CORNERS:]
            v = 3 if coord_type == "co" else 2
            orientation[:] = utils.get_orientation_array(coord, v, len(orientation), force_modulo=True)
        elif coord_type in ("cp", "ep", "pcp", "pep"):
            permutation = self.permutation[:NUM_CORNERS] if coord_type in ("cp", "pcp") else self.permutation[NUM_CORNERS:]
            if coord_type in ("cp", "ep"):
                if not isinstance(coord, int):
                    raise TypeError(f"coord must be int for coord_type '{coord_type}', not {type(coord).__name__}")
                permutation[:], permutation_parity = utils.get_permutation_array(coord, len(permutation), coord_type == "ep")
                if coord_type == "ep":
                    permutation += NUM_CORNERS
                if self.permutation_parity is None:
                    other_perm = self.permutation[NUM_CORNERS:] if coord_type == "cp" else self.permutation[:NUM_CORNERS]
                    if not np.any(other_perm == Cubie.NONE):
                        other_parity = utils.get_permutation_parity(other_perm)
                        self.permutation_parity = permutation_parity if permutation_parity == other_parity else other_parity
            else:
                orbits = [*Orbit.tetrads()] if coord_type == "pcp" else [*Orbit.slices()]
                if isinstance(coord, int):
                    coord_tuple = (coord,) + (NONE,) * (len(orbits) - 1)
                else:
                    coord_tuple = coord
                if len(coord_tuple) != len(orbits):
                    raise ValueError(f"coord tuple length must be {len(orbits)} for coord_type '{coord_type}' (got {len(coord_tuple)})")
                perm = np.full_like(permutation, Cubie.NONE)
                for coord, orbit in zip(coord_tuple, orbits):
                    if not isinstance(coord, int):
                        raise TypeError(f"coord tuple elements must be int, not {type(coord).__name__}")
                    if coord != NONE:
                        if coord < 0 or coord >= math.perm(len(perm), NUM_ORBIT_ELEMS):
                            raise ValueError(f"coord must be >= 0 and < {math.perm(len(perm), NUM_ORBIT_ELEMS)} (got {coord})")
                        partial_permuttion, combination = utils.get_partial_permutation_array(coord, NUM_ORBIT_ELEMS)
                        if np.any(perm[combination] != Cubie.NONE):
                            raise ValueError(f"invalid partial coordinates, overlapping detected (got {coord_tuple})")
                        perm[combination] = partial_permuttion + ORBIT_OFFSET[orbit]
                permutation[:] = perm
                if np.any(self.permutation == Cubie.NONE):
                    self.orientation = np.where(self.permutation == Cubie.NONE, NONE, self.orientation)  # TODO document
                    self.permutation_parity = None
                    permutation_parity = None
                else:
                    corner_parity = utils.get_permutation_parity(self.permutation[:NUM_CORNERS])
                    edge_parity = utils.get_permutation_parity(self.permutation[NUM_CORNERS:])
                    permutation_parity = corner_parity
                    if corner_parity == edge_parity:
                        self.permutation_parity = corner_parity
                    else:
                        if coord_type == "pcp":
                            self.permutation_parity = edge_parity  # TODO document that first corner
                        else:
                            warnings.warn("invalid cube parity")
                            self.permutation_parity = None
            if self.permutation_parity is not None and self.permutation_parity != permutation_parity:
                self.permutation[-2:] = self.permutation[[-1, -2]]
                if coord_type in ("cp", "pcp"):
                    self.permutation_parity = permutation_parity
            self.orientation = np.where((self.permutation != Cubie.NONE) & (self.orientation == NONE), 0, self.orientation)
        else:
            raise ValueError(f"coord_type must be one of 'co', 'eo', 'cp', 'ep', 'pcp', 'pep' (got '{coord_type}')")

    # TODO test
    def _set_coords(
            self,
            coords: tuple[int | tuple[int, ...], ...],
            partial_corner_perm: bool = False,
            partial_edge_perm: bool = False):
        """
        Set cube coordinates.

        Set the corner orientation, edge orientation,
        (partial) corner permutation and (partial) edge permutation coordinates.

        Parameters
        ----------
        coords : tuple of (int or tuple of int)
            Cube coordinates in the following order:
            corner orientation, edge orientation, (partial) corner permutation, (partial) edge permutation.
        partial_corner_perm : bool, optional
            If ``True``, sets the partial corner permutation coordinate,
            otherwise sets the normal corner permutation coordinate. Default is ``False``.
        partial_edge_perm : bool, optional
            If ``True``, sets the partial edge permutation coordinate,
            otherwise sets the normal edge permutation coordinate. Default is ``False``.

        See Also
        --------
        set_coord

        Examples
        --------
        >>> from cube_solver import Cube
        >>> cube = Cube()

        Set cube coordinates.

        >>> coords = (456, 673, 28179, 193381554)
        >>> cube.set_coords(coords)
        >>> cube.orientation
        array([0, 0, 1, 2, 1, 2, 2, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1])
        >>> cube.permutation
        array([ 5,  4,  0,  7,  1,  3,  6,  2, 12, 18, 10, 19,  9, 13, 14, 17, 16,
                8, 11, 15])

        Set cube coordinates with partial corner permutation and partial edge permutation.

        >>> coords = (456, 673, (1273, 391), (2633, 8640, 7262))
        >>> cube.set_coords(coords, partial_corner_perm=True, partial_edge_perm=True)
        >>> cube.orientation
        array([0, 0, 1, 2, 1, 2, 2, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1])
        >>> cube.permutation
        array([ 5,  4,  0,  7,  1,  3,  6,  2, 12, 18, 10, 19,  9, 13, 14, 17, 16,
                8, 11, 15])
        """
        # print("setter")
        # if not isinstance(coords, tuple):
        #     raise TypeError(f"coords must be tuple, not {type(coords).__name__}")
        # if not isinstance(partial_corner_perm, bool):
        #     raise TypeError(f"partial_corner_perm must be bool, not {type(partial_corner_perm).__name__}")
        # if not isinstance(partial_edge_perm, bool):
        #     raise TypeError(f"partial_edge_perm must be bool, not {type(partial_edge_perm).__name__}")

        self.set_coord("co", coords[0])
        self.set_coord("eo", coords[1])
        self.set_coord("pcp" if partial_corner_perm else "cp", coords[2])
        self.set_coord("pep" if partial_edge_perm else "ep", coords[3])


# # def generate_scramble(length: int = 25) -> str:
# #     assert length >= 1, "scramble length must be greater or equal than 1"

# #     scramble = []
# #     move = None
# #     for _ in range(length):
# #         move = np.random.choice(NEXT_MOVES[move])
# #         scramble.append(move)

# #     return " ".join(scramble)


# # def apply_move(cube: Cube, move: Move) -> Cube:  # TODO make move int
# #     """
# #     Return a copy of the cube with the move applied.

# #     Parameters
# #     ----------
# #     cube : Cube
# #         Cube object.
# #     move
# #         The move to apply.

# #     Returns
# #     -------
# #     Cube
# #         Cube object with move applied.
# #     """
# #     cube = cube.copy()
# #     cube.apply_move(move)
# #     return cube


# # def apply_maneuver(cube: Cube, maneuver: str) -> Cube:
# #     """
# #     Return a copy of the cube with a sequence of moves applied.

# #     Parameters
# #     ----------
# #     cube : Cube
# #         Cube object.
# #     maneuver
# #         The sequence of moves to apply.

# #     Returns
# #     -------
# #     Cube
# #         Cube object with sequence of moves applied.

# #     Examples
# #     --------
# #     >>> import cube_solver
# #     >>> from cube_solver import Cube
# #     >>> cube = Cube()
# #     >>> cube_solver.apply_maneuver(cube, "U F R")
# #     WWRWWROORGGYOOYOOYGGBGGYGGYWWWRRBRRBGOOWBBWBBRRBYYBYYO
# #     >>> cube  # the original cube remains unchanged
# #     WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY
# #     """
# #     cube = cube.copy()
# #     cube.apply_maneuver(maneuver)
# #     return cube
