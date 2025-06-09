"""Cube module."""
import numpy as np
# # from copy import deepcopy
from itertools import chain

# from cube_solver.cube import utils
from cube_solver.cube.enums import Axis, Color, Face, Cubie, Move


# TODO move to defs
SIZE = 3
NUM_CORNERS = 8
NUM_EDGES = 12

SWAP_CUBIE = {  # swap cubie along axis
    Axis.Y: [Axis.Y, Axis.X, Axis.Z],
    Axis.Z: [Axis.X, Axis.Z, Axis.Y],
    Axis.X: [Axis.Z, Axis.Y, Axis.X]
}
REPR_ORDER = [Face.UP, Face.LEFT, Face.FRONT, Face.RIGHT, Face.BACK, Face.DOWN]

# FACTORIAL = np.cumprod([1] + list(range(1, NUM_EDGES + 1)))
# CORNER_ORIENTATION_SIZE = 3 ** (NUM_CORNERS - 1)
# EDGE_ORIENTATION_SIZE = 2 ** (NUM_EDGES - 1)
# CORNER_PERMUTATION_SIZE = FACTORIAL[NUM_CORNERS]
# EDGE_PERMUTATION_SIZE = FACTORIAL[NUM_EDGES] // 2


class Cube:
    def __init__(self, scramble: str | None = None, repr: str | None = None, random_state: bool = False):
        """
        Create `Cube` object.

        Parameters
        ----------
        scramble : str or None, optional
            Initial scramble. If `None`, no scramble is applied.
        repr : str or None, optional
            Cube string representation. If not `None`, the `scramble` parameter is ignored and
            initializes the cube with the given string representation. See `Notes` for the representation format.
        random_state : bool, optional
            If `True`, the `scramble` and `repr` parameters are ignored and
            creates a `Cube` object with a uniform random state.

        Notes
        -----
        The `repr` parameter must contain characters from `{'W', 'G', 'R', 'Y', 'B', 'O'}`,
        representing the colors :attr:`Color.WHITE`, :attr:`Color.GREEN`, :attr:`Color.RED`,
        :attr:`Color.YELLOW`, :attr:`Color.BLUE`, and :attr:`Color.ORANGE`, respectively.
        The order of the `repr` parameter is::

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

        The default color scheme used for the `scramble` and `random_state` parameters is as follows
        (note: this may differ for the `repr` parameter):

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
        (0, 0, 0, 0)
        """
        if scramble is not None and not isinstance(scramble, str):
            raise TypeError(f"scramble must be str or None, not {type(scramble).__name__}")
        if repr is not None and not isinstance(repr, str):
            raise TypeError(f"repr must be str or None, not {type(repr).__name__}")
        if not isinstance(random_state, bool):
            raise TypeError(f"random_state must be bool, not {type(random_state).__name__}")

        # self._state: str
        # """Cached string representation."""
        # self._coords: tuple
        # """Cached cube coordinates."""
        self._cubies: np.ndarray
        """
        Cubie representation of the cube.
        Used to generate and parse the string representation.
        """
        self._scheme: dict[Face, Color]
        """
        Color shceme of the cube.
        Used to generate and parse the string representation.
        """
        self.orientation: np.ndarray
        """
        Orientation array.

        The `orientation` array contains the orientation values of the `8` corners and `12` edges of the cube.
        The first `8` elements represent the `corner` orientation values, and the remaining `12` elements represent the
        `edge` orientation values.

        A corner is correctly oriented when the `top` or `bottom` facelet of the corner piece matches eihter the `top` or
        `bottom` color of the cube. Corner orientation values are:

        * `0` if the corner is `correctly` oriented.
        * `1` if the corner is `twisted clockwise` relative to the correct orientation.
        * `2` if the corner is `twisted counter-clockwise` relative to the correct orientation.

        An edge is correctly oriented if, when placed in its correct position using only
        `U`, `D`, `R` and `L` face turns, it does not appear `flipped`. Edge orientation values are:

        * `0` if the edge is `correctly` oriented.
        * `1` if the edge is incorrectly oriented (i.e. `flipped`).
        """
        self.permutation: np.ndarray
        """
        Permutation array.

        The `permutation` array contains the permutation values of the `8` corners and `12` edges of the cube.
        The first `8` elements represent the `corner` permutation values, and the remaining `12` elements represent the
        `edge` permutation values.

        The `solved state` permutation goes from `0` to `7` for the corners and from `8` to `19`
        for the edges. The piece ordering in the solved state is:

        * Corners: [`UBL`, `UFR`, `DBR`, `DFL`, `UBR`, `UFL`, `DBL`, `DFR`]
        * Edges: [`UB`, `UF`, `DB`, `DF`, `UL`, `UR`, `DL`, `DR`, `BL`, `BR`, `FL`, `FR`]
        """
        self.reset()
        if random_state:
            self.set_random_state()
        elif repr is not None:
            self._parse_repr(repr)
        elif scramble is not None:
            self.apply_maneuver(scramble)

#     # @property
#     # def coords(self) -> tuple:
#     #     """
#     #     Cube coordinates.

#     #     Corner orientation, edge orientation,
#     #     corner permutation, and edge permutation coordinates.

#     #     Returns
#     #     -------
#     #     coords : tuple
#     #         Cube coordinates.

#     #     See Also
#     #     --------
#     #     get_coords()
#     #     set_coords()

#     #     Examples
#     #     --------
#     #     >>> from cube_solver import Cube
#     #     >>> cube = Cube("U F R")

#     #     Get cube coordinates.

#     #     >>> cube.coords
#     #     (456, 673, 28179, 193381554)

#     #     Set cube coordinates.

#     #     >>> cube.coords = (0, 0, 0, 0)
#     #     >>> cube.orientation
#     #     array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
#     #     >>> cube.permutation
#     #     array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
#     #            17, 18, 19])
#     #     """
#     #     if self._coords is None:
#     #         self._coords = self.get_coords()
#     #     return self._coords

#     # @coords.setter
#     # def coords(self, coords: tuple):
#     #     self._coords = coords
#     #     self.set_coords(self._coords)

# #     def copy(self) -> "Cube":
# #         """Return a copy of the Cube object."""
# #         return deepcopy(self)

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
        # self._state = None
        # self._coords = None
        self._cubies = np.full((SIZE,) * 3 + (3,), Color.NONE, dtype=int)
        self._scheme = {face: color for face, color in zip(Face.faces(), Color.colors())}
        self.orientation = np.zeros(NUM_CORNERS + NUM_EDGES, dtype=int)
        self.permutation = np.arange(NUM_CORNERS + NUM_EDGES, dtype=int)

    def _parse_repr(self, repr: str):
        """
        Parse a string representation into a cube state.

        Parameters
        ----------
        repr : str
            String representation of the cube.
        """
        if not isinstance(repr, str):
            raise TypeError(f"repr must be str, not {type(repr).__name__}")
        if len(repr) != len(REPR_ORDER)*SIZE*SIZE:
            raise ValueError(f"repr length must be {len(REPR_ORDER)*SIZE*SIZE}, got {len(repr)}")
        if "N" in repr:
            raise ValueError("invalid color character, got 'N'")

        face_repr = [repr[i:i+SIZE*SIZE] for i in range(0, len(REPR_ORDER)*SIZE*SIZE, SIZE*SIZE)]
        for face, _repr in zip(REPR_ORDER, face_repr):
            self._cubies[face._cubie_slice] = np.reshape([*map(Color.from_char, _repr)], shape=(SIZE, SIZE))

        # centers
        self._scheme = {}
        for center in Cubie.centers():
            self._scheme[Face.from_char(center.name)] = Color(self._cubies[center._index][center.axis])
        inv_scheme = {val: key for key, val in self._scheme.items()}

        if len(set(self._scheme.values())) != len([*Color.colors()]):
            raise ValueError(f"invalid cubie centers, got {[*self._scheme.values()]}")

        # corners and edges
        for slot in chain(Cubie.corners(), Cubie.edges()):
            faces = np.array([inv_scheme[color] for color in self._cubies[slot._index] if color != Color.NONE], dtype=Color)
            if slot.is_corner:
                if slot.axis == Axis.DIAG_M11:
                    faces = faces[SWAP_CUBIE[Axis.Y]]
                try:
                    self.orientation[slot] = [face.axis for face in faces].index(Axis.Y)
                except ValueError:
                    raise ValueError(f"invalid cubie faces, got {faces}")
            else:
                self.orientation[slot] = faces[0].axis > faces[1].axis
            self.permutation[slot] = Cubie.from_faces([*faces])

        if len(set(self.permutation)) != NUM_CORNERS + NUM_EDGES:
            raise ValueError(f"invalid permutation, got {self.permutation}")
        if np.sum(self.orientation[:NUM_CORNERS]) % 3 != 0 or np.sum(self.orientation[NUM_CORNERS:]) % 2 != 0:
            raise ValueError(f"invalid orientation, got {self.orientation}")

    def __repr__(self) -> str:
        """String representation of the `Cube` object."""
        cubies = np.full_like(self._cubies, Color.NONE)
        for face in Face.faces():
            cubies[face._cubie_slice] = self._scheme[face]

        # centers
        for slot in Cubie.centers():
            self._cubies[slot._index] = cubies[slot._index]

        # corners and edges
        for slot in chain(Cubie.corners(), Cubie.edges()):  # TODO what happens when permutation has -1?
            orientation = self.orientation[slot]
            permutation = Cubie(self.permutation[slot])
            cubie = cubies[permutation._index]
            if slot.is_corner:
                if slot.axis != permutation.axis:
                    cubie = cubie[SWAP_CUBIE[Axis.Y]]
                if orientation:
                    cubie = np.roll(cubie, orientation if slot.axis == Axis.DIAG_111 else -orientation)
            else:
                axes = {slot.axis, permutation.axis}
                if axes == {0, 2}:
                    cubie = np.roll(cubie, 1 if slot.axis != 2 else -1)
                elif axes == {1, 2}:
                    cubie = cubie[SWAP_CUBIE[Axis.Y]]
                elif axes == {0, 1}:
                    cubie = cubie[SWAP_CUBIE[Axis.X]]
                if orientation:
                    cubie = cubie[SWAP_CUBIE[slot.axis]]
            self._cubies[slot._index] = cubie

        return "".join([Color(color).char for color in np.ravel([self._cubies[face._cubie_slice] for face in REPR_ORDER])])

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
                if layer.axis == Axis.Z:
                    orientation = (orientation + [[1, 2, 1, 2], [1, 1, 1, 1]]) % ([3], [2])  # corners and edges
                elif layer.axis == Axis.X:
                    orientation[0] = (orientation[0] + [2, 1, 2, 1]) % 3  # corners
            self.orientation[layer.perm] = orientation
            self.permutation[layer.perm] = self.permutation[cubies]

        elif move.is_slice:
            shift = move.shifts[0]
            index = 2 if move.axis == Axis.Z else 0
            base_move = Move[move.axis.name + "1"].layers[index].name
            self.apply_move(Move[base_move + "W" + str(-shift % 4)])
            self.apply_move(Move[base_move + str(shift % 4)])

        elif move.is_wide:
            shift = move.shifts[0]
            mult = 1 if Move[move.axis.name + "1"].layers[0].name == move.name[0] else -1
            self.apply_move(Move[move.axis.name + str(mult * shift % 4)])
            self.apply_move(Move[Face.from_char(move.name[0]).opposite.char + str(shift % 4)])

        elif move.is_rotation:
            layers_perm = [layer.perm for layer in move.layers]
            layers_shifted = [np.roll(layer, shift, axis=1) for layer, shift in zip(layers_perm, move.shifts)]
            rotation = {center: center for center in Cubie.centers()}
            rotation.update({Cubie(key): Cubie(val) for key, val in zip(np.ravel(layers_shifted), np.ravel(layers_perm))})

            # centers
            scheme = self._scheme.copy()
            for center in Cubie.centers():
                self._scheme[Face.from_char(rotation[center].name)] = scheme[Face.from_char(center.name)]

            # corners and edges
            cubies = [*Cubie.corners()] + [*Cubie.edges()]
            slots = [rotation[cubie] for cubie in cubies]
            orientation = self.orientation[cubies]
            permutation = self.permutation[cubies]
            if move.shifts[0] % 2 == 1:
                is_corner = np.array([cubie.is_corner for cubie in cubies])
                cubie_axis = np.array([cubie.axis for cubie in cubies])
                perm_axis = np.array([Cubie(perm).axis for perm in permutation])
                corner_comp = edge_comp = [perm_axis, cubie_axis]
                if move.axis == Axis.X:
                    corner_comp = [cubie_axis, Axis.DIAG_111]
                    edge_comp = [Axis.X, Axis.X]
                elif move.axis == Axis.Y:
                    edge_comp = [Axis.Y, Axis.Y]
                elif move.axis == Axis.Z:
                    corner_comp = [cubie_axis, Axis.DIAG_M11]
                axis_comp = perm_axis != np.where(is_corner, corner_comp[0], edge_comp[0])
                condition = cubie_axis == np.where(is_corner, corner_comp[1], edge_comp[1])
                incr = np.where(condition, axis_comp, np.where(is_corner, -axis_comp.astype(int), ~axis_comp))
                orientation = (orientation + incr) % np.where(is_corner, 3, 2)
            self.orientation[slots] = orientation
            self.permutation[slots] = [rotation[perm] for perm in permutation]

        else:
            raise ValueError(f"invalid move, got {repr(move)}")

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
            self.apply_move(Move.from_str(move_str))

# #     def get_coord(self, coord_type: str) -> int | tuple:
# #         """
# #         Get cube coordinate value.

# #         Parameters
# #         ----------
# #         coord_type : {'co', 'eo', 'cp', 'ep', 'pcp', 'pep'}
# #             Get the specified cube coordinate.

# #             * 'co' means corner orientation.
# #             * 'eo' means edge orientation.
# #             * 'cp' means corner permutation.
# #             * 'ep' means edge permutation.
# #             * 'pcp' means partial corner permutation.
# #             * 'pep' means partial edge permutation.

# #         Returns
# #         -------
# #         coord : int or tuple
# #             Cube coordinate value. For partial coordinate values, a value of `-1`
# #             indicates an axis with no permutation values (e.g., `(-1, -1, 0)`),
# #             meaning the permutation values for these axes are set to `-1`.
# #             If only the value of the first axis is available (i.e., `(coord, -1, -1)`),
# #             returns an `int` representing that partial coordinate value.

# #         See Also
# #         --------
# #         utils.get_orientation_coord
# #         utils.get_permutation_coord
# #         utils.get_combination_coord
# #         utils.get_partial_permutation_coord

# #         Examples
# #         --------
# #         >>> from cube_solver import Cube
# #         >>> cube = Cube("U F R")

# #         Get corner coordinates.

# #         >>> cube.get_coord('co')
# #         456
# #         >>> cube.get_coord('cp')
# #         28179
# #         >>> cube.get_coord('pcp')
# #         (1273, 391)

# #         Get edge cordinates.

# #         >>> cube.get_coord('eo')
# #         673
# #         >>> cube.get_coord('ep')
# #         193381554
# #         >>> cube.get_coord('pep')
# #         (2633, 8640, 7262)
# #         """
# #         if not isinstance(coord_type, str):
# #             raise TypeError(f"coord_type must be str, not {type(coord_type).__name__}")

# #         if coord_type in ("co", "eo"):
# #             orientation = self.orientation[:NUM_CORNERS-1] if coord_type == "co" else self.orientation[NUM_CORNERS:-1]
# #             return utils.get_orientation_coord(orientation, 3 if coord_type == "co" else 2)

# #         if coord_type in ("cp", "ep", "pcp", "pep"):
# #             permutation = self.permutation[:NUM_CORNERS] if coord_type in ("cp", "pcp") else self.permutation[NUM_CORNERS:]
# #             if coord_type in ("cp", "ep"):
# #                 return utils.get_permutation_coord(permutation, coord_type == "ep")
# #             num_axes = 2 if coord_type == "pcp" else 3
# #             combinations = [np.where(np.array([AXIS[perm] for perm in permutation]) ==axis)[0] for axis in range(num_axes)]
# #             coord = tuple(utils.get_partial_permutation_coord(permutation, combination) for combination in combinations)
# #             if np.all([c == EMPTY for c in coord[1:]]):
# #                 return coord[0]
# #             return coord

# #         raise ValueError(f"coord_type must be one of 'co', 'eo', 'cp', 'ep', 'pcp', or 'pep' (got '{coord_type}')")

# #     def get_coords(self, partial_corner_perm: bool = False, partial_edge_perm: bool = False) -> tuple:
# #         """
# #         Get cube coordinates.

# #         Get the corner orientation, edge orientation,
# #         (partial) corner permutation and (partial) edge permutation coordinates.

# #         The permutation coordinates depend on the instance attributes
# #         `partial_corner_perm` and `partial_edge_perm`.

# #         Returns
# #         -------
# #         coords : tuple
# #             Cube coordinates in the following order:
# #             corner orientation, edge orientation, (partial) corner permutation, (partial) edge permutation.
# #         partial_corner_perm : bool
# #             If `True`,use the partial corner permutation coordinate,otherwise use the normal corner permutation coordinate.
# #         partial_edge_perm : bool
# #             If `True`, use the partial edge permutation coordinate, otherwise use the normal edge permutation coordinate.

# #         See Also
# #         --------
# #         utils.get_orientation_coord
# #         utils.get_permutation_coord
# #         utils.get_combination_coord
# #         utils.get_partial_permutation_coord

# #         Examples
# #         --------
# #         >>> from cube_solver import Cube
# #         >>> cube = Cube("U F R")

# #         Get cube coordinates.

# #         >>> cube.get_coords()
# #         (456, 673, 28179, 193381554)

# #         Get cube coordinates with partial corner permutation and partial edge permutation.

# #         >>> cube.partial_corner_perm = True
# #         >>> cube.partial_edge_perm = True
# #         >>> cube.get_coords()
# #         (456, 673, (1273, 391), (2633, 8640, 7262))
# #         """
# #         # print("getter")
# #         return (self.get_coord("co"), self.get_coord("eo"),
# #                 self.get_coord("pcp" if partial_corner_perm else "cp"),
# #                 self.get_coord("pep" if partial_edge_perm else "ep"))

#     def set_coord(self, coord_type: str, coord: int | tuple):
#         """
#         Set cube coordinate value.

#         Parameters
#         ----------
#         coord_type : {'co', 'eo', 'cp', 'ep', 'pcp', 'pep'}
#             Set the specified cube coordinate.

#             * 'co' means corner orientation.
#             * 'eo' means edge orientation.
#             * 'cp' means corner permutation.
#             * 'ep' means edge permutation.
#             * 'pcp' means partial corner permutation.
#             * 'pep' means partial edge permutation.
#         coord : int or tuple
#             Coordinate value or tuple of partial coordinate values.
#             For partial coordinate values, you may pass `-1`
#             to indicate an axis that should be ignored (e.g., `(0, -1, -1)`),
#             the permutation values for these axes will be set to `-1`.
#             If you pass an `int` as a partial coordinate value,
#             the value will be applied only to the first axis (i.e., `(coord, -1, -1)`).

#         See Also
#         --------
#         utils.get_orientation_array
#         utils.get_permutation_array
#         utils.get_combination_array
#         utils.get_partial_permutation_array

#         Examples
#         --------
#         >>> from cube_solver import Cube
#         >>> cube = Cube()

#         Set corner coordinates.

#         >>> cube.set_coord('co', 456)
#         >>> cube.orientation[:8]
#         array([0, 1, 2, 1, 2, 2, 0, 1])
#         >>> cube.set_coord('cp', 28179)
#         >>> cube.permutation[:8]
#         array([5, 4, 0, 7, 1, 3, 6, 2])
#         >>> cube.set_coord('pcp', (1273, 391))
#         >>> cube.permutation[:8]
#         array([5, 4, 0, 7, 1, 3, 6, 2])

#         Set edge coordinates.

#         >>> cube.set_coord('eo', 673)
#         >>> cube.orientation[8:]
#         array([0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0])
#         >>> cube.set_coord('ep', 193381554)
#         >>> cube.permutation[8:]
#         array([12, 18, 10, 19,  9, 13, 14, 17, 16,  8, 11, 15])
#         >>> cube.set_coord('pep', (2633, 8640, 7262))
#         >>> cube.permutation[8:]
#         array([12, 18, 10, 19,  9, 13, 14, 17, 16,  8, 11, 15])

#         Set some partial coordinates with `-1`.

#         >>> cube.set_coord('pcp', (1273, -1))      # same as cube.set_coord('pcp', 1273)
#         >>> cube.permutation[:8]
#         array([-1, -1,  0, -1,  1,  3, -1,  2])
#         >>> cube.set_coord('pep', (2633, -1, -1))  # same as cube.set_coord('pep', 2633)
#         >>> cube.permutation[8:]
#         array([-1, 18, -1, 19, -1, -1, -1, 17, 16, -1, -1, -1])
#         """
#         if not isinstance(coord_type, str):
#             raise TypeError(f"coord_type must be str, not {type(coord_type).__name__}")
#         if not isinstance(coord, (int, tuple)):
#             raise TypeError(f"coord must be int or tuple, not {type(coord).__name__}")

#         if coord_type in ("co", "eo"):
#             if isinstance(coord, int):
#                 v = 3 if coord_type == "co" else 2
#                 orientation = self.orientation[:NUM_CORNERS] if coord_type == "co" else self.orientation[NUM_CORNERS:]
#                 orientation[:-1] = utils.get_orientation_array(coord, v, len(orientation) - 1)
#                 orientation[-1] = -np.sum(orientation[:-1]) % v
#             else:
#                 raise ValueError(f"coord must be int for coord_type '{coord_type}', not {type(coord).__name__}")

#         elif coord_type in ("cp", "ep", "pcp", "pep"):
#             permutation = self.permutation[:NUM_CORNERS] if coord_type in ("cp", "pcp") else self.permutation[NUM_CORNERS:]
#             if coord_type in ("cp", "ep"):
#                 if isinstance(coord, int):
#                     if coord_type == "cp":
#                         permutation[:], corner_parity = utils.get_permutation_array(coord, NUM_CORNERS)
#                         edge_parity = utils.get_permutation_parity(self.permutation[NUM_CORNERS:])
#                     elif coord_type == "ep":
#                         corner_parity = utils.get_permutation_parity(self.permutation[:NUM_CORNERS])
#                         permutation[:], edge_parity = utils.get_permutation_array(coord, NUM_EDGES, even_parity=True)
#                         permutation += NUM_CORNERS
#                     if corner_parity != edge_parity:
#                         self.permutation[-2:] = self.permutation[[-1, -2]]
#                 else:
#                     raise ValueError(f"coord must be int for coord_type '{coord_type}', not {type(coord).__name__}")
#             else:  # TODO test this and print representation
#                 if isinstance(coord, int):
#                     coord = (coord, EMPTY, EMPTY)
#                 permutation[:] = EMPTY
#                 axis_offset = CORNER_AXIS_OFFSET if coord_type == "pcp" else EDGE_AXIS_OFFSET
#                 for axis in range(2 if coord_type == "pcp" else 3):
#                     if coord[axis] != EMPTY:
#                         partial_permutation, combination = utils.get_partial_permutation_array(coord[axis], NUM_AXIS_ELEMS)
#                         permutation[combination] = partial_permutation + axis_offset[axis]

#         else:
#             raise ValueError(f"coord_type must be one of 'co', 'eo', 'cp', 'ep', 'pcp', or 'pep' (got '{coord_type}')")

#         self._coords = ()

#     def set_coords(self, coords: tuple, partial_corner_perm: bool = False, partial_edge_perm: bool = False):
#         """
#         Set cube coordinates.

#         Set the corner orientation, edge orientation,
#         (partial) corner permutation and (partial) edge permutation coordinates.

#         Parameters
#         ----------
#         coords : tuple
#             Cube coordinates in the following order:
#             corner orientation, edge orientation, (partial) corner permutation, (partial) edge permutation.
#         partial_corner_perm : bool
#             If `True`, use the partial corner permutation coordinate, otherwise use the normal corner permutation coordinate.
#         partial_edge_perm : bool
#             If `True`, use the partial edge permutation coordinate, otherwise use the normal edge permutation coordinate.

#         See Also
#         --------
#         utils.get_orientation_array
#         utils.get_permutation_array
#         utils.get_combination_array
#         utils.get_partial_permutation_array

#         Examples
#         --------
#         >>> from cube_solver import Cube
#         >>> cube = Cube()

#         Set cube coordinates.

#         >>> coords = (456, 673, 28179, 193381554)
#         >>> cube.set_coords(coords)
#         >>> cube.orientation
#         array([0, 1, 2, 1, 2, 2, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0])
#         >>> cube.permutation
#         array([ 5,  4,  0,  7,  1,  3,  6,  2, 12, 18, 10, 19,  9, 13, 14, 17, 16,
#                 8, 11, 15])

#         Set cube coordinates with partial corner permutation and partial edge permutation.

#         >>> cube.partial_corner_perm = True
#         >>> cube.partial_edge_perm = True
#         >>> coords = (456, 673, (1273, 391), (2633, 8640, 7262))
#         >>> cube.set_coords(coords)
#         >>> cube.orientation
#         array([0, 1, 2, 1, 2, 2, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0])
#         >>> cube.permutation
#         array([ 5,  4,  0,  7,  1,  3,  6,  2, 12, 18, 10, 19,  9, 13, 14, 17, 16,
#                 8, 11, 15])
#         """
#         # print("setter")
#         if not isinstance(coords, tuple):
#             raise TypeError(f"coords must be tuple, not {type(coords).__name__}")
#         if not isinstance(partial_corner_perm, bool):
#             raise TypeError(f"partial_corner_perm must be bool, not {type(partial_corner_perm).__name__}")
#         if not isinstance(partial_edge_perm, bool):
#             raise TypeError(f"partial_edge_perm must be bool, not {type(partial_edge_perm).__name__}")

#         self.set_coord("co", coords[0])
#         self.set_coord("eo", coords[1])
#         self.set_coord("pcp" if partial_corner_perm else "cp", coords[2])
#         self.set_coord("pep" if partial_edge_perm else "ep", coords[3])

#     def set_random_state(self):
#         """
#         Set a uniform random state.

#         Sets random corner orientation, edge orientation,
#         corner permutation, and edge permutation coordinates.
#         """
#         self.set_coord("co", np.random.randint(CORNER_ORIENTATION_SIZE))
#         self.set_coord("eo", np.random.randint(EDGE_ORIENTATION_SIZE))
#         self.set_coord("cp", np.random.randint(CORNER_PERMUTATION_SIZE))
#         self.set_coord("ep", np.random.randint(EDGE_PERMUTATION_SIZE))


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


# if __name__ == "__main__":
#     cube = Cube("U5 F R", random_state=True)
