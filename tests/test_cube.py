import pytest
import numpy as np
from typing import Union, List

from cube_solver import Cube, Maneuver, apply_move, apply_maneuver
from cube_solver.cube import utils
from cube_solver.cube.enums import Axis, Orbit, Layer, Color, Face, Cubie, Move


def test_enums():
    # axis
    assert hasattr(Axis, "NONE")
    assert all(axis.is_cartesian for axis in Axis.cartesian_axes())
    assert all(axis.is_diagonal for axis in Axis.diagonal_axes())
    assert all(axis.is_edge for axis in Axis.edge_axes())
    assert len([*Axis.axes()]) == 13
    assert len([*Axis.cartesian_axes()]) == 3
    assert len([*Axis.diagonal_axes()]) == 4
    assert len([*Axis.edge_axes()]) == 6

    # orbit
    assert hasattr(Orbit, "NONE")
    assert all(orbit.is_slice for orbit in Orbit.slices())
    assert all(orbit.is_tetrad for orbit in Orbit.tetrads())
    assert len([*Orbit.orbits()]) == 5
    assert len([*Orbit.slices()]) == 3
    assert len([*Orbit.tetrads()]) == 2

    # layer
    assert hasattr(Layer, "NONE")
    assert Layer.NONE.char == "N"
    assert Layer.NONE.axis == Axis.NONE
    assert Layer.NONE.perm == [[Cubie.NONE]]
    assert all(layer.is_outer for layer in Layer.outers())
    assert all(layer.is_inner for layer in Layer.inners())
    with pytest.raises(TypeError, match=r"char must be str, not NoneType"):
        Layer.from_char(None)
    with pytest.raises(ValueError, match=r"invalid face character \(got 'None'\)"):
        Layer.from_char("None")
    assert all(layer == Layer.from_char(layer.char) for layer in Layer)
    assert len([*Layer.layers()]) == 9
    assert len([*Layer.outers()]) == 6
    assert len([*Layer.inners()]) == 3

    # color
    assert hasattr(Color, "NONE")
    assert Color.NONE.char == 'N'
    with pytest.raises(TypeError, match=r"char must be str, not NoneType"):
        Color.from_char(None)
    with pytest.raises(ValueError, match=r"invalid color character \(got 'None'\)"):
        Color.from_char("None")
    assert all(color == Color.from_char(color.char) for color in Color)
    assert len([*Color.colors()]) == 6

    # face
    assert hasattr(Face, "NONE")
    assert Face.NONE.char == 'N'
    assert Face.NONE.axis == Axis.NONE
    assert Face.UP.axis == Axis.Y
    assert Face.NONE.opposite == Face.NONE
    assert all([face.axis == face.opposite.axis for face in Face])
    assert Face.NONE._index == (slice(None), slice(None), slice(None))
    with pytest.raises(TypeError, match=r"char must be str, not NoneType"):
        Face.from_char(None)
    with pytest.raises(ValueError, match=r"invalid face character \(got 'None'\)"):
        Face.from_char("None")
    assert all(face == Face.from_char(face.char) for face in Face)
    assert len([*Face.faces()]) == 6

    # cubie
    assert hasattr(Cubie, "NONE")
    assert Cubie.NONE.axis == Axis.NONE
    assert Cubie.U.axis == Axis.Y
    assert Cubie.NONE.orbit == Orbit.NONE
    assert Cubie.U.orbit == Orbit.NONE
    assert Cubie.NONE.faces == [Face.NONE]
    assert Cubie.CORE.faces == []
    assert Cubie.UBL.faces == [Face.UP, Face.LEFT, Face.BACK]
    assert Cubie.UBR.faces == [Face.UP, Face.BACK, Face.RIGHT]
    assert Cubie.NONE._index == (1, 1, 1)
    assert all(cubie.is_corner for cubie in Cubie.corners())
    assert all(cubie.is_edge for cubie in Cubie.edges())
    assert all(cubie.is_center for cubie in Cubie.centers())
    with pytest.raises(TypeError, match=r"faces must be list, not NoneType"):
        Cubie.from_faces(None)
    with pytest.raises(TypeError, match=r"faces elements must be Face, not NoneType"):
        Cubie.from_faces([None])
    with pytest.raises(ValueError, match=r"invalid cubie faces \(got \[Face.NONE, Face.NONE\]\)"):
        Cubie.from_faces([Face.NONE, Face.NONE])
    with pytest.raises(ValueError, match=r"faces length must be at most 3 \(got 4\)"):
        Cubie.from_faces([Face.NONE, Face.NONE, Face.NONE, Face.NONE])
    assert all(cubie == Cubie.from_faces(cubie.faces) for cubie in Cubie)
    assert len([*Cubie.cubies()]) == 27
    assert len([*Cubie.corners()]) == 8
    assert len([*Cubie.edges()]) == 12
    assert len([*Cubie.centers()]) == 6

    # move
    assert hasattr(Move, "NONE")
    assert Move.NONE.string == ""
    assert Move.U1.string == "U"
    assert Move.U2.string == "U2"
    assert Move.U3.string == "U'"
    assert Move.X1.string == "x"
    assert Move.X2.string == "x2"
    assert Move.X3.string == "x'"
    assert Move.UW1.string == "Uw"
    assert Move.UW2.string == "Uw2"
    assert Move.UW3.string == "Uw'"
    assert Move.NONE.axis == Axis.NONE
    assert Move.X1.axis == Axis.X
    assert Move.U1.axis == Axis.Y
    assert Move.NONE.inverse == Move.NONE
    assert Move.U1.inverse == Move.U3
    assert Move.NONE.layers == [Layer.NONE]
    assert Move.X1.layers == [Layer.RIGHT, Layer.LEFT, Layer.MIDDLE]
    assert Move.FW1.layers == [Layer.FRONT, Layer.STANDING]
    assert Move.U1.layers == [Layer.UP]
    assert Move.NONE.shifts == [0]
    assert Move.X1.shifts == [1, -1, -1]
    assert Move.FW1.shifts == [1, 1]
    assert Move.RW1.shifts == [1, -1]
    assert Move.U1.shifts == [1]
    assert all(move.is_face for move in Move.face_moves())
    assert all(move.is_slice for move in Move.slice_moves())
    assert all(move.is_wide for move in Move.wide_moves())
    assert all(move.is_rotation for move in Move.rotations())
    with pytest.raises(TypeError, match=r"string must be str, not NoneType"):
        Move.from_string(None)
    with pytest.raises(ValueError, match=r"invalid move string \(got 'None'\)"):
        Move.from_string("None")
    assert all(move == Move.from_string(move.string) for move in Move)
    assert all(move == Move.from_string(move.string[0].lower() + move.string[2:]) for move in Move.wide_moves())
    assert len([*Move.moves()]) == 54
    assert len([*Move.face_moves()]) == 18
    assert len([*Move.slice_moves()]) == 9
    assert len([*Move.wide_moves()]) == 18
    assert len([*Move.rotations()]) == 9


def test_utils():
    # orientation array
    with pytest.raises(TypeError, match=r"coord must be int, not NoneType"):
        utils.get_orientation_array(None, None, None, None)
    with pytest.raises(TypeError, match=r"v must be int, not NoneType"):
        utils.get_orientation_array(-1, None, None, None)
    with pytest.raises(TypeError, match=r"n must be int, not NoneType"):
        utils.get_orientation_array(-1, 0, None, None)
    with pytest.raises(TypeError, match=r"force_modulo must be bool, not NoneType"):
        utils.get_orientation_array(-1, 0, 0, None)
    with pytest.raises(ValueError, match=r"v must be positive \(got 0\)"):
        utils.get_orientation_array(-1, 0, 0)
    with pytest.raises(ValueError, match=r"n must be positive \(got 0\)"):
        utils.get_orientation_array(-1, 1, 0)
    # foce_modulo = False
    with pytest.raises(ValueError, match=r"coord must be >= 0 and < 1 \(got -1\)"):
        utils.get_orientation_array(-1, 1, 1)
    assert np.all(utils.get_orientation_array(0, 1, 1) == [0])
    with pytest.raises(ValueError, match=r"coord must be >= 0 and < 65536 \(got 65536\)"):
        utils.get_orientation_array(65536, 4, 8)
    assert np.all(utils.get_orientation_array(58596, 4, 8) == [3, 2, 1, 0, 3, 2, 1, 0])
    # force_module = True
    with pytest.raises(ValueError, match=r"coord must be >= 0 and < 1 \(got -1\)"):
        utils.get_orientation_array(-1, 1, 1, True)
    assert np.all(utils.get_orientation_array(0, 1, 1, True) == [0])
    with pytest.raises(ValueError, match=r"coord must be >= 0 and < 16384 \(got 16384\)"):
        utils.get_orientation_array(16384, 4, 8, True)
    assert np.all(utils.get_orientation_array(14649, 4, 8, True) == [3, 2, 1, 0, 3, 2, 1, 0])

    # orientation coord
    with pytest.raises(TypeError, match=r"orientation must be ndarray, not NoneType"):
        utils.get_orientation_coord(None, None, None)
    with pytest.raises(TypeError, match=r"v must be int, not NoneType"):
        utils.get_orientation_coord(np.array([]), None, None)
    with pytest.raises(TypeError, match=r"is_modulo must be bool, not NoneType"):
        utils.get_orientation_coord(np.array([]), 0, None)
    with pytest.raises(ValueError, match=r"v must be positive \(got 0\)"):
        utils.get_orientation_coord(np.array([]), 0)
    with pytest.raises(ValueError, match=r"orientation length must be positive \(got 0\)"):
        utils.get_orientation_coord(np.array([]), 1)
    with pytest.raises(TypeError, match=r"orientation elements must be int, not object"):
        utils.get_orientation_coord(np.array([None]), 1)
    # is_modulo = False
    with pytest.raises(ValueError, match=r"orientation values must be >= 0 and < 1 \(got \[-1\]\)"):
        utils.get_orientation_coord(np.array([-1]), 1)
    assert utils.get_orientation_coord(np.array([0]), 1) == 0
    with pytest.raises(ValueError, match=r"orientation values must be >= 0 and < 4 \(got \[4\]\)"):
        utils.get_orientation_coord(np.array([4]), 4)
    coord, v, n = 58596, 4, 8
    array = utils.get_orientation_array(coord, v, n)
    assert utils.get_orientation_coord(array, v) == coord
    assert np.all(utils.get_orientation_array(utils.get_orientation_coord(array, v), v, n) == array)
    # is_modulo = True
    with pytest.raises(ValueError, match=r"orientation values must be >= 0 and < 1 \(got \[-1\]\)"):
        utils.get_orientation_coord(np.array([-1]), 1, True)
    assert utils.get_orientation_coord(np.array([0]), 1, True) == 0
    with pytest.raises(ValueError, match=r"orientation values must be >= 0 and < 4 \(got \[4\]\)"):
        utils.get_orientation_coord(np.array([4]), 4, True)
    with pytest.raises(ValueError, match=r"orientation has no total orientation 0 modulo 4 \(got 1\)"):
        utils.get_orientation_coord(np.array([1]), 4, True)
    coord, v, n = 14649, 4, 8
    array = utils.get_orientation_array(coord, v, n, True)
    assert utils.get_orientation_coord(array, v, True) == coord
    assert np.all(utils.get_orientation_array(utils.get_orientation_coord(array, v, True), v, n, True) == array)

    # permutation array
    with pytest.raises(TypeError, match=r"coord must be int, not NoneType"):
        utils.get_permutation_array(None, None, None)
    with pytest.raises(TypeError, match=r"n must be int, not NoneType"):
        utils.get_permutation_array(-1, None, None)
    with pytest.raises(TypeError, match=r"force_even_parity must be bool, not NoneType"):
        utils.get_permutation_array(-1, 0, None)
    # force_even_parity = False
    with pytest.raises(ValueError, match=r"n must be positive \(got 0\)"):
        utils.get_permutation_array(-1, 0)
    with pytest.raises(ValueError, match=r"coord must be >= 0 and < 1 \(got -1\)"):
        utils.get_permutation_array(-1, 1)
    permutation, parity = utils.get_permutation_array(0, 1)
    assert np.all(permutation == [0])
    assert parity is False
    with pytest.raises(ValueError, match=r"coord must be >= 0 and < 40320 \(got 40320\)"):
        utils.get_permutation_array(40320, 8)
    permutation, parity = utils.get_permutation_array(16702, 8)
    assert np.all(permutation == [3, 2, 1, 0, 7, 6, 4, 5])
    assert parity is True
    # force_even_parity = True
    with pytest.raises(ValueError, match=r"n must be > 1 \(got 1\)"):
        utils.get_permutation_array(-1, 1, True)
    with pytest.raises(ValueError, match=r"coord must be >= 0 and < 1 \(got -1\)"):
        utils.get_permutation_array(-1, 2, True)
    permutation, parity = utils.get_permutation_array(0, 2, True)
    assert np.all(permutation == [0, 1])
    assert parity is False
    with pytest.raises(ValueError, match=r"coord must be >= 0 and < 20160 \(got 20160\)"):
        utils.get_permutation_array(20160, 8, True)
    permutation, parity = utils.get_permutation_array(8351, 8, True)
    assert np.all(permutation == [3, 2, 1, 0, 7, 6, 5, 4])
    assert parity is False
    # no pre-computed factorial
    permutation, parity = utils.get_permutation_array(1000000001, 13)
    assert np.all(permutation == [2, 1, 0, 8, 10, 7, 9, 12, 4, 6, 11, 5, 3])
    assert parity is True

    # permutation coord
    with pytest.raises(TypeError, match=r"permutation must be ndarray, not NoneType"):
        utils.get_permutation_coord(None, None)
    with pytest.raises(TypeError, match=r"is_even_parity must be bool, not NoneType"):
        utils.get_permutation_coord(np.array([]), None)
    with pytest.raises(ValueError, match=r"permutation length must be positive \(got 0\)"):
        utils.get_permutation_coord(np.array([]))
    with pytest.raises(TypeError, match=r"permutation elements must be int, not object"):
        utils.get_permutation_coord(np.array([None]))
    # is_even_parity = False
    assert utils.get_permutation_coord(np.array([0])) == 0
    with pytest.raises(ValueError, match=r"permutation values must be different \(got \[0 0\]\)"):
        utils.get_permutation_coord(np.array([0, 0]))
    assert utils.get_permutation_coord(np.array([0, -1])) == 1
    coord, n = 16702, 8
    array = utils.get_permutation_array(coord, n)[0]
    assert utils.get_permutation_coord(array) == coord
    assert np.all(utils.get_permutation_array(utils.get_permutation_coord(array), n)[0] == array)
    # is_even_parity = True
    with pytest.raises(ValueError, match=r"permutation length must be > 1 \(got 1\)"):
        utils.get_permutation_coord(np.array([0]), True)
    with pytest.raises(ValueError, match=r"permutation values must be different \(got \[0 0\]\)"):
        utils.get_permutation_coord(np.array([0, 0]), True)
    with pytest.raises(ValueError, match=r"permutation has no even parity \(got \[1 0\]\)"):
        utils.get_permutation_coord(np.array([1, 0]), True)
    assert utils.get_permutation_coord(np.array([0, 1]), True) == 0
    assert utils.get_permutation_coord(np.array([0, 1, -1]), True) == 1
    coord, n = 8351, 8
    array = utils.get_permutation_array(coord, n, True)[0]
    assert utils.get_permutation_coord(array, True) == coord
    assert np.all(utils.get_permutation_array(utils.get_permutation_coord(array, True), n, True)[0] == array)
    # no pre-computed factorial
    coord, n = 1000000001, 13
    array = utils.get_permutation_array(coord, n)[0]
    assert utils.get_permutation_coord(array) == coord
    assert np.all(utils.get_permutation_array(utils.get_permutation_coord(array), n)[0] == array)

    # permutation parity
    with pytest.raises(TypeError, match=r"permutation must be ndarray, not NoneType"):
        utils.get_permutation_parity(None)
    with pytest.raises(ValueError, match=r"permutation length must be positive \(got 0\)"):
        utils.get_permutation_parity(np.array([]))
    with pytest.raises(TypeError, match=r"permutation elements must be int, not object"):
        utils.get_permutation_parity(np.array([None]))
    assert utils.get_permutation_parity(np.array([0])) is False
    with pytest.raises(ValueError, match=r"permutation values must be different \(got \[0 0\]\)"):
        utils.get_permutation_parity(np.array([0, 0]))
    assert utils.get_permutation_parity(np.array([0, 1])) is False
    assert utils.get_permutation_parity(np.array([-1, 0])) is False
    assert utils.get_permutation_parity(np.array([1, 0])) is True
    assert utils.get_permutation_parity(np.array([0, -1])) is True

    # combination array
    with pytest.raises(TypeError, match=r"coord must be int, not NoneType"):
        utils.get_combination_array(None, None)
    with pytest.raises(TypeError, match=r"n must be int, not NoneType"):
        utils.get_combination_array(-1, None)
    with pytest.raises(ValueError, match=r"n must be positive \(got 0\)"):
        utils.get_combination_array(-1, 0)
    with pytest.raises(ValueError, match=r"coord must be >= 0 \(got -1\)"):
        utils.get_combination_array(-1, 1)
    assert np.all(utils.get_combination_array(0, 1) == [0])
    assert np.all(utils.get_combination_array(450, 4) == [0, 1, 10, 11])
    # no pre-computed combination
    assert np.all(utils.get_combination_array(286, 3) == [0, 1, 13])
    assert np.all(utils.get_combination_array(450, 5) == [1, 6, 8, 9, 10])

    # combination coordinate
    with pytest.raises(TypeError, match=r"combination must be ndarray, not NoneType"):
        utils.get_combination_coord(None)
    with pytest.raises(ValueError, match=r"combination length must be positive \(got 0\)"):
        utils.get_combination_coord(np.array([]))
    with pytest.raises(TypeError, match=r"combination elements must be int, not object"):
        utils.get_combination_coord(np.array([None]))
    with pytest.raises(ValueError, match=r"combination values must be >= 0 \(got \[-1\]\)"):
        utils.get_combination_coord(np.array([-1]))
    assert utils.get_combination_coord(np.array([0])) == 0
    with pytest.raises(ValueError, match=r"combination values must be in increasing order \(got \[0 0\]\)"):
        utils.get_combination_coord(np.array([0, 0]))
    assert utils.get_combination_coord(np.array([0, 1])) == 0
    coord, n = 450, 4
    array = utils.get_combination_array(coord, n)
    assert utils.get_combination_coord(array) == coord
    assert np.all(utils.get_combination_array(utils.get_combination_coord(array), n) == array)
    # no pre-computed combination
    coord, n = 286, 3
    array = utils.get_combination_array(coord, n)
    assert utils.get_combination_coord(array) == coord
    assert np.all(utils.get_combination_array(utils.get_combination_coord(array), n) == array)
    coord, n = 450, 5
    array = utils.get_combination_array(coord, n)
    assert utils.get_combination_coord(array) == coord
    assert np.all(utils.get_combination_array(utils.get_combination_coord(array), n) == array)

    # partial permutation array
    with pytest.raises(TypeError, match=r"coord must be int, not NoneType"):
        utils.get_partial_permutation_array(None, None)
    with pytest.raises(TypeError, match=r"n must be int, not NoneType"):
        utils.get_partial_permutation_array(-1, None)
    with pytest.raises(ValueError, match=r"n must be positive \(got 0\)"):
        utils.get_partial_permutation_array(-1, 0)
    with pytest.raises(ValueError, match=r"coord must be >= 0 \(got -1\)"):
        utils.get_partial_permutation_array(-1, 1)
    permutation, combination = utils.get_partial_permutation_array(0, 1)
    assert np.all(permutation == [0])
    assert np.all(combination == [0])
    permutation, combination = utils.get_partial_permutation_array(450, 4)
    assert np.all(permutation == [3, 0, 1, 2])
    assert np.all(combination == [1, 2, 3, 6])
    # no pre-computed factorial
    permutation, combination = utils.get_partial_permutation_array(100000000000001, 13)
    assert np.all(permutation == [0, 7, 11, 3, 4, 2, 5, 12, 6, 9, 10, 8, 1])
    assert np.all(combination == [1, 2, 3, 5, 8, 9, 10, 11, 12, 13, 14, 17, 18])

    # partial permutation coord
    with pytest.raises(TypeError, match=r"permutation must be ndarray, not NoneType"):
        utils.get_partial_permutation_coord(None, None)
    with pytest.raises(TypeError, match=r"combination must be ndarray, not NoneType"):
        utils.get_partial_permutation_coord(np.array([]), None)
    with pytest.raises(ValueError, match=r"permutation length must be positive \(got 0\)"):
        utils.get_partial_permutation_coord(np.array([]), np.array([]))
    with pytest.raises(ValueError, match=r"permutation length and combination length must be the same \(got 1 != 0\)"):
        utils.get_partial_permutation_coord(np.array([None]), np.array([]))
    with pytest.raises(TypeError, match=r"permutation elements must be int, not object"):
        utils.get_partial_permutation_coord(np.array([None]), np.array([None]))
    with pytest.raises(TypeError, match=r"combination elements must be int, not object"):
        utils.get_partial_permutation_coord(np.array([0]), np.array([None]))
    with pytest.raises(ValueError, match=r"combination values must be >= 0 \(got \[-1\]\)"):
        utils.get_partial_permutation_coord(np.array([0]), np.array([-1]))
    assert utils.get_partial_permutation_coord(np.array([0]), np.array([0])) == 0
    with pytest.raises(ValueError, match=r"permutation values must be different \(got \[0 0\]\)"):
        utils.get_partial_permutation_coord(np.array([0, 0]), np.array([0, 0]))
    with pytest.raises(ValueError, match=r"combination values must be in increasing order \(got \[0 0\]\)"):
        utils.get_partial_permutation_coord(np.array([0, 1]), np.array([0, 0]))
    assert utils.get_partial_permutation_coord(np.array([0, 1]), np.array([0, 1])) == 0
    coord, n = 450, 4
    permutation, combination = utils.get_partial_permutation_array(coord, n)
    assert utils.get_partial_permutation_coord(permutation, combination) == coord
    perm, comb = utils.get_partial_permutation_array(utils.get_partial_permutation_coord(permutation, combination), n)
    assert np.all(perm == permutation)
    assert np.all(comb == combination)
    # no pre-computed factorial
    coord, n = 100000000000001, 13
    permutation, combination = utils.get_partial_permutation_array(coord, n)
    assert utils.get_partial_permutation_coord(permutation, combination) == coord
    perm, comb = utils.get_partial_permutation_array(utils.get_partial_permutation_coord(permutation, combination), n)
    assert np.all(perm == permutation)
    assert np.all(comb == combination)


def check_cube(cube: Cube, permutation_parity: Union[bool,  None], orientation: List[int], permutation: List[int]):
    # state
    assert np.all(cube.orientation == orientation)
    assert np.all(cube.permutation == permutation)
    assert cube.permutation_parity == permutation_parity


def test_cube():
    cube = Cube()
    assert cube.is_solved
    assert repr(cube) == "WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY"
    check_cube(cube, False,
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])

    # cube comparison
    assert cube != 0
    assert cube == Cube()

    # apply_maneuver
    with pytest.raises(TypeError, match=r"maneuver must be str, not NoneType"):
        cube.apply_maneuver(None)
    with pytest.raises(TypeError, match=r"cube must be Cube, not NoneType"):
        apply_maneuver(None, "")
    next_cube = apply_maneuver(cube, "U F2 R'")
    assert cube.is_solved
    assert not next_cube.is_solved
    assert repr(cube) == "WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY"
    assert repr(next_cube) == "WWBWWBYYOGGROOROOBGGWGGWRRYBRRBRROOGYOOYBBWBBWWGYYGYYR"
    check_cube(next_cube, True,
               [0, 2, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [5, 0, 1, 4, 2, 7, 6, 3, 12, 11, 10, 13, 9, 17, 14, 18, 16, 15, 19, 8])

    # scramble
    with pytest.raises(TypeError, match=r"scramble must be str or None, not int"):
        Cube(0)
    cube = Cube(Maneuver([Move.U1, Move.F2, Move.R3]))
    assert not cube.is_solved
    assert repr(cube) == "WWBWWBYYOGGROOROOBGGWGGWRRYBRRBRROOGYOOYBBWBBWWGYYGYYR"
    check_cube(cube, True,
               [0, 2, 2, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [5, 0, 1, 4, 2, 7, 6, 3, 12, 11, 10, 13, 9, 17, 14, 18, 16, 15, 19, 8])
    cube = Cube("U F2 R' D B2 L' M E2 S' Uw Fw2 Rw' Dw Bw2 Lw' u f2 r' d b2 l' x y2 z'")
    assert not cube.is_solved
    assert repr(cube) == "YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOYO"
    check_cube(cube, False,
               [0, 0, 2, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [4, 3, 1, 6, 7, 0, 2, 5, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])

    # repr
    with pytest.raises(TypeError, match=r"repr must be str or None, not int"):
        Cube(repr=0)
    with pytest.raises(ValueError, match=r"repr length must be 54 \(got 4\)"):
        Cube(repr="None")
    match = r"invalid string representation, setting undefined orientation and permutation values with -1"
    with pytest.warns(UserWarning, match=match):
        # original            Y
        cube = Cube(repr="YGWYNOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOYO")
        assert repr(cube) == "NGWNNNNWWNNNGRRWGBNBGBGBRRNRNOOONNONGRNNBWNNBWWNOWWONN"
        cube = Cube(repr=repr(cube))
    check_cube(cube, None,
               [-1, 0, -1, 0, 0, -1, 2, -1, 0, 0, -1, 0, -1, -1, 1, 0, 0, -1, 0, 0],
               [-1, 3, -1, 6, 7, -1, 2, -1, 18, 10, -1, 14, -1, -1, 19, 15, 11, -1, 16, 17])
    with pytest.warns(UserWarning, match=r"invalid corner orientation"):
        # original        Y        B                            O
        cube = Cube(repr="OGWYYOBWWYGRGRRWGBYBGBGBRRYRYOOOBGOGGRBYBWYRBWWROWWOYO")
        assert repr(cube) == "OGWYYOBWWYGRGRRWGBYBGBGBRRYRYOOOBGOGGRBYBWYRBWWROWWOYO"
        cube = Cube(repr=repr(cube))
    check_cube(cube, False,
               [1, 0, 2, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [4, 3, 1, 6, 7, 0, 2, 5, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])
    with pytest.warns(UserWarning, match=r"invalid edge orientation"):
        # original         G                                   R
        cube = Cube(repr="YRWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGGOYBWYRBWWROWWOYO")
        assert repr(cube) == "YRWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGGOYBWYRBWWROWWOYO"
        cube = Cube(repr=repr(cube))
    check_cube(cube, False,
               [0, 0, 2, 0, 0, 1, 2, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [4, 3, 1, 6, 7, 0, 2, 5, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])
    with pytest.warns(UserWarning, match=r"invalid corner permutation"):
        # original        Y        B                            O
        cube = Cube(repr="WGWYYOBWWOGRGRRWGBYBGBGBRRYRYOOOBGOGGRBYBWYRBWWROWWOYO")
        assert repr(cube) == "WGWYYOBWWOGRGRRWGBYBGBGBRRYRYOOOBGOGGRBYBWYRBWWROWWOYO"
        cube = Cube(repr=repr(cube))
    check_cube(cube, None,
               [0, 0, 2, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [2, 3, 1, 6, 7, 0, 2, 5, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])
    with pytest.warns(UserWarning, match=r"invalid edge permutation"):
        # original         G
        cube = Cube(repr="YBWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOYO")
        assert repr(cube) == "YBWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOYO"
        cube = Cube(repr=repr(cube))
    check_cube(cube, None,
               [0, 0, 2, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [4, 3, 1, 6, 7, 0, 2, 5, 16, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])
    with pytest.warns(UserWarning, match=r"invalid cube parity"):
        # original         G     W           B                 R
        cube = Cube(repr="YWWYYOBGWBGRGRRWGBYRGBGBRRYRYOOOBGOGGBOYBWYRBWWROWWOYO")
        assert repr(cube) == "YWWYYOBGWBGRGRRWGBYRGBGBRRYRYOOOBGOGGBOYBWYRBWWROWWOYO"
        cube = Cube(repr=repr(cube))
    check_cube(cube, None,
               [0, 0, 2, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [4, 3, 1, 6, 7, 0, 2, 5, 10, 18, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])
    cube = Cube(repr="YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOYO")
    assert repr(cube) == "YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOYO"
    cube = Cube(repr=repr(cube))
    check_cube(cube, False,
               [0, 0, 2, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [4, 3, 1, 6, 7, 0, 2, 5, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])

    # string
    assert str(cube) == """        ---------
        | Y G W |
        | Y Y O |
        | B W W |
---------------------------------
| B G R | Y B G | R Y O | G R O |
| G R R | B G B | O O B | Y B W |
| W G B | R R Y | G O G | Y R B |
---------------------------------
        | W W R |
        | O W W |
        | O Y O |
        ---------"""

    # get_coord
    with pytest.raises(TypeError, match=r"coord_name must be str, not NoneType"):
        cube.get_coord(None)
    with pytest.raises(ValueError, match=r"coord_name must be one of 'co', 'eo', 'cp', 'ep', 'pcp', 'pep' \(got 'None'\)"):
        cube.get_coord("None")
    assert cube.get_coord("co") == 167
    assert cube.get_coord("eo") == 48
    assert cube.get_coord("cp") == 22530
    assert cube.get_coord("ep") == 203841327
    assert cube.get_coord("pcp") == (668, 1011)
    assert cube.get_coord("pep") == (4551, 11176, 1202)

    # get_coords
    with pytest.raises(TypeError, match=r"partial_corner_perm must be bool, not NoneType"):
        cube.get_coords(None, None)
    with pytest.raises(TypeError, match=r"partial_edge_perm must be bool, not NoneType"):
        cube.get_coords(False, None)
    assert cube.get_coords(False, False) == (167, 48, 22530, 203841327)
    assert cube.get_coords(True, True) == (167, 48, (668, 1011), (4551, 11176, 1202))

    # random_state
    with pytest.raises(TypeError, match=r"random_state must be bool, not int"):
        Cube(random_state=0)
    cube = Cube(random_state=True)
    assert cube.coords != (0, 0, 0, 0)
    assert repr(cube) != "WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY"

    # reset
    cube.reset()
    assert repr(cube) == "WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY"
    check_cube(cube, False,
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])

    # set_coord
    with pytest.raises(TypeError, match=r"coord_name must be str, not NoneType"):
        cube.set_coord(None, None)
    with pytest.raises(TypeError, match=r"coord must be int or tuple, not NoneType"):
        cube.set_coord("None", None)
    with pytest.raises(ValueError, match=r"coord_name must be one of 'co', 'eo', 'cp', 'ep', 'pcp', 'pep' \(got 'None'\)"):
        cube.set_coord("None", ())
    # corner orientation
    with pytest.raises(TypeError, match=r"coord must be int for coord_name 'co', not tuple"):
        cube.set_coord("co", ())
    with pytest.raises(ValueError, match=r"coord must be >= 0 and < 2187 \(got -1\)"):
        cube.set_coord("co", -1)
    cube.set_coord("co", 167)
    assert cube.get_coord("co") == 167
    check_cube(cube, False,
               [0, 0, 2, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    # edge orientation
    with pytest.raises(TypeError, match=r"coord must be int for coord_name 'eo', not tuple"):
        cube.set_coord("eo", ())
    with pytest.raises(ValueError, match=r"coord must be >= 0 and < 2048 \(got -1\)"):
        cube.set_coord("eo", -1)
    cube.set_coord("eo", 48)
    assert cube.get_coord("eo") == 48
    check_cube(cube, False,
               [0, 0, 2, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    # corner permutation
    with pytest.raises(TypeError, match=r"coord must be int for coord_name 'cp', not tuple"):
        cube.set_coord("cp", ())
    with pytest.raises(ValueError, match=r"coord must be >= 0 and < 40320 \(got -1\)"):
        cube.set_coord("cp", -1)
    cube.set_coord("cp", 22530)
    assert cube.get_coord("cp") == 22530
    check_cube(cube, False,
               [0, 0, 2, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [4, 3, 1, 6, 7, 0, 2, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    # edge permutation
    with pytest.raises(TypeError, match=r"coord must be int for coord_name 'ep', not tuple"):
        cube.set_coord("ep", ())
    with pytest.raises(ValueError, match=r"coord must be >= 0 and < 239500800 \(got -1\)"):
        cube.set_coord("ep", -1)
    cube.set_coord("ep", 203841327)
    assert cube.get_coord("ep") == 203841327
    check_cube(cube, False,
               [0, 0, 2, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [4, 3, 1, 6, 7, 0, 2, 5, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])
    # partial corner permutation
    with pytest.raises(ValueError, match=r"coord tuple length must be 2 for coord_name 'pcp' \(got 0\)"):
        cube.set_coord("pcp", ())
    with pytest.raises(TypeError, match=r"coord tuple elements must be int, not NoneType"):
        cube.set_coord("pcp", (None, None))
    with pytest.raises(ValueError, match=r"coord must be >= 0 and < 1680 \(got -2\)"):
        cube.set_coord("pcp", (-2, -2))
    with pytest.raises(ValueError, match=r"invalid partial coordinates, overlapping detected \(got \(0, 0\)\)"):
        cube.set_coord("pcp", (0, 0))
    cube.set_coord("pcp", (-1, 1011))
    assert cube.get_coord("pcp") == (-1, 1011)
    check_cube(cube, None,
               [0, -1, -1, 0, 0, -1, -1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [4, -1, -1, 6, 7, -1, -1, 5, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])
    cube.set_coord("pcp", (668, -1))
    assert cube.get_coord("pcp") == 668
    check_cube(cube, None,
               [-1, 0, 0, -1, -1, 0, 0, -1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [-1, 3, 1, -1, -1, 0, 2, -1, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])
    cube.set_coord("pcp", (-1, -1))
    assert cube.get_coord("pcp") == -1
    check_cube(cube, None,
               [-1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [-1, -1, -1, -1, -1, -1, -1, -1, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])
    cube.set_coord("pcp", 668)
    assert cube.get_coord("pcp") == 668
    check_cube(cube, None,
               [-1, 0, 0, -1, -1, 0, 0, -1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [-1, 3, 1, -1, -1, 0, 2, -1, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])
    cube.set_coord("pcp", (668, 1011))
    assert cube.get_coord("pcp") == (668, 1011)
    check_cube(cube, False,
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [4, 3, 1, 6, 7, 0, 2, 5, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])
    # partial edge permutation
    with pytest.raises(ValueError, match=r"coord tuple length must be 3 for coord_name 'pep' \(got 0\)"):
        cube.set_coord("pep", ())
    with pytest.raises(TypeError, match=r"coord tuple elements must be int, not NoneType"):
        cube.set_coord("pep", (None, None, None))
    with pytest.raises(ValueError, match=r"coord must be >= 0 and < 11880 \(got -2\)"):
        cube.set_coord("pep", (-2, -2, -2))
    with pytest.raises(ValueError, match=r"invalid partial coordinates, overlapping detected \(got \(0, 0, 0\)\)"):
        cube.set_coord("pep", (0, 0, 0))
    cube.set_coord("pep", (-1, -1, 1202))
    assert cube.get_coord("pep") == (-1, -1, 1202)
    check_cube(cube, None,
               [0, 0, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0, -1, 1, -1, 0, -1, -1, -1, -1],
               [4, 3, 1, 6, 7, 0, 2, 5, -1, -1, 12, 14, -1, 13, -1, 15, -1, -1, -1, -1])
    cube.set_coord("pep", (-1, 11176, -1))
    assert cube.get_coord("pep") == (-1, 11176, -1)
    check_cube(cube, None,
               [0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, 0, -1, -1, -1, 0, 0],
               [4, 3, 1, 6, 7, 0, 2, 5, 18, -1, -1, -1, -1, -1, 19, -1, -1, -1, 16, 17])
    cube.set_coord("pep", (4551, -1, -1))
    assert cube.get_coord("pep") == 4551
    check_cube(cube, None,
               [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, -1, 0, -1, -1, -1, 0, 0, -1, -1],
               [4, 3, 1, 6, 7, 0, 2, 5, -1, 10, -1, -1, 9, -1, -1, -1, 11, 8, -1, -1])
    cube.set_coord("pep", (-1, -1, -1))
    assert cube.get_coord("pep") == -1
    check_cube(cube, None,
               [0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
               [4, 3, 1, 6, 7, 0, 2, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1])
    cube.set_coord("pep", 4551)
    assert cube.get_coord("pep") == 4551
    check_cube(cube, None,
               [0, 0, 0, 0, 0, 0, 0, 0, -1, 0, -1, -1, 0, -1, -1, -1, 0, 0, -1, -1],
               [4, 3, 1, 6, 7, 0, 2, 5, -1, 10, -1, -1, 9, -1, -1, -1, 11, 8, -1, -1])
    cube.set_coord("pep", (4551, 11176, 1202))
    assert cube.get_coord("pep") == (4551, 11176, 1202)
    check_cube(cube, False,
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [4, 3, 1, 6, 7, 0, 2, 5, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])

    # permutation parity
    cube.set_coord("pcp", 0)
    assert cube.get_coord("pcp") == 0
    check_cube(cube, None,
               [0, 0, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 1, 2, 3, -1, -1, -1, -1, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])
    cube.set_coord("ep", 0)
    assert cube.get_coord("pcp") == 0
    assert cube.get_coord("ep") == 0
    check_cube(cube, None,
               [0, 0, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 1, 2, 3, -1, -1, -1, -1, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    cube.set_coord("pcp", (0, 1657))
    assert cube.get_coord("pcp") == (0, 1657)
    assert cube.get_coord("ep") == 0
    check_cube(cube, True,
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 1, 2, 3, 4, 5, 7, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 18])
    with pytest.warns(UserWarning, match=r"invalid cube parity"):
        cube.set_coord("pep", (0, 11856, 1656))
    assert cube.get_coord("pcp") == (0, 1657)
    assert cube.get_coord("pep") == (0, 11856, 1656)
    check_cube(cube, None,
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 1, 2, 3, 4, 5, 7, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    cube.set_coord("ep", 0)
    assert cube.get_coord("pcp") == (0, 1657)
    assert cube.get_coord("ep") == 0
    check_cube(cube, True,
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 1, 2, 3, 4, 5, 7, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 18])

    # set_coords
    with pytest.raises(TypeError, match=r"coords must be tuple, not NoneType"):
        cube.set_coords(None, None, None)
    with pytest.raises(TypeError, match=r"partial_corner_perm must be bool, not NoneType"):
        cube.set_coords((), None, None)
    with pytest.raises(TypeError, match=r"partial_edge_perm must be bool, not NoneType"):
        cube.set_coords((), False, None)
    with pytest.raises(ValueError, match=r"coords tuple length must be 4 \(got 0\)"):
        assert cube.set_coords((), False, False)
    cube.set_coords((167, 48, 22530, 203841327), False, False)
    assert repr(cube) == "WGYWWRBYYBGOGOOYGBWBGBGBOOWOWRRRBGRGGORWBYWOBYYORYYRWR"
    check_cube(cube, False,
               [0, 0, 2, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [4, 3, 1, 6, 7, 0, 2, 5, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])
    cube.set_coords((167, 48, (668, 1011), (4551, 11176, 1202)), True, True)
    assert repr(cube) == "WGYWWRBYYBGOGOOYGBWBGBGBOOWOWRRRBGRGGORWBYWOBYYORYYRWR"
    check_cube(cube, False,
               [0, 0, 2, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [4, 3, 1, 6, 7, 0, 2, 5, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])

    # coords
    cube.coords = (167, 48, 22530, 203841327)
    assert cube.coords == (167, 48, 22530, 203841327)
    assert repr(cube) == "WGYWWRBYYBGOGOOYGBWBGBGBOOWOWRRRBGRGGORWBYWOBYYORYYRWR"
    check_cube(cube, False,
               [0, 0, 2, 0, 0, 1, 2, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [4, 3, 1, 6, 7, 0, 2, 5, 18, 10, 12, 14, 9, 13, 19, 15, 11, 8, 16, 17])

    # apply_move
    with pytest.raises(TypeError, match=r"move must be Move, not NoneType"):
        cube.apply_move(None)
    cube.reset()
    cube.apply_move(Move.NONE)
    assert repr(cube) == "WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY"
    check_cube(cube, False,
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    cube.apply_move(Move.F1)
    assert repr(cube) == "WWWWWWOOOOOYOOYOOYGGGGGGGGGWRRWRRWRRBBBBBBBBBRRRYYYYYY"
    check_cube(cube, True,
               [0, 1, 0, 1, 0, 2, 0, 2, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1],
               [0, 5, 2, 7, 4, 3, 6, 1, 8, 18, 10, 19, 12, 13, 14, 15, 16, 17, 11, 9])
    cube.apply_move(Move.X1)
    assert repr(cube) == "GGGGGGGGGYYYOOOOOORRRYYYYYYWWWRRRRRROOOWWWWWWBBBBBBBBB"
    check_cube(cube, True,
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [5, 4, 2, 3, 0, 1, 6, 7, 12, 13, 10, 11, 9, 8, 14, 15, 16, 17, 18, 19])
    cube.set_coords((0, 0, 0, 0), True, True)
    cube.apply_move(Move.NONE)
    assert repr(cube) == "GGNNGNNGGONNNONNNONYYNYNYYNRNNNRNNNRNWWNWNWWNBBNNBNNBB"
    check_cube(cube, None,
               [0, 0, 0, 0, -1, -1, -1, -1, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1],
               [0, 1, 2, 3, -1, -1, -1, -1, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1, -1])
    cube.apply_move(Move.F1)
    assert repr(cube) == "GGNNGNONNONBNOBNNNYNNYYYNNYNNNGRNGNRNWWNWNWWNNNRNBNNBB"
    check_cube(cube, None,
               [0, -1, 0, -1, -1, 2, -1, 2, 0, -1, 0, -1, -1, -1, -1, -1, -1, -1, 1, 1],
               [0, -1, 2, -1, -1, 3, -1, 1, 8, -1, 10, -1, -1, -1, -1, -1, -1, -1, 11, 9])
    cube.apply_move(Move.X1)
    assert repr(cube) == "YNNYYYNNYBBNNONONNNNRNBNNBBGGNNRNRNNNNONGNNGGNWWNWNWWN"
    check_cube(cube, None,
               [0, 0, -1, -1, -1, -1, 0, 0, -1, -1, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1],
               [5, 4, -1, -1, -1, -1, 6, 7, -1, -1, 10, 11, 9, 8, -1, -1, -1, -1, -1, -1])
    with pytest.raises(TypeError, match=r"cube must be Cube, not NoneType"):
        apply_move(None, Move.NONE)
    cube.reset()
    next_cube = apply_move(cube, Move.NONE)
    assert repr(cube) == "WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY"
    assert repr(next_cube) == "WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY"
    next_cube = apply_move(cube, Move.F1)
    assert repr(cube) == "WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY"
    assert repr(next_cube) == "WWWWWWOOOOOYOOYOOYGGGGGGGGGWRRWRRWRRBBBBBBBBBRRRYYYYYY"
    next_cube = apply_move(cube, Move.X1)
    assert repr(cube) == "WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY"
    assert repr(next_cube) == "GGGGGGGGGOOOOOOOOOYYYYYYYYYRRRRRRRRRWWWWWWWWWBBBBBBBBB"


def test_maneuver():
    with pytest.raises(TypeError, match=r"moves must be str or list, not NoneType"):
        Maneuver(None, None)
    with pytest.raises(TypeError, match=r"reduce must be bool, not NoneType"):
        Maneuver([None], None)
    with pytest.raises(TypeError, match="moves list elements must be Move, not NoneType"):
        Maneuver([None], True)
    maneuver = Maneuver("")
    assert maneuver != 0
    assert maneuver != [None]
    assert maneuver != "None"
    assert maneuver == ""
    assert maneuver == []
    assert maneuver == Maneuver("")
    assert maneuver == Maneuver([])
    assert maneuver.moves == ()
    maneuver = Maneuver([Move.NONE])
    assert maneuver == ""
    assert maneuver == []
    assert maneuver == Maneuver("")
    assert maneuver == Maneuver([])
    assert maneuver.moves == ()
    maneuver = Maneuver("U F2 R'")
    assert maneuver == "U F2 R'"
    assert maneuver == [Move.U1, Move.F2, Move.R3]
    assert maneuver == Maneuver("U F2 R'")
    assert maneuver == Maneuver([Move.U1, Move.F2, Move.R3])
    assert maneuver.moves == (Move.U1, Move.F2, Move.R3)
    maneuver = Maneuver([Move.U1, Move.NONE, Move.F2, Move.NONE, Move.R3])
    assert maneuver == "U F2 R'"
    assert maneuver == [Move.U1, Move.F2, Move.R3]
    assert maneuver == Maneuver("U F2 R'")
    assert maneuver == Maneuver([Move.U1, Move.F2, Move.R3])
    assert maneuver.moves == (Move.U1, Move.F2, Move.R3)

    # equivalent
    assert Maneuver("U") == "U2 U'"
    assert Maneuver("Fw2") != "B2"
    assert Maneuver("f2") == "z2 B2"
    assert Maneuver("M'") != "R' L"
    assert Maneuver("M'") == "x R' L"
    assert Maneuver("z") != ""
    assert Maneuver("z") == "F S B'"
    assert Maneuver("z'") == "f' B"

    # reduce
    assert str(Maneuver("U U U", False)) == "U U U"
    assert str(Maneuver("U U U")) == "U'"
    assert str(Maneuver("U U U2")) == ""
    assert str(Maneuver("U D U'")) == "D"
    assert str(Maneuver("U D Uw'")) == "Dw"
    assert str(Maneuver("U D' E'")) == "y"
    assert str(Maneuver("U D' y'")) == "E"
    assert str(Maneuver("U D E2")) == "U D E2"
    maneuver = Maneuver("U D E2 U")
    assert len(maneuver) == 2
    assert maneuver == "D Uw2"
    maneuver = Maneuver("B B' U2 U' U' F2 F2 B F2 F2 U U' D2 D' D' U2 D' R U2 U2 R2 R' R' R' L' R L2 L\
                         R L2 F B' B F2 B' F' B' F' F' L2 L' F' U' U U2 U L2 D2 D D L R' L R B2 F2 F2 B'\
                         B' F' B2 B F2 F2 F2 B F2 F2 U U' D2 D' D' U2 D' U U U L2 L2 U' D U2 U2 U D' U")
    assert len(maneuver) == 11
    assert maneuver == "B U2 D' R2 B2 L F' U' F U2 D'"
    assert str(Maneuver("x' Rw y Uw' z2 Fw2 y2 Uw2 Dw2 y' Dw' y Dw Dw2 y2")) == "L D' B2 Uw2"

    # inverse
    maneuver = Maneuver("U F2 R' D B2 L' M E2 S' Uw Fw2 Rw' Dw Bw2 Lw' u f2 r' d b2 l' x y2 z'")
    cube = Cube(maneuver)
    assert repr(cube) == "YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOYO"
    cube.apply_maneuver(maneuver.inverse)
    assert repr(cube) == "WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY"
    assert maneuver == maneuver.inverse.inverse

    # container operations
    maneuver = Maneuver("U F2 R'")
    assert len(maneuver) == 3
    with pytest.raises(TypeError, match=r"Maneuver indices must be int or slice, not NoneType"):
        maneuver[None]
    with pytest.raises(IndexError, match=r"Maneuver index out of range"):
        maneuver[3]
    with pytest.raises(IndexError, match=r"Maneuver index out of range"):
        maneuver[-4]
    assert maneuver[0] == Move.U1
    assert maneuver[-1] == Move.R3
    assert maneuver[:2] == "U F2"
    assert maneuver[::-1] == "R' F2 U"
    moves = [Move.U1, Move.F2, Move.R3]
    for i, move in enumerate(maneuver):
        assert move == moves[i]
    for i, move in enumerate(reversed(maneuver)):
        assert move == moves[2-i]
    assert None not in maneuver
    assert Move.U1 in maneuver
    assert Move.U2 not in maneuver
    assert "U" in maneuver
    assert "U2" not in maneuver

    # numeric operations
    maneuver = Maneuver("U F2 R'")
    # negation
    assert str(-maneuver) == "R F2 U'"  # same as inverse
    # addition
    assert str(maneuver + "D B2 L'") == "U F2 R' D B2 L'"
    assert str("D B2 L'" + maneuver) == "D B2 L' U F2 R'"
    assert str(maneuver + [Move.L1, Move.R1, Move.U3]) == "U F2 L U'"
    assert str([Move.L1, Move.R1, Move.U3] + maneuver) == "L R F2 R'"
    # subtraction
    assert str(maneuver - "D B2 L'") == "U F2 R' L B2 D'"
    assert str("D B2 L'" - maneuver) == "D B2 L' R F2 U'"
    assert str(maneuver - [Move.L1, Move.R1, Move.U3]) == "U F2 R' U R' L'"
    assert str([Move.L1, Move.R1, Move.U3] - maneuver) == "L R U' R F2 U'"
    # integer multiplication
    assert str(maneuver * 2) == "U F2 R' U F2 R'"
    assert str(2 * maneuver) == "U F2 R' U F2 R'"
    # conjugation
    assert str(maneuver * "D B2 L'") == "U F2 R' D B2 L' R F2 U'"
    assert str("D B2 L'" * maneuver) == "D B2 L' U F2 R' L B2 D'"
    assert str(maneuver * [Move.L1, Move.R1, Move.U3]) == "U F2 L U' R F2 U'"
    assert str([Move.L1, Move.R1, Move.U3] * maneuver) == "L R F2 R' U R' L'"
    # commutator
    assert str(maneuver @ "D B2 L'") == "U F2 R' D B2 L' R F2 U' L B2 D'"
    assert str("D B2 L'" @ maneuver) == "D B2 L' U F2 R' L B2 D' R F2 U'"
    assert str(maneuver @ [Move.L1, Move.R1, Move.U3]) == "U F2 L U' R F2 R' L'"
    assert str([Move.L1, Move.R1, Move.U3] @ maneuver) == "L R F2 R' U L' F2 U'"

    # random
    with pytest.raises(TypeError, match=r"length must be int, not NoneType"):
        Maneuver.random(None)
    with pytest.raises(ValueError, match=r"length must be >= 0 \(got -1\)"):
        Maneuver.random(-1)
    maneuver = Maneuver.random(0)
    assert maneuver == ""
    maneuver = Maneuver.random(100)
    assert len(maneuver) == 100
    for first, second, third in zip(iter(maneuver[:-2]), iter(maneuver[1:-1]), iter(maneuver[2:])):
        assert first.axis != second.axis or second.axis != third.axis
