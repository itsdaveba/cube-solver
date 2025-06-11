#!/usr/bin/env python

"""Tests for `cube` package."""

import pytest
import numpy as np

from cube_solver import Cube
from cube_solver.cube import utils
from cube_solver.cube.enums import Axis, Layer, Color, Face, Cubie, Move


@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_enums(response):
    # axis
    assert hasattr(Axis, "NONE")
    assert all(axis.is_edge for axis in Axis.edge_axes())
    assert all(axis.is_corner for axis in Axis.corner_axes())
    assert len([*Axis.axes()]) == 5
    assert len([*Axis.edge_axes()]) == 3
    assert len([*Axis.corner_axes()]) == 2

    # layer
    assert hasattr(Layer, "NONE")
    assert Layer.NONE.axis == Axis.NONE
    assert Layer.NONE.perm == [[Cubie.NONE]]
    assert all(layer.is_outer for layer in Layer.outers())
    assert all(layer.is_inner for layer in Layer.inners())
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
    assert Face.NONE._cubie_slice == (slice(None), slice(None), slice(None), slice(None))
    with pytest.raises(TypeError, match=r"char must be str, not NoneType"):
        Face.from_char(None)
    with pytest.raises(ValueError, match=r"invalid face character \(got 'None'\)"):
        Face.from_char("None")
    assert all(face == Face.from_char(face.char) for face in Face)
    assert len([*Face.faces()]) == 6

    # cubie
    assert hasattr(Cubie, "NONE")
    assert Cubie.NONE.axis == Axis.NONE
    assert Cubie.NONE._index == (1, 1, 1)
    assert Cubie.NONE.faces == [Face.NONE]
    assert Cubie.CORE.faces == []
    assert Cubie.UBL.faces == [Face.UP, Face.LEFT, Face.BACK]
    assert Cubie.UBR.faces == [Face.UP, Face.BACK, Face.RIGHT]
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
    assert Move.NONE.string == "NONE"
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
    assert Move.NONE.layers == [Layer.NONE]
    assert Move.X1.layers == [Layer.R, Layer.L, Layer.M]
    assert Move.FW1.layers == [Layer.F, Layer.S]
    assert Move.U1.layers == [Layer.U]
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


def test_utils(response):
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
    permutation, parity = utils.get_permutation_array(16702, 8, True)
    assert np.all(permutation == [6, 4, 2, 1, 7, 3, 5, 0])
    assert parity is False
    # no pre-computed factorial
    permutation, parity = utils.get_permutation_array(16702, 13)
    assert np.all(permutation == [0, 1, 2, 3, 4, 8, 7, 6, 5, 12, 11, 9, 10])
    assert parity is True

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
    assert np.all(utils.get_combination_array(450, 5) == [1, 6, 8, 9, 10])

    # partial permutation array
    with pytest.raises(TypeError, match=r"coord must be int, not NoneType"):
        utils.get_partial_permutation_array(None, None)
    with pytest.raises(TypeError, match=r"n must be int, not NoneType"):
        utils.get_partial_permutation_array(-1, None)
    with pytest.raises(ValueError, match=r"n must be positive \(got 0\)"):
        utils.get_partial_permutation_array(-1, 0)
    with pytest.raises(ValueError, match=r"coord must be >= 0 \(got -1\)"):
        utils.get_partial_permutation_array(-1, 1)
    partial_permutation, combination = utils.get_partial_permutation_array(0, 1)
    assert np.all(partial_permutation == [0])
    assert np.all(combination == [0])
    partial_permutation, combination = utils.get_partial_permutation_array(450, 4)
    assert np.all(partial_permutation == [3, 0, 1, 2])
    assert np.all(combination == [1, 2, 3, 6])
    partial_permutation, combination = utils.get_partial_permutation_array(450, 13)
    assert np.all(partial_permutation == [0, 1, 2, 3, 4, 5, 6, 10, 11, 12, 7, 8, 9])
    assert np.all(combination == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])


def test_cube(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    cube = Cube()
    assert repr(cube) == "WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY"

    cube = Cube("U F2 R' D B2 L' M E2 S' Uw Fw2 Rw' Dw Bw2 Lw' u f2 r' d b2 l' x y2 z'")
    assert repr(cube) == "YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOYO"
    with pytest.raises(TypeError, match="scramble must be str or None, not int"):
        Cube(0)
    with pytest.raises(TypeError, match="maneuver must be str, not NoneType"):
        cube.apply_maneuver(None)
    with pytest.raises(TypeError, match="move must be Move, not NoneType"):
        cube.apply_move(None)
    with pytest.raises(ValueError, match="invalid move, got Move.NONE"):
        cube.apply_move(Move.NONE)

    cube = Cube(repr="YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOYO")
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
    with pytest.raises(TypeError, match="repr must be str or None, not int"):
        Cube(repr=0)
    with pytest.raises(ValueError, match="repr length must be 54, got 4"):
        Cube(repr="None")
    with pytest.raises(ValueError, match="invalid color character, got 'N'"):
        Cube(repr="YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOYN")
    match = r"invalid cubie centers, got \[Color.YELLOW, Color.GREEN, Color.ORANGE, Color.YELLOW, Color.BLUE, Color.RED\]"
    with pytest.raises(ValueError, match=match):
        Cube(repr="YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROYWOYO")
    with pytest.raises(ValueError, match=r"invalid cubie faces, got \[Face.RIGHT, Face.LEFT\]"):
        Cube(repr="YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOOO")
    with pytest.raises(ValueError, match=r"invalid cubie faces, got \[Face.UP, Face.FRONT, Face.RIGHT\]"):
        Cube(repr="YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWORBWWROWWOOY")
    with pytest.raises(ValueError, match=r"invalid cubie faces, got \[Face.FRONT, Face.FRONT, Face.FRONT\]"):
        Cube(repr="YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWGRBWWROWWOOG")
    with pytest.raises(ValueError, match=r"invalid cubie permutation"):
        Cube(repr="YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWRRBWWROWWOYY")
    with pytest.raises(ValueError, match=r"invalid cubie orientation"):
        Cube(repr="YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYYBWWROWWORO")

    cube = Cube(random_state=True)
    assert cube.coords != (0, 0, 0, 0)
    with pytest.raises(TypeError, match="random_state must be bool, not int"):
        Cube(random_state=0)
