#!/usr/bin/env python

"""Tests for `cube_solver` package."""

import pytest

from cube_solver import Cube
from cube_solver.cube.enums import Axis, Layer, Color, Face, Cubie, Move


@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_enums(response):
    assert hasattr(Axis, "NONE")
    assert len([*Axis.axes()]) == 5

    assert hasattr(Layer, "NONE")
    assert Layer.NONE.axis == Axis.NONE
    assert Layer.NONE.perm == [[Cubie.NONE]]
    assert all(layer.is_outer for layer in Layer.outers())
    assert all(layer.is_inner for layer in Layer.inners())
    assert len([*Layer.layers()]) == 9
    assert len([*Layer.outers()]) == 6
    assert len([*Layer.inners()]) == 3

    assert hasattr(Color, "NONE")
    assert Color.NONE.char == 'N'
    assert all(color == Color.from_char(color.char) for color in Color)
    with pytest.raises(TypeError, match="char must be str, not NoneType"):
        Color.from_char(None)
    with pytest.raises(ValueError, match="invalid color character, got 'None'"):
        Color.from_char("None")
    assert len([*Color.colors()]) == 6

    assert hasattr(Face, "NONE")
    assert Face.NONE.char == 'N'
    assert Face.NONE.axis == Axis.NONE
    assert Face.NONE.opposite == Face.NONE
    assert Face.NONE._cubie_slice == (slice(None), slice(None), slice(None), slice(None))
    assert all(face == Face.from_char(face.char) for face in Face)
    with pytest.raises(TypeError, match="char must be str, not NoneType"):
        Face.from_char(None)
    with pytest.raises(ValueError, match="invalid face character, got 'None'"):
        Face.from_char("None")
    assert len([*Face.faces()]) == 6

    assert hasattr(Cubie, "NONE")
    assert Cubie.NONE.axis == Axis.NONE
    assert Cubie.NONE._index == (slice(None), slice(None), slice(None))
    assert Cubie.NONE.faces == [Face.NONE]
    assert all(cubie.is_corner for cubie in Cubie.corners())
    assert all(cubie.is_edge for cubie in Cubie.edges())
    assert all(cubie.is_center for cubie in Cubie.centers())
    assert all(cubie == Cubie.from_faces(cubie.faces) for cubie in Cubie)
    with pytest.raises(TypeError, match="faces must be list, not NoneType"):
        Cubie.from_faces(None)
    with pytest.raises(TypeError, match="faces elements must be Face, not NoneType"):
        Cubie.from_faces([None])
    with pytest.raises(ValueError, match=r"invalid cubie faces, got \[Face.NONE, Face.NONE\]"):
        Cubie.from_faces([Face.NONE, Face.NONE])
    with pytest.raises(ValueError, match="faces length must be at most 3, got 4"):
        Cubie.from_faces([Face.NONE, Face.NONE, Face.NONE, Face.NONE])
    assert len([*Cubie.cubies()]) == 27
    assert len([*Cubie.corners()]) == 8
    assert len([*Cubie.edges()]) == 12
    assert len([*Cubie.centers()]) == 6

    assert hasattr(Move, "NONE")
    assert Move.NONE.string == "NONE"
    assert Move.NONE.axis == Axis.NONE
    assert Move.U1.axis == Axis.Y
    assert Move.X1.axis == Axis.X
    assert Move.NONE.layers == [Layer.NONE]
    assert Move.X1.layers == [Layer.R, Layer.L, Layer.M]
    assert Move.FW1.layers == [Layer.F, Layer.S]
    assert Move.NONE.shifts == [0]
    assert Move.X1.shifts == [1, -1, -1]
    assert Move.FW1.shifts == [1, 1]
    assert all(move.is_face for move in Move.face_moves())
    assert all(move.is_slice for move in Move.slice_moves())
    assert all(move.is_wide for move in Move.wide_moves())
    assert all(move.is_rotation for move in Move.rotations())
    assert all(move == Move.from_string(move.string) for move in Move)
    assert all(move == Move.from_string(move.string[0].lower() + move.string[2:]) for move in Move.wide_moves())
    with pytest.raises(TypeError, match="string must be str, not NoneType"):
        Move.from_string(None)
    with pytest.raises(ValueError, match="invalid move string, got 'None'"):
        Move.from_string("None")
    assert len([*Move.moves()]) == 54
    assert len([*Move.face_moves()]) == 18
    assert len([*Move.slice_moves()]) == 9
    assert len([*Move.wide_moves()]) == 18
    assert len([*Move.rotations()]) == 9


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
