#!/usr/bin/env python

"""Tests for `cube_solver` package."""

import pytest

from cube_solver import Cube
from cube_solver.constants import OPPOSITE_FACE


@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    cube = Cube()
    assert str(cube) == \
        "        ---------\n        | W W W |\n        | W W W |\n        | W W W |\n---------------------------------\n| O O"\
        + " O | G G G | R R R | B B B |\n| O O O | G G G | R R R | B B B |\n| O O O | G G G | R R R | B B B |\n--------------"\
        + "-------------------\n        | Y Y Y |\n        | Y Y Y |\n        | Y Y Y |\n        ---------"

    cube = Cube("B L B L F' L U' R' D2 U' B2 F R' B' F D U F2 D U' L2 R D F B", representation="face")
    assert str(cube) == \
        "        ---------\n        | B Y B |\n        | B W R |\n        | R B W |\n---------------------------------\n| O Y"\
        + " Y | B O B | O G W | R G Y |\n| R O W | R G Y | O R B | W B Y |\n| O G R | Y W W | R G O | W B Y |\n--------------"\
        + "-------------------\n        | G O G |\n        | O Y W |\n        | G R G |\n        ---------"

    cube.reset()
    cube.apply_maneuver("D F' U B2 R' B' F D2 U' L U B' R D' U F D R' F' U2 F' L F2 R2 D2")
    assert repr(cube) == "WOYWWWBBRBBROOWRGWYRGOGBGRGYROYROYROBYRBBGBYWOYOWYGGGW"

    cube = Cube("R2 L U F' R D U F U2 R2 U2 R2 U2 L2 U' R' F2 R' D2 U B' D2 R' B' U2", representation="array")
    assert str(cube) == \
        "        ---------\n        | Y B O |\n        | R W R |\n        | Y Y G |\n---------------------------------\n| B Y"\
        + " R | G B O | Y W G | W O O |\n| G O G | W G O | W R G | O B R |\n| G B R | W Y B | W B R | B Y R |\n--------------"\
        + "-------------------\n        | B G O |\n        | R Y W |\n        | W O Y |\n        ---------"

    cube.reset()
    cube.apply_maneuver("U D F2 R D F' D B L F R2 F2 D2 R F' L2 U2 D' F2 L2 U2 F2 L' D2 R'")
    assert repr(cube) == "OYBBWRYBWGRBOOYWWWRWORGWOWRGGWRROYYBRGYBBGOOGBOGGYBRYY"

    cube = Cube("F' U B R2 D B' R2 F2 B2 L B D2 F' L2 F L2 D B' D L' F' L F' U F2", representation="cubie")
    assert str(cube) == \
        "        ---------\n        | B B W |\n        | W W G |\n        | G Y G |\n---------------------------------\n| Y R"\
        + " Y | R O R | W R B | O Y R |\n| B O O | W G G | W R W | B B R |\n| G Y Y | O G B | R R B | O O W |\n--------------"\
        + "-------------------\n        | G O W |\n        | G Y Y |\n        | O B Y |\n        ---------"

    cube.reset()
    cube.apply_maneuver("B R' L F2 L B2 L' B R B R D2 R' L2 U2 L2 U B U' D B L D2 F' D2")
    assert repr(cube) == "BWRGWYGRRYYOGOYORWYGWBGWGOWGOGBRBBBBYGRRBOWRBOWRWYOYYO"

    cube = Cube("B2 R L' F' D F' D' U B D' F2 L2 R B2 D2 L' U2 L U2 R L' U' R2 D U", representation="coord")
    assert str(cube) == \
        "        ---------\n        | G W O |\n        | Y W O |\n        | O R O |\n---------------------------------\n| W R"\
        + " Y | G B G | W W B | Y G R |\n| G O W | B G O | Y R G | Y B O |\n| R R R | B B O | B W W | B O G |\n--------------"\
        + "-------------------\n        | Y Y W |\n        | G Y R |\n        | Y B R |\n        ---------"

    cube.reset()
    cube.apply_maneuver("U2 B' R U2 B R' D' U F' D B' D2 U' F B2 L' R2 B' L' D2 F' B' U' B2 U2")
    assert repr(cube) == "OOGYWWYYWBBRBOBOWGGGRRGBORBBGRORYYOBWGWRBWOOYWGRRYWGYY"

    assert cube.get_coords() == (1468, 1043, 2717, 206101212)
    cube.set_coords((1607, 604, 7987, 173635732))
    assert repr(cube) == "RBROWWBRWGGRYORYYRYYOWGGWROGBGWROGBWYOWWBBOOOBGYGYRBYB"

    assert cube.get_coords(partial_edge=True) == (1607, 604, 7987, (1013, 7126, 9749))
    cube.set_coords((1440, 578, 31234, (4639, 8061, 6317)), partial_edge=True)
    assert repr(cube) == "OYROWBYRRYYRROGBGBGBWWGOYBWGOBWRGGBGWRBOBGYRORWOYYYWWO"

    cube.set_coords((-1, -1, -1, (-1, -1, -1)), partial_edge=True)
    assert cube.get_coords(partial_edge=True) == (2186, 2047, 40319, (-1, -1, -1))

    cube.reset()
    cube.apply_maneuver("D2 U F U2 D2 L' F' B U L F B2 R2 L' D' U2 B2 R2 L2 F' U2 D B' R' B'")
    coords = cube.get_coords()
    cube.reset()
    cube.set_coords(coords)
    assert repr(cube) == "WYRBWGOGBRYGOOBBWOWRRWGYYRYYOGGRWORBYOGOBBWWWBYGRYBRGO"

    cube.reset()
    cube.apply_maneuver("L' U F2 R' L2 B L2 U' F R2 F2 B2 R2 B' U F' B D' U L U2 R F' B' D2")
    coords = cube.get_coords(partial_edge=True)
    cube.reset()
    cube.set_coords(coords, partial_edge=True)
    assert repr(cube) == "GBYRWOBRGRGYBOBBWROBWWGOGYGOGBYRWYYWRYYGBOBOOWRORYGWWR"

    scramble = Cube.generate_scramble(1000).split()
    assert len(scramble) == 1000
    for move, next_move in zip(scramble[:-1], scramble[1:]):
        assert move[0] != next_move[0]
        if move[0] in "DBL":
            assert OPPOSITE_FACE[move[0]] != next_move[0]
