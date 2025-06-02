#!/usr/bin/env python

"""Tests for `cube_solver` package."""

import pytest

from cube_solver import Cube, generate_scramble
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

    cube = Cube("B2 R L' F' D F' D' U B D' F2 L2 R B2 D2 L' U2 L U2 R L' U' R2 D U")
    assert str(cube) == \
        "        ---------\n        | G W O |\n        | Y W O |\n        | O R O |\n---------------------------------\n| W R"\
        + " Y | G B G | W W B | Y G R |\n| G O W | B G O | Y R G | Y B O |\n| R R R | B B O | B W W | B O G |\n--------------"\
        + "-------------------\n        | Y Y W |\n        | G Y R |\n        | Y B R |\n        ---------"

    cube.reset()
    cube.apply_maneuver("U2 B' R U2 B R' D' U F' D B' D2 U' F B2 L' R2 B' L' D2 F' B' U' B2 U2")
    assert repr(cube) == "OOGYWWYYWBBRBOBOWGGGRRGBORBBGRORYYOBWGWRBWOOYWGRRYWGYY"
    assert cube.get_coords() == [1468, 1043, 2717, 412202425]

    cube.set_coords([1607, 604, 7987, 347271465])
    assert repr(cube) == "RBROWWBRWGGRYORYYRYYOWGGWROGBGWROGBWYOWWBBOOOBGYGYRBYB"
    assert cube.get_coords(partial_corner=True, partial_edge=True) == [1607, 604, [1110, 559], [1013, 7126, 9749]]

    cube.set_coords([1440, 578, 31234, [4639, 8061, 6317]])
    assert repr(cube) == "OYROWBYRRYYRROGBGBGBWWGOYBWGOBWRGGBGWRBOBGYRORWOYYYWWO"

    cube.set_coords([-1, -1, -1, [-1, -1, -1]])
    assert cube.get_coords(partial_edge=True) == [2186, 2047, 40319, [-1, -1, -1]]

    cube.reset()
    cube.apply_maneuver("D2 U F U2 D2 L' F' B U L F B2 R2 L' D' U2 B2 R2 L2 F' U2 D B' R' B'")
    coords = cube.get_coords()
    cube.reset()
    cube.set_coords(coords)
    assert repr(cube) == "WYRBWGOGBRYGOOBBWOWRRWGYYRYYOGGRWORBYOGOBBWWWBYGRYBRGO"

    cube.reset()
    cube.apply_maneuver("L' U F2 R' L2 B L2 U' F R2 F2 B2 R2 B' U F' B D' U L U2 R F' B' D2")
    coords = cube.get_coords(partial_corner=True, partial_edge=True)
    cube.reset()
    cube.set_coords(coords)
    assert repr(cube) == "GBYRWOBRGRGYBOBBWROBWWGOGYGOGBYRWYYWRYYGBOBOOWRORYGWWR"

    scramble = generate_scramble(1000).split()
    assert len(scramble) == 1000
    for move, next_move in zip(scramble[:-1], scramble[1:]):
        assert move[0] != next_move[0]
        if move[0] in "DBL":
            assert OPPOSITE_FACE[move[0]] != next_move[0]
