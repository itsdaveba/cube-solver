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
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string
    cube = Cube(size=2)
    assert str(cube) == \
        "      -------\n      | W W | \n      | W W | \n-------------------------\n| O O | G G | R R | B B | \n| O O | G G | "\
        + "R R | B B | \n-------------------------\n      | Y Y | \n      | Y Y | \n      -------"
    cube = Cube()
    assert str(cube) == \
        "        ---------\n        | W W W | \n        | W W W | \n        | W W W | \n---------------------------------\n| "\
        + "O O O | G G G | R R R | B B B | \n| O O O | G G G | R R R | B B B | \n| O O O | G G G | R R R | B B B | \n--------"\
        + "-------------------------\n        | Y Y Y | \n        | Y Y Y | \n        | Y Y Y | \n        ---------"

    cube = Cube("B L B L F' L U' R' D2 U' B2 F R' B' F D U F2 D U' L2 R D F B")
    assert str(cube) == \
        "        ---------\n        | B Y B | \n        | B W R | \n        | R B W | \n---------------------------------\n| "\
        + "O Y Y | B O B | O G W | R G Y | \n| R O W | R G Y | O R B | W B Y | \n| O G R | Y W W | R G O | W B Y | \n--------"\
        + "-------------------------\n        | G O G | \n        | O Y W | \n        | G R G | \n        ---------"

    cube.reset()
    cube.apply_maneuver("D F' U B2 R' B' F D2 U' L U B' R D' U F D R' F' U2 F' L F2 R2 D2")
    assert repr(cube) == "WOYWWWBBRBBROOWRGWYRGOGBGRGYROYROYROBYRBBGBYWOYOWYGGGW"

    cube.reset()
    cube.apply_maneuver("U D F2 R D F' D B L F R2 F2 D2 R F' L2 U2 D' F2 L2 U2 F2 L' D2 R'")
    assert repr(cube) == "OYBBWRYBWGRBOOYWWWRWORGWOWRGGWRROYYBRGYBBGOOGBOGGYBRYY"

    cube = Cube("F' U B R2 D B' R2 F2 B2 L B D2 F' L2 F L2 D B' D L' F' L F' U F2", representation="cubie")
    assert str(cube) == \
        "        ---------\n        | B B W | \n        | W W G | \n        | G Y G | \n---------------------------------\n| "\
        + "Y R Y | R O R | W R B | O Y R | \n| B O O | W G G | W R W | B B R | \n| G Y Y | O G B | R R B | O O W | \n--------"\
        + "-------------------------\n        | G O W | \n        | G Y Y | \n        | O B Y | \n        ---------"

    cube.reset()
    cube.apply_maneuver("B R' L F2 L B2 L' B R B R D2 R' L2 U2 L2 U B U' D B L D2 F' D2")
    assert repr(cube) == "BWRGWYGRRYYOGOYORWYGWBGWGOWGOGBRBBBBYGRRBOWRBOWRWYOYYO"

    scramble = Cube.generate_scramble(1000).split()
    for move, next_move in zip(scramble[:-1], scramble[1:]):
        assert move[0] != next_move[0]
        if move[0] in "UFR":
            assert OPPOSITE_FACE[move[0]] != next_move[0]
