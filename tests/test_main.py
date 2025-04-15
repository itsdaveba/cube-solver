#!/usr/bin/env python

"""Tests for `cube_solver` package."""

import pytest

from cube_solver import main
from cube_solver.constants import OPPOSITE


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
    cube = main.Cube(size=2)
    assert str(cube) == \
        "      -------\n      | W W | \n      | W W | \n-------------------------\n| O O | G G | R R | B B | \n| O O | G G | "\
        + "R R | B B | \n-------------------------\n      | Y Y | \n      | Y Y | \n      -------"
    cube = main.Cube()
    assert str(cube) == \
        "        ---------\n        | W W W | \n        | W W W | \n        | W W W | \n---------------------------------\n| "\
        + "O O O | G G G | R R R | B B B | \n| O O O | G G G | R R R | B B B | \n| O O O | G G G | R R R | B B B | \n--------"\
        + "-------------------------\n        | Y Y Y | \n        | Y Y Y | \n        | Y Y Y | \n        ---------"

    cube.apply_maneuver("D F' U B2 R' B' F D2 U' L U B' R D' U F D R' F' U2 F' L F2 R2 D2")
    assert repr(cube) == "WOYWWWBBRBBROOWRGWYRGOGBGRGYROYROYROBYRBBGBYWOYOWYGGGW"

    cube = main.Cube("B R' L F2 L B2 L' B R B R D2 R' L2 U2 L2 U B U' D B L D2 F' D2")
    assert repr(cube) == "BWRGWYGRRYYOGOYORWYGWBGWGOWGOGBRBBBBYGRRBOWRBOWRWYOYYO"

    cube = main.Cube("U D F2 R D F' D B L F R2 F2 D2 R F' L2 U2 D' F2 L2 U2 F2 L' D2 R'")
    assert repr(cube) == "OYBBWRYBWGRBOOYWWWRWORGWOWRGGWRROYYBRGYBBGOOGBOGGYBRYY"

    scramble = main.Cube.generate_scramble(1000).split()
    for move, next_move in zip(scramble[:-1], scramble[1:]):
        assert move[0] != next_move[0]
        if move[0] in "ULF":
            assert OPPOSITE[move[0]] != next_move[0]

    cube = main.Cube("R L2 F'")
    assert cube.solve() == "F L2 R'"

    cube = main.Cube("L F' L")
    assert cube.solve() == "L' F L'"

    cube = main.Cube("D L2 U'")
    assert cube.solve() == "U L2 D'"
