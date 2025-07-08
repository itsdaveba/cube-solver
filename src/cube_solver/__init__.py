"""Top-level package for Cube Solver."""

__author__ = """Dave Barragan"""
__email__ = 'itsdaveba@gmail.com'
__version__ = '1.1.0'

from .cube import Cube, apply_move, apply_maneuver, Move, Maneuver
from .solver import BaseSolver, DummySolver, Korf, Thistlethwaite, Kociemba

__all__ = ["Cube", "apply_move", "apply_maneuver", "Move", "Maneuver",
           "BaseSolver", "DummySolver", "Korf", "Thistlethwaite", "Kociemba"]
