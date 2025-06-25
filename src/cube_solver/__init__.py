"""Top-level package for Cube Solver."""

__author__ = """Dave Barragan"""
__email__ = 'itsdaveba@gmail.com'
__version__ = '1.0.1'

from .cube import Cube, Move, Maneuver, apply_move, apply_maneuver
from .solver import DummySolver, Korf, Thistlethwaite, Kociemba

__all__ = ["Cube", "apply_move", "apply_maneuver", "Move", "Maneuver",
           "DummySolver", "Korf", "Thistlethwaite", "Kociemba"]
