"""Cube definitions."""
import numpy as np

NONE = -1

SIZE = 2  #: Cube size.
NUM_DIMS = 3  #: Number of dimensions.
NUM_CORNERS = 8  #: Number of corners.

# precomputed factorials and combinations
FACTORIAL = np.cumprod([1] + [*range(1, NUM_CORNERS + 1)])

# coordinate sizes
CORNER_ORIENTATION_SIZE = (3 ** (NUM_CORNERS - 2))
"""Number of possible corner orientations. ``3 ^ 6``"""
CORNER_PERMUTATION_SIZE = FACTORIAL[NUM_CORNERS - 1].item()
"""Number of possible corner permutations. ``7!``"""

NUM_CUBE_POSITIONS = CORNER_ORIENTATION_SIZE * CORNER_PERMUTATION_SIZE
"""Number of all possible cube positions. ``3 ^ 6 * 7!``"""
