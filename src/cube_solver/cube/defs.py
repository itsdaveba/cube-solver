"""Cube definitions."""
import numpy as np

NONE = -1

SIZE = 3  #: Cube size.
NUM_DIMS = 3  #: Number of dimensions.
NUM_CORNERS = 8  #: Number of corners.
NUM_EDGES = 12  #: Number of edges.
NUM_ORBIT_ELEMS = 4

# precomputed factorials and combinations
FACTORIAL = np.cumprod([1] + [*range(1, NUM_EDGES + 1)])
COMBINATION = np.zeros((NUM_EDGES + 1, NUM_ORBIT_ELEMS + 1), dtype=int)
COMBINATION[:, 0] = 1
for i in range(1, NUM_ORBIT_ELEMS + 1):
    COMBINATION[i:, i] = COMBINATION[i-1:-1, i-1].cumsum()

# coordinate sizes
CORNER_ORIENTATION_SIZE = (3 ** (NUM_CORNERS - 1))
"""Number of possible corner orientations. ``3 ^ 7``"""
EDGE_ORIENTATION_SIZE = 2 ** (NUM_EDGES - 1)
"""Number of possible edge orientations. ``2 ^ 11``"""
CORNER_PERMUTATION_SIZE = FACTORIAL[NUM_CORNERS].item()
"""Number of possible corner permutations. ``8!``"""
EDGE_PERMUTATION_SIZE = FACTORIAL[NUM_EDGES].item() // 2
"""Number of possible edge permutations. ``12! / 2``"""

PARTIAL_CORNER_PERMUTATION_SIZE = FACTORIAL[NUM_CORNERS].item() // FACTORIAL[NUM_CORNERS - NUM_ORBIT_ELEMS].item()
PARTIAL_EDGE_PERMUTATION_SIZE = FACTORIAL[NUM_EDGES].item() // FACTORIAL[NUM_EDGES - NUM_ORBIT_ELEMS].item()

NUM_CUBE_POSITIONS = CORNER_ORIENTATION_SIZE * EDGE_ORIENTATION_SIZE * CORNER_PERMUTATION_SIZE * EDGE_PERMUTATION_SIZE
"""Number of all possible cube positions. ``3 ^ 7 * 2 ^ 11 * 8! * 12! / 2``"""
