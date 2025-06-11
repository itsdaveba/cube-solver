import numpy as np

NONE = -1

SIZE = 3
NUM_CORNERS = 8
NUM_EDGES = 12
NUM_AXIS_ELEMS = 4

# TODO print to see how many times this is run
# precomputed factorials and combinations
FACTORIAL = np.cumprod([1] + list(range(1, NUM_EDGES + 1)))
COMBINATION = np.zeros((np.max([NUM_CORNERS, NUM_EDGES]) + 1, NUM_AXIS_ELEMS + 1), dtype=int)
COMBINATION[:, 0] = 1
for i in range(1, NUM_AXIS_ELEMS + 1):
    COMBINATION[i:, i] = COMBINATION[i-1:-1, i-1].cumsum()

# coordinate sizes
CORNER_ORIENTATION_SIZE = 3 ** (NUM_CORNERS - 1)
EDGE_ORIENTATION_SIZE = 2 ** (NUM_EDGES - 1)
CORNER_PERMUTATION_SIZE = FACTORIAL[NUM_CORNERS]
EDGE_PERMUTATION_SIZE = FACTORIAL[NUM_EDGES] // 2
