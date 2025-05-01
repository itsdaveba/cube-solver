import numpy as np

# visual representation
SIZE = 3  # 3x3 cube
COLORS = "WGRYBO"  # up, front, right, down, back, left
FACES = "UFRDBL"  # up, front, right, down, back, left
REPR_ORDER = [0, 5, 1, 2, 4, 3]  # up, left, front, right, back, down - for repr() and str()

# scramble
MOVE_COUNT_STR = ["'", "", "2"]  # example: U0 -> U', U1 -> U, U2 -> U2
OPPOSITE_FACE = {face: opp for face, opp in zip("UFRDBL", "DBLUFR")}
NEXT_BASE_MOVES = {face: set(FACES) - {face} - ({OPPOSITE_FACE[face]} if face in "DBL" else {None}) for face in FACES}
NEXT_BASE_MOVES.update({None: FACES})
ALL_MOVES = [face + count_str for face in FACES for count_str in MOVE_COUNT_STR]
ALL_MOVES_INDEX = {move: i for i, move in enumerate(ALL_MOVES)}
NEXT_MOVES = {move: [m + cs for m in NEXT_BASE_MOVES[move[0]] for cs in MOVE_COUNT_STR] for move in ALL_MOVES}
NEXT_MOVES.update({None: ALL_MOVES})

# face representation
FACE_MOVES = {
    "U": [([0, 0, 0, 0], [0, 0, 2, 2], [0, 2, 2, 0]),  # corners
          ([1, 5, 4, 2], [0, 0, 0, 0], [0, 0, 0, 0]),
          ([1, 5, 4, 2], [0, 0, 0, 0], [2, 2, 2, 2]),
          ([0, 0, 0, 0], [0, 1, 2, 1], [1, 2, 1, 0]),  # edges
          ([1, 5, 4, 2], [0, 0, 0, 0], [1, 1, 1, 1])],
    "F": [([0, 2, 3, 5], [2, 0, 0, 2], [0, 0, 2, 2]),
          ([0, 2, 3, 5], [2, 2, 0, 0], [2, 0, 0, 2]),
          ([1, 1, 1, 1], [0, 0, 2, 2], [0, 2, 2, 0]),
          ([0, 2, 3, 5], [2, 1, 0, 1], [1, 0, 1, 2]),
          ([1, 1, 1, 1], [0, 1, 2, 1], [1, 2, 1, 0])],
    "R": [([0, 4, 3, 1], [0, 2, 0, 0], [2, 0, 2, 2]),
          ([0, 4, 3, 1], [2, 0, 2, 2], [2, 0, 2, 2]),
          ([2, 2, 2, 2], [0, 0, 2, 2], [0, 2, 2, 0]),
          ([0, 4, 3, 1], [1, 1, 1, 1], [2, 0, 2, 2]),
          ([2, 2, 2, 2], [0, 1, 2, 1], [1, 2, 1, 0])],
    "D": [([1, 2, 4, 5], [2, 2, 2, 2], [0, 0, 0, 0]),
          ([1, 2, 4, 5], [2, 2, 2, 2], [2, 2, 2, 2]),
          ([3, 3, 3, 3], [0, 0, 2, 2], [0, 2, 2, 0]),
          ([1, 2, 4, 5], [2, 2, 2, 2], [1, 1, 1, 1]),
          ([3, 3, 3, 3], [0, 1, 2, 1], [1, 2, 1, 0])],
    "B": [([0, 5, 3, 2], [0, 2, 2, 0], [0, 0, 2, 2]),
          ([0, 5, 3, 2], [0, 0, 2, 2], [2, 0, 0, 2]),
          ([4, 4, 4, 4], [0, 0, 2, 2], [0, 2, 2, 0]),
          ([0, 5, 3, 2], [0, 1, 2, 1], [1, 0, 1, 2]),
          ([4, 4, 4, 4], [0, 1, 2, 1], [1, 2, 1, 0])],
    "L": [([0, 1, 3, 4], [0, 0, 0, 2], [0, 0, 0, 2]),
          ([0, 1, 3, 4], [2, 2, 2, 0], [0, 0, 0, 2]),
          ([5, 5, 5, 5], [0, 0, 2, 2], [0, 2, 2, 0]),
          ([0, 1, 3, 4], [1, 1, 1, 1], [0, 0, 0, 2]),
          ([5, 5, 5, 5], [0, 1, 2, 1], [1, 2, 1, 0])]
}

# array representation
ARRAY_MOVES = {
    "U": [[0, 2, 8, 6], [9, 45, 36, 18], [11, 47, 38, 20], [1, 5, 7, 3], [10, 46, 37, 19]],
    "F": [[6, 18, 29, 53], [8, 24, 27, 47], [9, 11, 17, 15], [7, 21, 28, 50], [10, 14, 16, 12]],
    "R": [[2, 42, 29, 11], [8, 36, 35, 17], [18, 20, 26, 24], [5, 39, 32, 14], [19, 23, 25, 21]],
    "D": [[15, 24, 42, 51], [17, 26, 44, 53], [27, 29, 35, 33], [16, 25, 43, 52], [28, 32, 34, 30]],
    "B": [[0, 51, 35, 20], [2, 45, 33, 26], [36, 38, 44, 42], [1, 48, 34, 23], [37, 41, 43, 39]],
    "L": [[0, 9, 27, 44], [6, 15, 33, 38], [45, 47, 53, 51], [3, 12, 30, 41], [46, 50, 52, 48]]
}

# cubie representation
AXES = {face: axis for face, axis in zip(FACES, [0, 1, 2, 0, 1, 2])}  # up, front, right, down, back, left
SWAP = [[0, 2, 1], [2, 1, 0,], [1, 0, 2]]  # swap cubie along axis
CUBIE_IDX = {
    "U": (0, slice(None), slice(None)),
    "F": (slice(None), 2, slice(None)),
    "R": (slice(None), slice(None, None, -1), 2),
    "D": (2, slice(None, None, -1), slice(None)),
    "B": (slice(None), 0, slice(None, None, -1)),
    "L": (slice(None), slice(None), 0),
}

# coord representation
NUM_CORNERS = 8
NUM_EDGES = 12
EMPTY = -1

CORNER_CYCLE = [0 if i < 4 else 1 for i in range(8)]
EDGE_AXIS = {i: 2 if i < 12 else 1 if i < 16 else 0 for i in range(8, 20)}
EDGE_AXIS.update({EMPTY: EMPTY})
NUM_EDGES_AXIS = 4
EDGE_AXIS_OFFSET = [16, 12, 8]

PERM_EDGES_AXIS = 24
COMB_EDGES_AXIS = np.zeros((NUM_EDGES_AXIS, NUM_EDGES), dtype=int)
COMB_EDGES_AXIS[0] = range(12)
for i in range(1, NUM_EDGES_AXIS):
    COMB_EDGES_AXIS[i, i:] = COMB_EDGES_AXIS[i-1, i-1:-1].cumsum()

COORD_MOVES = {
    "U": [[0, 4, 1, 5], [8, 13, 9, 12]],
    "F": [[1, 7, 3, 5], [9, 19, 11, 18]],
    "R": [[1, 4, 2, 7], [13, 17, 15, 19]],
    "D": [[2, 6, 3, 7], [10, 14, 11, 15]],
    "B": [[0, 6, 2, 4], [8, 16, 10, 17]],
    "L": [[0, 5, 3, 6], [12, 18, 14, 16]]
}
COORD_CUBIE_INDEX = [
    (0, 0, 0), (0, 2, 2), (2, 0, 2), (2, 2, 0), (0, 0, 2), (0, 2, 0), (2, 0, 0), (2, 2, 2),  # corners
    (0, 0, 1), (0, 2, 1), (2, 0, 1), (2, 2, 1), (0, 1, 0), (0, 1, 2),  # edges
    (2, 1, 0), (2, 1, 2), (1, 0, 0), (1, 0, 2), (1, 2, 0), (1, 2, 2)
]

# solver with coord representation
SOLVED_PARTIAL_EDGE_PERMUTATION = (11856, 1656, 0)
SOLVED_PARTIAL_COORD = (0, 0, 0, SOLVED_PARTIAL_EDGE_PERMUTATION)
SOLVED_REPR = "".join([COLORS[r] * SIZE * SIZE for r in REPR_ORDER])

CORNER_ORIENTATION_SIZE = 3 ** (NUM_CORNERS - 1)
EDGE_ORIENTATION_SIZE = 2 ** (NUM_EDGES - 1)
CORNER_PERMUTATION_SIZE = 8 * 7 * 6 * 5 * 4 * 3 * 2
EDGE_PERMUTATION_SIZE = 12 * 11 * 10 * 9 * 8 * 7 * 6 * 5 * 4 * 3
PARTIAL_EDGE_PERMUTATION_SIZE = 12 * 11 * 10 * 9
COORDS_SIZES = [CORNER_ORIENTATION_SIZE, EDGE_ORIENTATION_SIZE, CORNER_PERMUTATION_SIZE, PARTIAL_EDGE_PERMUTATION_SIZE]
