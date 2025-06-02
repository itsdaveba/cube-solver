import numpy as np

# visual representation
SIZE = 3  # 3x3 cube
COLORS = "WGRYBO"  # up, front, right, down, back, left
FACES = "UFRDBL"  # up, front, right, down, back, left
REPR_ORDER = [0, 5, 1, 2, 4, 3]  # up, left, front, right, back, down - for repr() and str()

# scramble
MOVE_COUNT_STR = ["'", "", "2"]  # example: U0 -> U', U1 -> U, U2 -> U2
OPPOSITE_FACE = {face: opp for face, opp in zip("UFRDBL", "DBLUFR")}
NEXT_BASE_MOVES = {face: [f for f in FACES if f != face and (f != OPPOSITE_FACE[face] or face in "UFR")] for face in FACES}
NEXT_BASE_MOVES.update({None: FACES})
ALL_MOVES = [face + count_str for face in FACES for count_str in MOVE_COUNT_STR]
ALL_MOVES_INDEX = {move: i for i, move in enumerate(ALL_MOVES)}
NEXT_MOVES = {None: ALL_MOVES}
for move in ALL_MOVES:
    next_moves = [m + cs for m in NEXT_BASE_MOVES[move[0]] for cs in MOVE_COUNT_STR]
    NEXT_MOVES.update({move: next_moves})

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

NUM_AXIS_ELEMS = 4
NUM_PERM_AXIS = 24
CORNER_AXIS_OFFSET = [0, 4]
EDGE_AXIS_OFFSET = [16, 12, 8]
FACTORIAL = np.cumprod(range(1, NUM_EDGES + 1))

COORD_MOVES = {
    0: [[0, 4, 1, 5], [8, 13, 9, 12]],
    1: [[1, 7, 3, 5], [9, 19, 11, 18]],
    2: [[1, 4, 2, 7], [13, 17, 15, 19]],
    3: [[2, 6, 3, 7], [10, 14, 11, 15]],
    4: [[0, 6, 2, 4], [8, 16, 10, 17]],
    5: [[0, 5, 3, 6], [12, 18, 14, 16]]
}
COORD_CUBIE_INDEX = [
    (0, 0, 0), (0, 2, 2), (2, 0, 2), (2, 2, 0), (0, 0, 2), (0, 2, 0), (2, 0, 0), (2, 2, 2),  # corners
    (0, 0, 1), (0, 2, 1), (2, 0, 1), (2, 2, 1), (0, 1, 0), (0, 1, 2),  # edges
    (2, 1, 0), (2, 1, 2), (1, 0, 0), (1, 0, 2), (1, 2, 0), (1, 2, 2)
]

# solver with coord representation
SOLVED_PARTIAL_CORNER_PERMUTATION = (0, 1656)
SOLVED_PARTIAL_EDGE_PERMUTATION = (11856, 1656, 0)
SOLVED_PARTIAL_COORDS = (0, 0, SOLVED_PARTIAL_CORNER_PERMUTATION, SOLVED_PARTIAL_EDGE_PERMUTATION)
SOLVED_REPR = "".join(COLORS[r] * SIZE * SIZE for r in REPR_ORDER)

CORNER_ORIENTATION_SIZE = 3 ** (NUM_CORNERS - 1)
EDGE_ORIENTATION_SIZE = 2 ** (NUM_EDGES - 1)
CORNER_PERMUTATION_SIZE = FACTORIAL[NUM_CORNERS - 1]
EDGE_PERMUTATION_SIZE = FACTORIAL[NUM_EDGES - 1] // 2
PARTIAL_CORNER_PERMUTATION_SIZE = FACTORIAL[NUM_CORNERS - 1] // FACTORIAL[NUM_CORNERS - NUM_AXIS_ELEMS - 1]
PARTIAL_EDGE_PERMUTATION_SIZE = FACTORIAL[NUM_EDGES - 1] // FACTORIAL[NUM_EDGES - NUM_AXIS_ELEMS - 1]
COORDS_SIZES = [CORNER_ORIENTATION_SIZE, EDGE_ORIENTATION_SIZE, [CORNER_PERMUTATION_SIZE, PARTIAL_CORNER_PERMUTATION_SIZE], [EDGE_PERMUTATION_SIZE, PARTIAL_EDGE_PERMUTATION_SIZE]]
COORDS_SIZES = [CORNER_ORIENTATION_SIZE, EDGE_ORIENTATION_SIZE, CORNER_PERMUTATION_SIZE, EDGE_PERMUTATION_SIZE]


# thistlethwaite
NUM_PHASES = 4
NUM_THREADS = 6

PHASE_MOVES = []
PHASE_NEXT_MOVES = []
DOUBLE_STR = MOVE_COUNT_STR[2]
RESTRICT = ["", "FB", "FRBL", FACES]
for i in range(NUM_PHASES):
    phase_moves = [face + count_str for face in FACES for count_str in (DOUBLE_STR if face in RESTRICT[i] else MOVE_COUNT_STR)]
    phase_next_moves = {None: phase_moves}
    for move in phase_moves:
        next_moves = [m + cs for m in NEXT_BASE_MOVES[move[0]] for cs in (DOUBLE_STR if m in RESTRICT[i] else MOVE_COUNT_STR)]
        phase_next_moves.update({move: next_moves})
    PHASE_MOVES.append(phase_moves)
    PHASE_NEXT_MOVES.append(phase_next_moves)

COMBINATION_8C4 = FACTORIAL[8 - 1] // FACTORIAL[4 - 1]
PHASE_TABLE_SIZES = [[(EDGE_ORIENTATION_SIZE,)]]
PHASE_TABLE_SIZES.append([(CORNER_ORIENTATION_SIZE, PARTIAL_EDGE_PERMUTATION_SIZE // NUM_PERM_AXIS)])
PHASE_TABLE_SIZES.append([(COMBINATION_8C4, COMBINATION_8C4, NUM_THREADS)])
PHASE_TABLE_SIZES.append([(NUM_PERM_AXIS, NUM_PERM_AXIS // NUM_THREADS, NUM_PERM_AXIS, NUM_PERM_AXIS, NUM_PERM_AXIS // 2)])

THREAD_PERM_GROUP = [1, 0, 4, 5, 2, 3]
THREAD_SELECTOR = np.full((NUM_THREADS, NUM_THREADS), EMPTY, dtype=int)
filter = np.eye(NUM_THREADS, dtype=bool)
for i in range(NUM_THREADS):
    THREAD_SELECTOR[filter] = i
    filter = filter[THREAD_PERM_GROUP] if i % 2 == 0 else np.rot90(filter, 2)
CORNER_THREAD = np.full((NUM_THREADS, NUM_THREADS), EMPTY, dtype=int)
for i, y in enumerate(THREAD_SELECTOR):
    CORNER_THREAD[y, range(NUM_THREADS)] = i
CORNER_THREAD = np.vstack((CORNER_THREAD, CORNER_THREAD[THREAD_PERM_GROUP]))
CORNER_THREAD = np.hstack((CORNER_THREAD, CORNER_THREAD[:, THREAD_PERM_GROUP]))
CORNER_THREAD = np.vstack((CORNER_THREAD, np.flipud(CORNER_THREAD)))
CORNER_THREAD = np.hstack((CORNER_THREAD, np.fliplr(CORNER_THREAD)))
