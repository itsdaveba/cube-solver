COLORS = "WGRYBO"  # up, front, right, down, back, left
FACES = "UFRDBL"  # up, front, right, down, back, left
AXES = {face: axis for face, axis in zip(FACES, [0, 1, 2, 0, 1, 2])}  # up, front, right, down, back, left
OPPOSITE_FACE = {face: opp for face, opp in zip("UFRDBL", "DBLUFR")}
MOVE_COUNT_STR = ["'", "", "2"]  # example: U0 -> U', U1 -> U, U2 -> U2
REPR_ORDER = [0, 5, 1, 2, 4, 3]  # up, left, front, right, back, down - for repr() and str()

FACE_MOVES = {
    "U": [([0, 0, 0, 0], [0, 0, 2, 2], [0, 2, 2, 0]),  # corner
          ([1, 5, 4, 2], [0, 0, 0, 0], [0, 0, 0, 0]),
          ([1, 5, 4, 2], [0, 0, 0, 0], [2, 2, 2, 2]),
          ([0, 0, 0, 0], [0, 1, 2, 1], [1, 2, 1, 0]),  # edge
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

ARRAY_MOVES = {
    "U": [[0, 2, 8, 6], [9, 45, 36, 18], [11, 47, 38, 20], [1, 5, 7, 3], [10, 46, 37, 19]],
    "F": [[6, 18, 29, 53], [8, 24, 27, 47], [9, 11, 17, 15], [7, 21, 28, 50], [10, 14, 16, 12]],
    "R": [[2, 42, 29, 11], [8, 36, 35, 17], [18, 20, 26, 24], [5, 39, 32, 14], [19, 23, 25, 21]],
    "D": [[15, 24, 42, 51], [17, 26, 44, 53], [27, 29, 35, 33], [16, 25, 43, 52], [28, 32, 34, 30]],
    "B": [[0, 51, 35, 20], [2, 45, 33, 26], [36, 38, 44, 42], [1, 48, 34, 23], [37, 41, 43, 39]],
    "L": [[0, 9, 27, 44], [6, 15, 33, 38], [45, 47, 53, 51], [3, 12, 30, 41], [46, 50, 52, 48]]
}

CUBIE_IDX = {
    "U": (0, slice(None), slice(None)),
    "F": (slice(None), 2, slice(None)),
    "R": (slice(None), slice(None, None, -1), 2),
    "D": (2, slice(None, None, -1), slice(None)),
    "B": (slice(None), 0, slice(None, None, -1)),
    "L": (slice(None), slice(None), 0),
}

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
