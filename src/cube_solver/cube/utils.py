import math
import numpy as np

from .defs import FACTORIAL, COMBINATION


# def get_orientation_coord(orientation: np.ndarray, v: int) -> int:
#     """
#     Get orientation coordinate.

#     Given the length of the `orientaion` array is `n`,
#     the number of possible orientation values is `v`,
#     and the array contains values between `0` and `v - 1`,
#     returns a unique number between `0` and `n ^ v - 1`
#     for every possible orientation.

#     Parameters
#     ----------
#     orientation : ndarray
#         Array containing orientation values between `0` and `v - 1`.
#     v : int
#         Number of possible orientation values.

#     Returns
#     -------
#     coord : int
#         Orientation coordinate.

#     Examples
#     --------
#     >>> from cube_solver.cube.utils import get_orientation_coord
#     >>> import numpy as np
#     >>> orientation = np.array([2, 1, 0])
#     >>> get_orientation_coord(orientation, 3)
#     21
#     """
#     if not isinstance(orientation, np.ndarray):
#         raise TypeError(f"orientation must be ndarray, not {type(orientation).__name__}")
#     if not isinstance(v, int):
#         raise TypeError(f"v must be int, not {type(v).__name__}")

#     coord = np.array([0])[0]
#     for o in orientation:
#         coord *= v
#         coord += o
#     return coord.item()


# def get_permutation_coord(permutation: np.ndarray, parity: bool = False) -> int:
#     """
#     Get permutation coordinate.

#     Given the length of the `permutation` array is `n`,
#     and the array contains all values between `0` and `n - 1`,
#     returns a unique number between `0` and `n! - 1`
#     for every posible permutation.

#     Parameters
#     ----------
#     permutation : ndarray
#         Array containing all values between `0` and `n - 1`.

#     Returns
#     -------
#     coord : int
#         Permutation coordinate.

#     Examples
#     --------
#     >>> from cube_solver.cube.utils import get_permutation_coord
#     >>> import numpy as np
#     >>> permutation = np.array([2, 1, 0])
#     >>> get_permutation_coord(permutation)
#     5
#     """
#     if not isinstance(permutation, np.ndarray):
#         raise TypeError(f"permutation must be ndarray, not {type(permutation).__name__}")

#     coord = np.array([0])[0]
#     num_elems = len(permutation)
#     for i in range(num_elems - (1 if not parity else 2)):
#         coord *= num_elems - i
#         coord += np.sum(permutation[i] > permutation[i+1:])
#     return coord.item()


# def get_combination_coord(combination: np.ndarray) -> int:
#     """
#     Get combination coordinate.

#     Given the length of the `combination` array is `n`,
#     the maximum value of the array is `m - 1`,
#     and the array contains different values in increasing order between `0` and `m - 1`,
#     returns a unique number between `0` and `C(m , n) - 1`
#     for every possible combination.

#     Parameters
#     ----------
#     combination : ndarray
#         Array containing different values in increasing order between `0` and `m`.

#     Returns
#     -------
#     coord : int
#         Combination coordinate.

#     Examples
#     --------
#     >>> from cube_solver.cube.utils import get_combination_coord
#     >>> import numpy as np
#     >>> combination = np.array([3, 4, 5])
#     >>> get_combination_coord(combination)
#     19
#     """
#     if not isinstance(combination, np.ndarray):
#         raise TypeError(f"combination must be ndarray, not {type(combination).__name__}")

#     try:
#         return np.sum(COMBINATION[range(len(combination)), combination]).item()
#     except IndexError:
#         return np.sum([math.comb(c, i + 1) for i, c in enumerate(combination)]).item()


# def get_partial_permutation_coord(permutation: np.ndarray, combination: np.ndarray) -> int:
#     """
#     Get partial permutation coordinate.

#     Given the length of the `permutation` array is `m`,
#     the length of the `combination` array is `n` with `n <= m`,
#     the `permutation` array contains all values between `0` and `m - 1`,
#     and the `combination` array contains different values in increasing order between `0` and `m - 1`,
#     returns a unique number between `0` and `P(m, n) - 1`
#     for every possible partial permutation of the `permutation` array elements
#     defined by the `combination` array indexes.
#     If the `combination` array es empty, returns `-1`.

#     Parameters
#     ----------
#     permutation : ndarray
#         Array containing all values between `0` and `m - 1`.
#     combination : ndarray
#         Array containing different values in increasing order between `0` and `m - 1`.

#     Returns
#     -------
#     coord : int
#         Partial permutation coordinate.

#     Examples
#     --------
#     >>> from cube_solver.cube.utils import get_partial_permutation_coord
#     >>> import numpy as np
#     >>> permutation = np.array([3, 4, 5, 2, 1, 0])
#     >>> combination = np.array([3, 4, 5])
#     >>> get_partial_permutation_coord(permutation, combination)
#     119
#     """
#     if not isinstance(permutation, np.ndarray):
#         raise TypeError(f"permutation must be ndarray, not {type(permutation).__name__}")
#     if not isinstance(combination, np.ndarray):
#         raise TypeError(f"combination must be ndarray, not {type(combination).__name__}")
#     if not combination.size:
#         return -1

#     coord = get_permutation_coord(permutation[combination])
#     try:
#         return (coord + FACTORIAL[len(combination)-1] * get_combination_coord(combination)).item()
#     except IndexError:
#         return coord + math.factorial(len(combination)) * get_combination_coord(combination)


def get_orientation_array(coord: int, v: int, n: int, force_modulo: bool = False) -> np.ndarray:
    """
    Get orientation array.

    Given the `coord` number between `0` and `v ^ n - 1`
    where `v` is the number of possible orientation values,
    returns the `orientation` array of length `n`
    with orientation values between `0` and `v - 1`.

    If `force_modulo` is `True`, given the `coord` number
    between `0` and `v ^ (n - 1) - 1`, returns the
    `orientation` array with total orientation `0` modulo `v`.

    Parameters
    ----------
    coord : int
        Orientation coordinate value between `0` and `v ^ n - 1`.
        If `force_modulo` is `True`, orientation coordinate value
        between `0` and `v ^ (n - 1) - 1`.
    v : int
        Number of possible orientation values.
    n : int
        Length of orientation array.
    force_modulo : bool, optional
        If `True`, returns the `orientation` array
        with total orientation `0` modulo `v`.

    Returns
    -------
    array : ndarray
        Orientation array.

    Examples
    --------
    >>> from cube_solver.cube.utils import get_orientation_array
    >>> get_orientation_array(6, 3, 3)
    array([0, 2, 0])
    >>> get_orientation_array(6, 3, 3, force_modulo=True)
    array([2, 0, 1])
    """
    if not isinstance(coord, int):
        raise TypeError(f"coord must be int, not {type(coord).__name__}")
    if not isinstance(v, int):
        raise TypeError(f"v must be int, not {type(v).__name__}")
    if not isinstance(n, int):
        raise TypeError(f"n must be int, not {type(n).__name__}")
    if not isinstance(force_modulo, bool):  # TODO add value errors with lims
        raise TypeError(f"force_modulo must be bool, not {type(force_modulo).__name__}")
    if v <= 0:
        raise ValueError(f"v must be positive (got {v})")
    if n <= 0:
        raise ValueError(f"n must be positive (got {n})")
    upper_lim = v ** n
    if force_modulo:
        upper_lim //= v
    if coord < 0 or coord >= upper_lim:
        raise ValueError(f"coord must be >= 0 and < {upper_lim} (got {coord})")

    orientation = np.zeros(n, dtype=int)
    for i in range(len(orientation) - (2 if force_modulo else 1), -1, -1):
        coord, orientation[i] = divmod(coord, v)
    if force_modulo:
        orientation[-1] = -np.sum(orientation[:-1]) % v
    return orientation


def get_permutation_array(coord: int, n: int, force_even_parity: bool = False) -> tuple[np.ndarray, bool]:
    """
    Get permutation array and permutation parity.

    Given the `coord` number between `0` and `n! - 1`,
    returns the `permutation` array of length `n`
    with permutation values between `0` and `n - 1`.
    The `permutation parity` is `True` if the `permutation` array
    has `odd` parity, `False` if it has `even` parity.

    If `force_even_parity` is `True`, given the `coord` number
    between `0` and `n! / 2 - 1`, returns the
    `permutation` array with `even` parity.

    Parameters
    ----------
    coord : int
        Permutation coordinate value between `0` and `n! - 1`.
        If `force_even_parity` is `True`, permutation coordinate value
        between `0` and `n! / 2 - 1`.
    n : int
        Length of permutation array.
    force_even_parity : bool, optional
        If `True`, returns the `permutation` array with `even` parity.

    Returns
    -------
    tuple of (ndarray, bool)
        Permutation array and permutation parity.

    Examples
    --------
    >>> from cube_solver.cube.utils import get_permutation_array
    >>> get_permutation_array(2, 3)
    (array([1, 0, 2]), True)
    >>> get_permutation_array(2, 3, force_even_parity=True)
    (array([2, 0, 1]), False)
    """
    if not isinstance(coord, int):
        raise TypeError(f"coord must be int, not {type(coord).__name__}")
    if not isinstance(n, int):
        raise TypeError(f"n must be int, not {type(n).__name__}")
    if not isinstance(force_even_parity, bool):
        raise TypeError(f"force_even_parity must be bool, not {type(force_even_parity).__name__}")
    if not force_even_parity and n <= 0:
        raise ValueError(f"n must be positive (got {n})")
    if force_even_parity and n <= 1:
        raise ValueError(f"n must be > 1 (got {n})")
    try:
        upper_lim = FACTORIAL[n].item()
    except IndexError:
        upper_lim = math.factorial(n)
    if force_even_parity:
        upper_lim //= 2
    if coord < 0 or coord >= upper_lim:
        raise ValueError(f"coord must be >= 0 and < {upper_lim} (got {coord})")

    permutation = np.zeros(n, dtype=int)
    if force_even_parity:
        permutation[-1] = 1
    permutation_parity = 0
    for i in range(n - (3 if force_even_parity else 2), -1, -1):
        coord, permutation[i] = divmod(coord, n - i)
        permutation[i+1:] += permutation[i+1:] >= permutation[i]
        permutation_parity += permutation[i]
    if force_even_parity and bool(permutation_parity % 2):
        permutation[-2:] = permutation[[-1, -2]]
    permutation_parity = False if force_even_parity else bool(permutation_parity % 2)
    return permutation, permutation_parity


# def get_permutation_parity(permutation: np.ndarray) -> bool:
#     """
#     Get permutation parity.

#     Given the `permutation` array, returns `True` if it has `odd`
#     parity, or `False` if it has `even` parity.

#     Parameters
#     ----------
#     permutation : ndarray
#         Permutation array to check.

#     Returns
#     -------
#     parity : bool
#         Permutation parity.

#     Examples
#     --------
#     >>> from cube_solver.cube.utils import get_permutation_parity
#     >>> get_permutation_parity([2, 1, 0])
#     True
#     >>> get_permutation_parity([2, 0, 1])
#     False
#     """
#     n = len(permutation)
#     permutation_parity = np.zeros_like(permutation, dtype=int)
#     coord = get_permutation_coord(permutation)
#     for i in range(n - 2, -1, -1):
#         coord, permutation_parity[i] = divmod(coord, n - i)
#     return bool(np.sum(permutation_parity) % 2)


def get_combination_array(coord: int, n: int) -> np.ndarray:
    """
    Get combination array.

    Given the `coord` number between `0` and `C(m, n) - 1`
    where `m - 1` is the maximum value of the `combination` array,
    returns the `combination` array of length `n`
    with combination values between `0` and `m - 1`.

    Parameters
    ----------
    coord : int
        Combination coordinate value between `0` and `C(m, n) - 1`
        where `m - 1` is the maximum value of the `combination` array.
    n : int
        Length of combination array.

    Returns
    -------
    array : ndarray
        Combination array.

    Examples
    --------
    >>> from cube_solver.cube.utils import get_combination_array
    >>> get_combination_array(6, 3)
    array([1, 2, 4])
    """
    if not isinstance(coord, int):
        raise TypeError(f"coord must be int, not {type(coord).__name__}")
    if not isinstance(n, int):
        raise TypeError(f"n must be int, not {type(n).__name__}")
    if n <= 0:
        raise ValueError(f"n must be positive (got {n})")
    if coord < 0:
        raise ValueError(f"coord must be >= 0 (got {coord})")

    if n >= COMBINATION.shape[1] or coord >= COMBINATION[-1, n]:
        def comb(n: int, k: int) -> int: return math.comb(n, k)
    else:
        def comb(n: int, k: int) -> int: return COMBINATION[n, k].item()

    i = n - 1
    m = 0
    while coord >= comb(m, n):
        m += 1
    combination = np.zeros(n, dtype=int)
    for c in range(m - 1, 0, -1):
        if coord >= comb(c, i + 1):
            coord -= comb(c, i + 1)
            combination[i] = c
            i -= 1
            if i < 0:
                break

    return combination


def get_partial_permutation_array(coord: int, n: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Get partial permutation array and combination array.

    Given the `coord` number between `0` and `P(m, n) - 1`
    where `m - 1` is the maximum value of the `combination` array with `n <= m`,
    returns the `partial permutation` array of length `n`
    with partial permutation values between `0` and `n - 1`,
    and the `combination` array of length `n`
    with combination values between `0` and `m - 1`
    defining the partial permutation indexes.

    Parameters
    ----------
    coord : int
        Partial permutation coordinate value between `0` and `P(m, n) - 1`
        where `m - 1` is the maximum value of the `combination` array.
    n : int
        Length of partial permutation array and combination array.

    Returns
    -------
    arrays : tuple of ndarray
        Partial permutation array and combination array.

    Examples
    --------
    >>> from cube_solver.cube.utils import get_partial_permutation_array
    >>> get_partial_permutation_array(38, 3)
    (array([1, 0, 2]), array([1, 2, 4]))
    """
    if not isinstance(coord, int):
        raise TypeError(f"coord must be int, not {type(coord).__name__}")
    if not isinstance(n, int):
        raise TypeError(f"n must be int, not {type(n).__name__}")
    if n <= 0:
        raise ValueError(f"n must be positive (got {n})")
    if coord < 0:
        raise ValueError(f"coord must be >= 0 (got {coord})")

    try:
        comb_coord, perm_coord = divmod(coord, FACTORIAL[n].item())
    except IndexError:
        comb_coord, perm_coord = divmod(coord, math.factorial(n))
    combination = get_combination_array(comb_coord, n)
    permutation = get_permutation_array(perm_coord, n)[0]
    return permutation, combination
