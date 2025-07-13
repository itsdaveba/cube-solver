"""Cube utils module."""
import math
import numpy as np
from typing import Tuple

from .defs import FACTORIAL


def get_orientation_array(coord: int, v: int, n: int, force_modulo: bool = False) -> np.ndarray:
    """
    Get orientation array.

    Given an orientation coordinate number ``coord`` with a value between
    ``0`` and ``v ^ n - 1``, where ``v`` is the number of possible orientation values,
    returns a unique orientation array of length ``n`` with values between
    ``0`` and ``v - 1`` for every possible orientation coordinate value.

    If ``force_modulo`` is ``True``, given an orientation coordinate number with a value between
    ``0`` and ``v ^ (n - 1) - 1``, returns a unique orientation array with total orientation ``0 modulo v``.

    Parameters
    ----------
    coord : int
        Orientation coordinate number with a value between ``0`` and ``v ^ n - 1``,
        or between ``0`` and ``v ^ (n - 1) - 1`` if ``force_modulo`` is ``True``.
    v : int
        Number of possible orientation values.
    n : int
        Length of the orientation array.
    force_modulo : bool, optional
        If ``True``, returns a unique orientation array
        with total orientation ``0 modulo v``. Default is ``False``.

    Returns
    -------
    orientation : ndarray
        Orientation array of length ``n`` with values between ``0`` and ``v - 1``.

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
    if not isinstance(force_modulo, bool):
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
    for i in range(n - (2 if force_modulo else 1), -1, -1):
        coord, orientation[i] = divmod(coord, v)
    if force_modulo:
        orientation[-1] = -np.sum(orientation[:-1]) % v
    return orientation


def get_orientation_coord(orientation: np.ndarray, v: int, is_modulo: bool = False) -> int:
    """
    Get orientation coordinate number.

    Given an ``orientation`` array of length ``n`` with values between
    ``0`` and ``v - 1``, where ``v`` is the number of possible orientation values,
    returns a unique orientation coordinate number with a value between
    ``0`` and ``v ^ n - 1`` for every possible orientation array.

    If ``is_modulo`` is ``True``, given an orientation array with total orientation ``0 modulo v``,
    returns a unique orientation coordinate number with a value between ``0`` and ``v ^ (n - 1) - 1``.

    Parameters
    ----------
    orientation : ndarray
        Orientation array of length ``n`` with values between ``0`` and ``v - 1``.
    v : int
        Number of possible orientation values.
    is_modulo : bool, optional
        If ``True``, returns a unique orientation coordinate number
        for every possible orientation array with total orientation ``0 modulo v``.
        Default is ``False``.

    Returns
    -------
    coord : int
        Orientation coordinate number with a value between ``0`` and ``v ^ n - 1``,
        or between ``0`` and ``v ^ (n - 1) - 1`` if ``is_modulo`` is ``True``.

    Examples
    --------
    >>> import numpy as np
    >>> from cube_solver.cube.utils import get_orientation_coord
    >>> get_orientation_coord(np.array([0, 2, 0]), 3)
    6
    >>> get_orientation_coord(np.array([2, 0, 1]), 3, is_modulo=True)
    6
    """
    if not isinstance(orientation, np.ndarray):
        raise TypeError(f"orientation must be ndarray, not {type(orientation).__name__}")
    if not isinstance(v, int):
        raise TypeError(f"v must be int, not {type(v).__name__}")
    if not isinstance(is_modulo, bool):
        raise TypeError(f"is_modulo must be bool, not {type(is_modulo).__name__}")
    if v <= 0:
        raise ValueError(f"v must be positive (got {v})")
    n = len(orientation)
    if n <= 0:
        raise ValueError(f"orientation length must be positive (got {n})")
    if not np.issubdtype(orientation.dtype, np.integer):
        raise TypeError(f"orientation elements must be int, not {orientation.dtype}")
    if np.any((orientation < 0) | (orientation >= v)):
        raise ValueError(f"orientation values must be >= 0 and < {v} (got {orientation})")
    if is_modulo and np.sum(orientation) % v != 0:
        raise ValueError(f"orientation has no total orientation 0 modulo {v} (got {np.sum(orientation) % v})")

    coord = np.array([0])[0]
    for o in orientation[slice(n - 1 if is_modulo else n)]:
        coord *= v
        coord += o
    return coord.item()


def get_permutation_array(coord: int, n: int) -> Tuple[np.ndarray, bool]:
    """
    Get permutation array and permutation parity.

    Given a permutation coordinate number ``coord`` with a value between
    ``0`` and ``n! - 1``, returns a unique permutation array of length ``n``
    with values between ``0`` and ``n - 1`` for every possible permutation coordinate value.
    Also returns the permutation parity: ``True`` if the permutation array
    has ``odd`` parity, ``False`` if it has ``even`` parity.

    Parameters
    ----------
    coord : int
        Permutation coordinate number with a value between ``0`` and ``n! - 1``,
        or between ``0`` and ``n! / 2 - 1`` if ``force_even_parity`` is ``True``.
    n : int
        Length of the permutation array.

    Returns
    -------
    permutation : ndarray
        Permutation array of length ``n`` with values between ``0`` and ``n - 1``.
    parity : bool
        Permutation parity: ``True`` if the permutation array has ``odd`` parity,
        ``False`` if it has ``even`` parity.

    Examples
    --------
    >>> from cube_solver.cube.utils import get_permutation_array
    >>> get_permutation_array(2, 3)
    (array([1, 0, 2]), True)
    """
    if not isinstance(coord, int):
        raise TypeError(f"coord must be int, not {type(coord).__name__}")
    if not isinstance(n, int):
        raise TypeError(f"n must be int, not {type(n).__name__}")
    if n <= 0:
        raise ValueError(f"n must be positive (got {n})")
    try:
        upper_lim = FACTORIAL[n].item()
    except IndexError:
        upper_lim = math.factorial(n)
    if coord < 0 or coord >= upper_lim:
        raise ValueError(f"coord must be >= 0 and < {upper_lim} (got {coord})")

    permutation = np.zeros(n, dtype=int)
    permutation_parity = 0
    for i in range(n - 2, -1, -1):
        coord, permutation[i] = divmod(coord, n - i)
        permutation[i+1:] += permutation[i+1:] >= permutation[i]
        permutation_parity += permutation[i]
    permutation_parity = bool(permutation_parity % 2)
    return permutation, permutation_parity


def get_permutation_coord(permutation: np.ndarray) -> int:
    """
    Get permutation coordinate number.

    Given a ``permutation`` array of length ``n`` with ``n`` different values,
    returns a unique permutation coordinate number with a value between
    ``0`` and ``n! - 1`` for every possible permutation array.

    Parameters
    ----------
    permutation : ndarray
        Permutation array of length ``n`` with ``n`` different values.

    Returns
    -------
    coord : int
        Permutation coordinate number with a value between ``0`` and ``n! - 1``,
        or between ``0`` and ``n! / 2 - 1`` if ``is_even_parity`` is ``True``.

    Examples
    --------
    >>> import numpy as np
    >>> from cube_solver.cube.utils import get_permutation_coord
    >>> get_permutation_coord(np.array([1, 0, 2]))
    2
    """
    if not isinstance(permutation, np.ndarray):
        raise TypeError(f"permutation must be ndarray, not {type(permutation).__name__}")
    n = len(permutation)
    if n <= 0:
        raise ValueError(f"permutation length must be positive (got {n})")
    if not np.issubdtype(permutation.dtype, np.integer):
        raise TypeError(f"permutation elements must be int, not {permutation.dtype}")
    if len(set(permutation)) != n:
        raise ValueError(f"permutation values must be different (got {permutation})")

    coord = np.array([0])[0]
    for i in range(n - 1):
        coord *= n - i
        coord += np.sum(permutation[i] > permutation[i+1:])
    return coord.item()


def get_permutation_parity(permutation: np.ndarray) -> bool:
    """
    Get permutation parity.

    Given a ``permutation`` array of length ``n`` with ``n`` different values,
    returns the permutation parity: ``True`` if the permutation array
    has ``odd`` parity, ``False`` if it has ``even`` parity.

    Parameters
    ----------
    permutation : ndarray
        Permutation array of length ``n`` with ``n`` different values.

    Returns
    -------
    parity : bool
        Permutation parity: ``True`` if the permutation array
        has ``odd`` parity, ``False`` if it has ``even`` parity.

    Examples
    --------
    >>> import numpy as np
    >>> from cube_solver.cube.utils import get_permutation_parity
    >>> get_permutation_parity(np.array([1, 0, 2]))
    True
    >>> get_permutation_parity(np.array([2, 0, 1]))
    False
    """
    coord = get_permutation_coord(permutation)
    _, permutation_parity = get_permutation_array(coord, len(permutation))
    return permutation_parity
