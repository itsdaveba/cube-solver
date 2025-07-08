cube\_solver.cli module
=======================

Cube Solver CLI
---------------

**Usage**:

.. code:: console

   $ cube [OPTIONS] COMMAND [ARGS]...

**Options**:

- ``--help``: Show this message and exit.

**Commands**:

- ``maneuver``: Apply a sequence of moves to a cube.
- ``scramble``: Generate a random scramble.
- ``solve``: Solve a cube.

``cube maneuver``
-----------------

Apply a sequence of moves to a cube.

Accepts the following move types:

- Face moves (e.g. `U`, `F2`, `R'`).
- Slice moves (e.g. `M`, `E2`, `S'`).
- Wide moves (e.g. `Uw`, `Fw2`, `Rw'`` or `u`, `f2`, `r'`).
- Rotations (e.g. `x`, `y2`, `z'`).

**Usage**:

.. code:: console

   $ cube maneuver [OPTIONS] [MOVES]

**Arguments**:

- ``[MOVES]``: Sequence of moves.

**Options**:

- ``-c, --cube TEXT``: Starting cube string representation.
- ``--help``: Show this message and exit.

``cube scramble``
-----------------

Generate a random scramble.

**Usage**:

.. code:: console

   $ cube scramble [OPTIONS]

**Options**:

- ``-l, --length INTEGER``: Scramble length. `[default: 25]`
- ``--wca``: Scramble following WCA rules (uses the Kociemba solver).
- ``-v, --verbose``: Show cube layout and cube string representation. `[default: 0]`
- ``--help``: Show this message and exit.

``cube solve``
--------------

Solve a cube.

The cube string representation must contain characters from `{W, G, R, Y, B, O}`,
representing the colors ``WHITE``, ``GREEN``, ``RED``, ``YELLOW``, ``BLUE``, and ``ORANGE``, respectively.
The face order of the string representation is: ``UP``, ``LEFT``, ``FRONT``, ``RIGHT``, ``BACK``, ``DOWN``.

Example:

The string representation of a cube with ``WHITE`` on the ``UP`` face and ``GREEN``
on the ``FRONT`` face, after the scramble ``R U R' U'``, is: ``WWOWWGWWGBOOOOOOOOGGYGGWGGGRRWBRRWRRBRRBBBBBBYYRYYYYYY``

**Usage**:

.. code:: console

   $ cube solve [OPTIONS] [CUBE]

**Arguments**:

- ``[CUBE]``: Cube string representation.

**Options**:

- ``-a, --algorithm [kociemba|thistle]``: Solver algorithm. `[default: kociemba]`
- ``-l, --length INTEGER``: Maximum solution length.
- ``-t, --timeout INTEGER``: Maximum time in seconds.
- ``-s, --scramble TEXT``: Cube scramble.
- ``-r, --random``: Solve a random cube.
- ``-o, --optimal``: Find the optimal solution.
- ``-v, --verbose``: Show cube layout and cube string representation. `[default: 0]`
- ``--help``: Show this message and exit.
