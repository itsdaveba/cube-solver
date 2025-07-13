=====
Usage
=====

After installation (see :doc:`installation guide <installation>`), you can use the ``cube`` command straight away:

.. code-block:: console

    $ cube --help

To perform a maneuver to a cube, use the ``maneuver`` subcommand:

.. code-block:: console

    $ cube maneuver --help       # maneuver subcommand help
    $ cube maneuver "R U R' U'"  # apply maneuver
          -------
          | W O |
          | W G |
    -------------------------
    | B O | G Y | R W | B R |
    | O O | G G | W R | B B |
    -------------------------
          | Y R |
          | Y Y |
          -------
    Cube: WOWGBOOOGYGGRWWRBRBBYRYY

To generate a scramble, use the ``scramble`` subcommand:

.. code-block:: console

    $ cube scramble --help       # scramble subcommand help
    $ cube scramble              # scramble of length 15
    B' R U2 L2 U B' R2 B D R U' F L F' R'

    $ cube scramble --length 20  # scramble of length 20
    B2 D' F' D R F2 U' F2 D2 R' B' L D' L2 B' D F' D2 F2 U2

    $ cube scramble --wca        # scramble following WCA rules (uses the Kociemba solver)
    U' R2 U F2 U' F2 U F2 U R2 U' R' F' U2 F2

To solve a cube, use the ``solve`` subcommand.
The first time you solve a cube with a specific algorithm,
the required tables will be generated. This process takes around 1 minute.

The cube string representation must contain characters from `{'W', 'G', 'R', 'Y', 'B', 'O'}`,
representing the colors ``WHITE``, ``GREEN``, ``RED``, ``YELLOW``, ``BLUE``, and ``ORANGE``, respectively.
The order of the string representation is:

.. code-block:: console

            =========
            | 01 02 |
            ---------
            | 03 04 |
    =================================
    | 05 06 | 09 10 | 13 14 | 17 18 |
    ---------------------------------
    | 07 08 | 11 12 | 15 16 | 19 20 |
    =================================
            | 21 22 |
            ---------
            | 24 23 |
            =========


.. code-block:: console

    $ cube solve --help                                             # solve subcommand help
    $ cube solve RWWRWOOOGGGWYBRROBBBYGYY                           # solve cube representation
    R U R' U'

    $ cube solve --scramble "U R U' R'"                             # solve scramble
    R U R' U'

    $ cube solve --random --verbose                                 # solve random cube
          -------
          | B W |
          | R W |
    -------------------------
    | W Y | B G | R R | B O |
    | Y G | O G | R Y | B G |
    -------------------------
          | W Y |
          | O O |
          -------
    Cube: BWRWWYYGBGOGRRRYBOBGWYOO
    Solution: U F' R' F R' F2 R U' (8)

    $ cube solve --random --verbose --verbose --algorithm optimal   # Optimal algorithm (default)
          -------
          | B R |
          | G O |
    -------------------------
    | W R | W B | Y Y | B O |
    | W W | G G | R Y | O R |
    -------------------------
          | O Y |
          | B G |
          -------
    Cube: BRGOWRWWWBGGYYRYBOOROYBG
    Solution: ["U2 F2 U' F R' F' R'"] (7)

    $ cube solve --random --verbose --verbose --algorithm kociemba  # Kociemba algorithm
          -------
          | R B |
          | W B |
    -------------------------
    | G R | B R | Y O | W W |
    | W O | Y G | R Y | O O |
    -------------------------
          | B Y |
          | G G |
          -------
    Cube: RBWBGRWOBRYGYORYWWOOBYGG
    Solution: ["U2 F' U F", "U F2 U' R2 U R2 U' R2 U F2"] (14)

    $ cube solve --random --verbose --algorithm kociemba --optimal  # find the optimal solution
          -------
          | O W |
          | O Y |
    -------------------------
    | Y Y | G B | R R | B B |
    | W R | G G | R W | O O |
    -------------------------
          | W Y |
          | G B |
          -------
    Cube: OWOYYYWRGBGGRRRWBBOOWYGB
    INFO: Solution: L B2 L U2 L L2 U L2 U2 B2 U' L2 (12)
    INFO: Solution: L B2 L U2 L' U L2 U2 B2 U' L2 (11)
    INFO: Solution: B U L U2 L2 B2 L B L2 U2 (10)
    INFO: Solution: B2 U L' B' U B' U B (8)
    Optimal: F2 U R' F' U F' U F (8)

To use **Cube Solver** in a Python project:

.. code-block:: python

    from cube_solver import Cube, Maneuver, Kociemba

    scramble = Maneuver.random()
    print(f"Scramble: {scramble}")

    cube = Cube(scramble)
    print(cube)
    print(f"Cube: {repr(cube)}")

    solver = Kociemba()
    solution = solver.solve(cube)
    assert solution is not None
    print(f"Solution: {solution} ({len(solution)})")
