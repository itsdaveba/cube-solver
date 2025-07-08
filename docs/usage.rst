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
            ---------
            | W W O |
            | W W G |
            | W W G |
    ---------------------------------
    | B O O | G G Y | R R W | B R R |
    | O O O | G G W | B R R | B B B |
    | O O O | G G G | W R R | B B B |
    ---------------------------------
            | Y Y R |
            | Y Y Y |
            | Y Y Y |
            ---------
    Cube: WWOWWGWWGBOOOOOOOOGGYGGWGGGRRWBRRWRRBRRBBBBBBYYRYYYYYY

To generate a scramble, use the ``scramble`` subcommand:

.. code-block:: console

    $ cube scramble --help       # scramble subcommand help
    $ cube scramble              # scramble of length 25
    L' U' B' L' D2 R' D B L D R2 F D2 F L2 F2 B' R' F' B U' R L2 D B2

    $ cube scramble --length 30  # scramble of length 30
    F' L2 B' R D' L D2 F2 B L' F B' U' L' U' F2 B2 D' R B R D' L2 D L2 B2 L' U' L F'

    $ cube scramble --wca        # scramble following WCA rules (uses the Kociemba solver)
    R2 U' F2 U' R2 B2 U2 F2 D' B2 U' F2 D2 R2 B' R' D' L' R B D' R' B2

To solve a cube, use the ``solve`` subcommand.
The first time you solve a cube with a specific algorithm,
the required tables will be generated. This process takes around 5 minutes.

The cube string representation must contain characters from `{'W', 'G', 'R', 'Y', 'B', 'O'}`,
representing the colors ``WHITE``, ``GREEN``, ``RED``, ``YELLOW``, ``BLUE``, and ``ORANGE``, respectively.
The order of the string representation is:

.. code-block:: console

               ------------
               | 01 02 03 |
               | 04 05 06 |
               | 07 08 09 |
    ---------------------------------------------
    | 10 11 12 | 19 20 21 | 28 29 30 | 37 38 39 |
    | 13 14 15 | 22 23 24 | 31 32 33 | 40 41 42 |
    | 16 17 18 | 25 26 27 | 34 35 36 | 43 44 45 |
    ---------------------------------------------
               | 46 47 48 |
               | 49 50 51 |
               | 52 53 54 |
               ------------


.. code-block:: console

    $ cube solve --help                                                  # solve subcommand help
    $ cube solve RGWWWWWWRWOOOOOOOOGGGGGWGGWYBBRRRRRRORBBBBBBBYYGYYYYYY  # solve cube representation
    R U R' U'

    $ cube solve --scramble "U R U' R'"                                  # solve scramble
    R U R' U'

    $ cube solve --random --verbose                                      # solve random cube
            ---------
            | W O Y |
            | W W G |
            | O R Y |
    ---------------------------------
    | R O W | B G G | O W B | R G G |
    | O O B | O G G | Y R W | R B Y |
    | B R W | B R O | W W G | R Y Y |
    ---------------------------------
            | R B G |
            | Y Y B |
            | O B Y |
            ---------
    Cube: WOYWWGORYROWOOBBRWBGGOGGBROOWBYRWWWGRGGRBYRYYRBGYYBOBY
    Solution: U2 D2 F' D' L' B' U2 B2 R U F2 D B2 R2 D F2 B2 R2 D F2 L2 U2 L2 U' (24)

    $ cube solve --random --verbose --verbose --algorithm kociemba       # Kociemba algorithm (default)
            ---------
            | W G W |
            | B W G |
            | O R O |
    ---------------------------------
    | O Y G | W B B | Y R G | R W B |
    | G O Y | R G B | O R W | R B Y |
    | G Y Y | R O B | R B G | R O O |
    ---------------------------------
            | B G W |
            | O Y W |
            | Y W Y |
            ---------
    Cube: WGWBWGOROOYGGOYGYYWBBRGBROBYRGORWRBGRWBRBYROOBGWOYWYWY
    Solution: ["D L2 B L B U2 D2 B' L", "U' R2 D' R2 B2 D2 R2 U' F2 D2 R2 D' R2 D' B2"] (24)

    $ cube solve --random --verbose --verbose --algorithm thistle        # Thistlethwaite algorithm
            ---------
            | B Y R |
            | B W B |
            | W B B |
    ---------------------------------
    | O O R | B W O | W R W | G R Y |
    | G O W | R G G | R R W | G B O |
    | G Y Y | G O Y | R O O | W Y O |
    ---------------------------------
            | R Y B |
            | B Y W |
            | Y G G |
            ---------
    Cube: BYRBWBWBBOORGOWGYYBWORGGGOYWRWRRWROOGRYGBOWYORYBBYWYGG
    Solution: ["B' R2 D F", "U R U R' L' U2 D L", "U' R2 U F2 L2 F2 B2 U", 'R2 F2 D2 F2 R2 U2 R2'] (27)

    $ cube solve --scramble "L2 U R D' B2 D2 F B D" --optimal --verbose  # find the optimal solution
            ---------
            | B O Y |
            | W W G |
            | B O G |
    ---------------------------------
    | W G Y | O B R | W R B | R G O |
    | W O Y | O G R | W R B | R B O |
    | G W O | W B R | W Y Y | G G O |
    ---------------------------------
            | G Y B |
            | B Y R |
            | Y Y R |
            ---------
    Cube: BOYWWGBOGWGYWOYGWOOBROGRWBRWRBWRBWYYRGORBOGGOGYBBYRYYR
    INFO: Solution: D' F' B' U2 F2 D L L2 F2 D2 L2 F2 U D L2 B2 D L2 (18)
    INFO: Solution: D' F' B' U2 F2 D L' F2 D2 L2 F2 U D L2 B2 D L2 (17)
    INFO: Solution: D' F' B' D2 B2 D R R2 U' L2 (10)
    INFO: Solution: D' F' B' D2 B2 D R' U' L2 (9)
    Optimal: D' F' B' D2 B2 D R' U' L2 (9)

    $ cube solve --random --optimal --verbose --timeout 10               # stop search after 10 seconds
            ---------
            | B W G |
            | O W Y |
            | R R R |
    ---------------------------------
    | W B G | Y G B | W O W | O B O |
    | W O B | Y G W | R R R | Y B G |
    | Y Y W | R W O | G R Y | R G O |
    ---------------------------------
            | G O Y |
            | G Y B |
            | B O B |
            ---------
    Cube: BWGOWYRRRWBGWOBYYWYGBYGWRWOWOWRRRGRYOBOYBGRGOGOYGYBBOB
    INFO: Solution: U2 R' B U2 L' F' U F2 L R2 U L2 U F2 B2 U2 L2 D L2 B2 D (21)
    INFO: Solution: U B' L2 F D' R L F' R U' B' U' B2 D2 R2 B2 D F2 U2 F2 (20)
    Suboptimal: U B' L2 F D' R L F' R U' B' U' B2 D2 R2 B2 D F2 U2 F2 (20)

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
    assert solution == scramble.inverse
    print(f"Solution: {solution} ({len(solution)})")
