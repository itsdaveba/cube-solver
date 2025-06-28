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

To generate a scramble, use the ``scramble`` subcommand:

.. code-block:: console

    $ cube scramble --help       # scramble subcommand help
    $ cube scramble              # scramble of length 25
    $ cube scramble --length 30  # scramble of length 30
    $ cube scramble --wca        # scramble following WCA rules

To solve a cube, use the ``solve`` subcommand.
The first time you solve a cube with a specific algorithm,
the required tables will be generated. This process takes around 5 minutes.

The cube string representation must contain characters from `{'W', 'G', 'R', 'Y', 'B', 'O'}`,
representing the colors ``WHITE``, ``GREEN``, ``RED``, ``YELLOW``, ``BLUE``, and ``ORANGE``, respectively.
The order of the string representation is::

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
    $ cube solve --scramble "U R U' R'"                                  # solve scramble
    $ cube solve --random                                                # solve random cube
    $ cube solve --random --algorithm thistle                            # solver algorithm
    $ cube solve --scramble "L2 U R D' B2 L' D2 F B D" --optimal         # find the optimal solution
    $ cube solve --random --optimal --verbose --timeout 10               # stop search after 10 seconds

To use **Cube Solver** in a Python project:

.. code-block:: python

    from cube_solver import Cube, Kociemba

    cube = Cube(random_state=True)
    print(cube)
    print("Cube:", repr(cube))

    solver = Kociemba()
    solution = solver.solve(cube)
    print("Solution:", solution)
