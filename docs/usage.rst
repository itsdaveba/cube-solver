=====
Usage
=====

After installation (see :doc:`installation guide <installation>`), you can use the ``cube`` command straight away:

.. code-block:: console

    $ cube --help

To generate a scramble, use the ``scramble`` subcommand:

.. code-block:: console

    $ cube scramble --help
    $ cube scramble              # random scramble of length 25
    $ cube scramble --length 30  # random scramble of specific length

To solve a cube using the **Thistlethwaite** algorithm, use the ``solve`` subcommand.
The first time you solve a cube, it will generate the required tables, which takes around 5 minutes:

.. code-block:: console

    $ cube solve --help
    $ cube solve "U F R"      # solve specific scramble
    $ cube solve --random     # solve random scramble of length 25
    $ cube solve --random 30  # solve random scramble of specific length

To use **Cube Solver** in a Python project:

.. code-block:: python

    from cube_solver import Cube, Solver

    scramble = Cube.generate_scramble()
    print("Scramble:", scramble)

    cube = Cube(scramble)
    print(cube)

    solver = Solver(transition_tables=True, pruning_tables=True)
    solution = solver.thistlethwaite(cube)
    print("Solution:", solution)
