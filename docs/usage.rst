=====
Usage
=====

After installation (see :doc:`installation guide <installation>`), you can use the ``cube`` command straight away:

.. code-block:: console

    $ cube --help

To generate a scramble, use the ``scramble`` subcommand:

.. code-block:: console

    $ cube scramble

To solve a cube using the **Thistlethwaite** algorithm, use the ``solve`` subcommand.
The first time you solve a cube, it will generate the required tables, which takes around 5 minutes:

.. code-block:: console

    $ cube solve -r

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
