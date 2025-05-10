=====
Usage
=====

After intallation (see :doc:`installation guide <installation>`), you can use the ``cube`` command straight away::

    $ cube --help
    $ cube scramble

The first time you solve a cube, it will generate the required tables, which takes around 3 minutes::

    $ cube solve -r

To use Cube Solver in a Python project:

.. code-block:: python

    from cube_solver import Cube, Solver

    scramble = Cube.generate_scramble()
    print("Scramble:", scramble)

    cube = Cube(scramble)
    print(cube)

    solver = Solver(transition_tables=True, pruning_tables=True)
    solution = solver.thistlethwaite(cube)
    print("Solution:", solution)
