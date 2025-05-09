=====
Usage
=====

You can use the command `cube` strightaway::

    $ cube --help
    $ cube scramble

The first time you solve a cube, it will generate required tables at it will take around::

    $ cube solve -r




To use Cube Solver in a project

.. code-block:: python

    from cube_solver import Cube, Solver

    scramble = Cube.generate_scramble()
    print("Scramble:", scramble)

    cube = Cube(scramble)
    print(cube)

    solver = Solver()
    solution = solver.thistlethwaite(cube)
    print("Solution:", solution)
