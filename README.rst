===========
Cube Solver
===========


.. image:: https://img.shields.io/pypi/v/cube-solver.svg
        :target: https://pypi.python.org/pypi/cube-solver

.. image:: https://readthedocs.org/projects/cube-solver/badge/?version=latest
        :target: https://cube-solver.readthedocs.io/en/latest/?version=latest
        :alt: Documentation Status




Rubik's Cube Solver


* Free software: MIT License
* Documentation: https://cube-solver.readthedocs.io.


Features
--------

* Command line interface
* Transition and pruning tables
* Thistlethwaite solver algorithm


============
Installation
============


Stable release
--------------

To install Cube Solver, run this command in your terminal:

.. code-block:: console

    pip install cube-solver

This is the preferred method to install Cube Solver, as it will always install the most recent stable release.


=====
Usage
=====

After intallation, you can use the ``cube`` command straight away:

.. code-block:: console

    cube --help

.. code-block:: console

    cube scramble

The first time you solve a cube, it will generate the required tables, which takes around 3 minutes:

.. code-block:: console

    cube solve -r

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


=======
Credits
=======

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
