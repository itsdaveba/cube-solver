import numpy
from pathlib import Path
from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="cube_solver.csolver",
            sources=[*Path("src/cmodules").iterdir()],
            include_dirs=["src/include/", numpy.get_include()]
        ),
    ]
)
