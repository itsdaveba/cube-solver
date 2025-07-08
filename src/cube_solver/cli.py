"""Console script for cube_solver."""
from enum import Enum
from cube_solver import Cube, Maneuver, BaseSolver, Thistlethwaite, Kociemba

import typer
from typing import Union
from rich.console import Console
from typing_extensions import Annotated


app = typer.Typer(add_completion=False)
console = Console()


class Algorithm(str, Enum):
    KOCIEMBA = "kociemba"
    THISTLETHWAITE = "thistle"


ALGS = {
    Algorithm.KOCIEMBA: Kociemba,
    Algorithm.THISTLETHWAITE: Thistlethwaite
}


@app.command()
def maneuver(moves: Annotated[str, typer.Argument(help="Sequence of moves.")] = "",
             cube: Annotated[str, typer.Option("--cube", "-c", help="Starting cube string representation.")] = ""):
    """
    Apply a sequence of moves to a cube.

    Accepts the following move types:
    \n\n* Face moves (e.g. U, F2, R').
    \n\n* Slice moves (e.g. M, E2, S').
    \n\n* Wide moves (e.g. Uw, Fw2, Rw' or u, f2, r').
    \n\n* Rotations (e.g. x, y2, z').
    """
    _cube = Cube(repr=cube) if cube else Cube()
    _cube.apply_maneuver(moves)
    console.print(_cube)
    console.print(f"Cube: {repr(_cube)}")


@app.command()
def scramble(length: Annotated[int, typer.Option("--length", "-l", show_envvar=False, help="Scramble length.")] = 25,
             wca: Annotated[bool, typer.Option("--wca",
                                               help="Scramble following WCA rules (uses the Kociemba solver).")] = False,
             verbose: Annotated[int, typer.Option("--verbose", "-v", count=True,
                                                  help="Show cube layout and cube string representation.")] = 0):
    """Generate a random scramble."""
    if wca:
        solver = Kociemba()
        while True:
            cube = Cube(random_state=True)
            solution = solver.solve(cube, length)
            assert isinstance(solution, Maneuver)
            scramble = solution.inverse
            if solver.solve(cube, 1) is None:  # pragma: no cover
                break
    else:
        scramble = Maneuver.random(length)
    if verbose:
        console.print(f"Scramble: {scramble}")
        cube = Cube(scramble)
        console.print(cube)
        console.print(f"Cube: {repr(cube)}")
    else:
        console.print(f"{scramble}")


@app.command()
def solve(cube: Annotated[str, typer.Argument(help="Cube string representation.")] = "",
          algorithm: Annotated[Algorithm, typer.Option("--algorithm", "-a", show_envvar=False,
                                                       help="Solver algorithm.", show_choices=True)] = Algorithm.KOCIEMBA,
          length: Annotated[Union[int, None], typer.Option("--length", "-l", show_envvar=False,
                                                           help="Maximum solution length.")] = None,
          timeout: Annotated[Union[int, None], typer.Option("--timeout", "-t", show_envvar=False,
                                                            help="Maximum time in seconds.")] = None,
          scramble: Annotated[str, typer.Option("--scramble", "-s", show_envvar=False, help="Cube scramble.")] = "",
          random: Annotated[bool, typer.Option("--random", "-r", help="Solve a random cube.")] = False,
          optimal: Annotated[bool, typer.Option("--optimal", "-o", help="Find the optimal solution.")] = False,
          verbose: Annotated[int, typer.Option("--verbose", "-v", count=True,
                                               help="Show cube layout and cube string representation.")] = 0):
    """
    Solve a cube.

    The cube string representation must contain characters from {W, G, R, Y, B, O},
    \n\nrepresenting the colors WHITE, GREEN, RED, YELLOW, BLUE, and ORANGE, respectively.
    \n\nThe face order of the string representation is: UP, LEFT, FRONT, RIGHT, BACK, DOWN.

    \n\nExample:
    \n\nThe string representation of a cube with WHITE on the UP face and GREEN on the FRONT face,
    \n\nafter the scramble R U R' U', is:
    \n\nWWOWWGWWGBOOOOOOOOGGYGGWGGGRRWBRRWRRBRRBBBBBBYYRYYYYYY
    """
    solver: BaseSolver = ALGS[algorithm]()
    if not cube and not scramble and not random:
        msg = "Must provide either the 'cube' argument, the '--scramble' / '-s' option, or the '--random' / '-r' option."
        console.print(msg)
        raise typer.Exit(1)
    if cube and scramble:
        console.print("The '--scramble' / '-s' option cannot be used with the 'cube' argument.")
        raise typer.Exit(1)
    if cube and random:
        console.print("The '--random' / '-r' option cannot be used with the 'cube' argument.")
        raise typer.Exit(1)
    if scramble and random:
        console.print("The '--random' / '-r' option cannot be used with the '--scramble' / '-s' option.")
        raise typer.Exit(1)
    _cube = Cube(repr=cube) if cube else Cube(scramble) if scramble else Cube(random_state=True)
    if verbose:
        console.print(_cube)
        console.print(f"Cube: {repr(_cube)}")
        solution = solver.solve(_cube, length, optimal, timeout, verbose)
        if solution is not None:
            length = len(solution) if isinstance(solution, Maneuver) else sum(len(sol) for sol in solution)
            if optimal:
                if solver.terminated:
                    console.print(f"Suboptimal: {solution} ({length})")
                else:
                    console.print(f"Optimal: {solution} ({length})")
            else:
                console.print(f"Solution: {solution} ({length})")
        else:
            console.print(f"Solution: {solution}")
    else:
        solution = solver.solve(_cube, length, optimal, timeout)
        console.print(f"{solution}")


if __name__ == "__main__":
    app()
