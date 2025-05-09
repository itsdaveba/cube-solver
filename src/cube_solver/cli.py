"""Console script for cube_solver."""
from cube_solver import Cube, Solver

import typer
from rich.console import Console
from typing_extensions import Annotated

app = typer.Typer(no_args_is_help=True, add_completion=False)
console = Console()


@app.command()
def scramble(length: Annotated[int, typer.Option("--length", "-l", help="Scramble length.")] = 25):
    """Generate a random scramble."""
    scramble = Cube.generate_scramble(length)
    console.print(scramble)
    cube = Cube(scramble)
    console.print(cube)


@app.command()
def solve(scramble: Annotated[str, typer.Argument()] = None, random: Annotated[bool, typer.Option("--random", "-r", help="Generate a random scramble.")] = False):
    """Solve a cube."""
    if (scramble is None) != random:
        console.print("You must provide either 'scramble' or the '--random' option.")
        return

    if random:
        scramble = Cube.generate_scramble()
    console.print("Scramble:", scramble)
    cube = Cube(scramble)
    console.print(cube)
    solver = Solver(transition_tables=True, pruning_tables=True)
    solution = solver.thistlethwaite(cube)
    console.print("Solution:", solution)


if __name__ == "__main__":
    app()
