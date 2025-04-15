"""Console script for cube_solver."""
from cube_solver import Cube

import typer
from typing_extensions import Annotated
from rich.console import Console

app = typer.Typer(no_args_is_help=True, add_completion=False)
console = Console()


@app.command()
def solve(scramble: Annotated[str, typer.Argument()] = None, size: int = 3, scramble_length: int = None):
    """Initialize a new cube."""
    if scramble_length is not None:
        scramble = Cube.generate_scramble(scramble_length)
    console.print("Scramble:", scramble)
    cube = Cube(scramble, size)
    console.print(cube)
    solution = cube.solve()
    console.print("Solution:", solution)


@app.command()
def scramble(length: int = 25):
    """Generate a random scramble."""
    scramble = Cube.generate_scramble(length)
    console.print(scramble)
    cube = Cube(scramble)
    console.print(cube)


if __name__ == "__main__":
    app()
