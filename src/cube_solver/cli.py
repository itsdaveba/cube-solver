"""Console script for cube_solver."""
from cube_solver import Cube

import typer
from typing_extensions import Annotated
from rich.console import Console

app = typer.Typer(no_args_is_help=True, add_completion=False)
console = Console()


@app.command()
def init(scramble: Annotated[str, typer.Argument()] = None, size: int = 3):
    """Console script for cube-solver."""
    cube = Cube(scramble, size)
    console.print(cube)


@app.command()
def scramble(length: int = 25):
    """Generate a random scramble."""
    scramble = Cube.generate_scramble(length)
    console.print(scramble)


if __name__ == "__main__":
    app()
