"""Console script for cube_solver."""
from cube_solver import Cube

import typer
from typing_extensions import Annotated
from rich.console import Console

app = typer.Typer()
console = Console()


@app.command()
def main(scramble: Annotated[str, typer.Argument()] = None, size: int = 3):
    """Console script for cube-solver."""
    cube = Cube(scramble, size)
    console.print(cube)


if __name__ == "__main__":
    app()
