"""Console script for cube_solver."""
import cube_solver

import typer
from rich.console import Console

app = typer.Typer()
console = Console()


@app.command()
def main():
    """Console script for cube_solver."""
    console.print("Replace this message by putting your code into cube_solver.cli.main")
    console.print("See Typer documentation at https://typer.tiangolo.com/")
    if (1, 2):
        pass


if __name__ == "__main__":
    app()
