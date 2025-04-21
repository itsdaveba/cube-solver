"""Console script for cube_solver."""
from cube_solver import Cube, Solver

import typer
from rich.console import Console

app = typer.Typer(no_args_is_help=True, add_completion=False)
console = Console()


@app.command()
def solve(scramble: str, size: int = 3, representation: str = "face"):
    """Solve a cube."""
    console.print("Scramble:", scramble)
    cube = Cube(scramble, size, representation)
    console.print(cube)
    solver = Solver(cube)
    solution = solver.solve()
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
