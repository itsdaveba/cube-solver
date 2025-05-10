"""Console script for cube_solver."""
from cube_solver import Cube, Solver

import click
import typer
from rich.console import Console
from typing_extensions import Annotated

app = typer.Typer(no_args_is_help=True, add_completion=False, rich_markup_mode=None)
console = Console()


@app.command()
def scramble(length: Annotated[int, typer.Option("--length", "-l", help="Scramble length.  [default: 25]", show_default=False)] = 25):
    """Generate a random scramble."""
    scramble = Cube.generate_scramble(length)
    console.print(scramble)
    cube = Cube(scramble)
    console.print(cube)

@click.command(help="Solve a cube.")
@click.argument("scramble", required=False)
@click.option("-r", "--random", "length", is_flag=False, flag_value=25, default=0, help="Generate a random scramble of the specified length.  [default: 25]")
def solve(scramble, length):
    if (scramble is None) != bool(length):
        console.print("You must provide either 'scramble' or the '-r' / '--random' option.")
        return
    if length:
        scramble = Cube.generate_scramble(length)
    console.print("Scramble:", scramble)
    cube = Cube(scramble)
    console.print(cube)
    solver = Solver(transition_tables=True, pruning_tables=True)
    solution = solver.thistlethwaite(cube)
    console.print("Solution:", solution)

@app.callback()
def callback():
    """
    Cube Solver
    """


app = typer.main.get_command(app)
app.add_command(solve, "solve")

if __name__ == "__main__":
    app()
