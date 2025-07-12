from typer.testing import CliRunner

from cube_solver import Maneuver
from cube_solver.cli import app

runner = CliRunner()


def test_cli():
    # maneuver
    result = runner.invoke(app, ["maneuver"])
    assert result.exit_code == 0
    assert "WWWWOOOOGGGGRRRRBBBBYYYY" in result.output
    result = runner.invoke(app, ["maneuver", "U F2 R' D B2 L' x y2 z'"])
    assert result.exit_code == 0
    assert "RYBRBYBWOGORYGGOOYWWBWRG" in result.output
    result = runner.invoke(app, ["maneuver", "-c", "RYBRBYBWOGORYGGOOYWWBWRG", "U R2 F' R U2 F'"])
    assert result.exit_code == 0
    assert "OOOOBBBBYYYYGGGGWWWWRRRR" in result.output

    # scramble
    result = runner.invoke(app, ["scramble"])
    assert result.exit_code == 0
    assert len(Maneuver(result.output)) == 15
    result = runner.invoke(app, ["scramble", "-l", "0", "-v"])
    assert result.exit_code == 0
    assert "WWWWOOOOGGGGRRRRBBBBYYYY" in result.output
    result = runner.invoke(app, ["scramble", "-l", "20"])
    assert result.exit_code == 0
    assert len(Maneuver(result.output)) == 20
    result = runner.invoke(app, ["scramble", "--wca"])
    assert result.exit_code == 0
    assert len(Maneuver(result.output)) <= 15

    # solve
    result = runner.invoke(app, ["solve"])
    assert result.exit_code == 1
    assert "cube" in result.output
    assert "--scramble" in result.output
    assert "--random" in result.output
    result = runner.invoke(app, ["solve", "RYBRBYBWOGORYGGOOYWWBWRG", "-s", "U F2 R' D B2 L' x y2 z'"])
    assert result.exit_code == 1
    assert "cube" in result.output
    assert "--scramble" in result.output
    assert "--random" not in result.output
    result = runner.invoke(app, ["solve", "RYBRBYBWOGORYGGOOYWWBWRG", "-r"])
    assert result.exit_code == 1
    assert "cube" in result.output
    assert "--scramble" not in result.output
    assert "--random" in result.output
    result = runner.invoke(app, ["solve", "-s", "U F2 R' D B2 L' x y2 z'", "-r"])
    assert result.exit_code == 1
    assert "cube" not in result.output
    assert "--scramble" in result.output
    assert "--random" in result.output
    result = runner.invoke(app, ["solve", "RYBRBYBWOGORYGGOOYWWBWRG"])
    assert result.exit_code == 0
    assert result.output == "U R2 F' R U2 F'\n"
    result = runner.invoke(app, ["solve", "-s", "U F2 R' D B2 L' x y2 z'"])
    assert result.exit_code == 0
    assert result.output == "U R2 F' R U2 F'\n"
    result = runner.invoke(app, ["solve", "-r", "-v"])
    assert result.exit_code == 0
    assert "Solution:" in result.output
    assert "(" in result.output and ")" in result.output
    result = runner.invoke(app, ["solve", "-r", "-v", "-l", "0"])
    assert result.exit_code == 0
    assert "Solution: None" in result.output
    result = runner.invoke(app, ["solve", "-r", "-a", "kociemba", "-o", "-v"])
    assert result.exit_code == 0
    assert "Optimal:" in result.output
