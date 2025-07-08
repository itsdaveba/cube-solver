from typer.testing import CliRunner

from cube_solver import Maneuver
from cube_solver.cli import app

runner = CliRunner()


def test_cli():
    # maneuver
    result = runner.invoke(app, ["maneuver"])
    assert result.exit_code == 0
    assert "WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY" in result.output
    result = runner.invoke(app, ["maneuver", "U F2 R' D B2 L' M E2 S' Uw Fw2 Rw' Dw Bw2 Lw' u f2 r' d b2 l' x y2 z'"])
    assert result.exit_code == 0
    assert "YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOYO" in result.output
    result = runner.invoke(app, ["maneuver", "-c", "YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOYO",
                                 "U2 F2 L' B2 D' L' B U B D2 R2 F2 U L2 B2 U R2 F2 U B2 U' L2 B2"])
    assert result.exit_code == 0
    assert "YYYYYYYYYRRRRRRRRRGGGGGGGGGOOOOOOOOOBBBBBBBBBWWWWWWWWW" in result.output

    # scramble
    result = runner.invoke(app, ["scramble"])
    assert result.exit_code == 0
    assert len(Maneuver(result.output)) == 25
    result = runner.invoke(app, ["scramble", "-l", "0", "-v"])
    assert result.exit_code == 0
    assert "WWWWWWWWWOOOOOOOOOGGGGGGGGGRRRRRRRRRBBBBBBBBBYYYYYYYYY" in result.output
    result = runner.invoke(app, ["scramble", "-l", "30"])
    assert result.exit_code == 0
    assert len(Maneuver(result.output)) == 30
    result = runner.invoke(app, ["scramble", "--wca"])
    assert result.exit_code == 0
    assert len(Maneuver(result.output)) <= 25

    # solve
    result = runner.invoke(app, ["solve"])
    assert result.exit_code == 1
    assert "cube" in result.output
    assert "--scramble" in result.output
    assert "--random" in result.output
    result = runner.invoke(app, ["solve", "YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOYO",
                                 "-s", "U F2 R' D B2 L' M E2 S' Uw Fw2 Rw' Dw Bw2 Lw' u f2 r' d b2 l' x y2 z'"])
    assert result.exit_code == 1
    assert "cube" in result.output
    assert "--scramble" in result.output
    assert "--random" not in result.output
    result = runner.invoke(app, ["solve", "YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOYO", "-r"])
    assert result.exit_code == 1
    assert "cube" in result.output
    assert "--scramble" not in result.output
    assert "--random" in result.output
    result = runner.invoke(app, ["solve", "-s", "U F2 R' D B2 L' M E2 S' Uw Fw2 Rw' Dw Bw2 Lw' u f2 r' d b2 l' x y2 z'", "-r"])
    assert result.exit_code == 1
    assert "cube" not in result.output
    assert "--scramble" in result.output
    assert "--random" in result.output
    result = runner.invoke(app, ["solve", "YGWYYOBWWBGRGRRWGBYBGBGBRRYRYOOOBGOGGROYBWYRBWWROWWOYO"])
    assert result.exit_code == 0
    assert result.output == "U2 F2 L' B2 D' L' B U B D2 R2 F2 U L2 B2 U R2 F2 U B2 U' L2 B2\n"
    result = runner.invoke(app, ["solve", "-s", "U F2 R' D B2 L' M E2 S' Uw Fw2 Rw' Dw Bw2 Lw' u f2 r' d b2 l' x y2 z'"])
    assert result.exit_code == 0
    assert result.output == "U2 F2 L' B2 D' L' B U B D2 R2 F2 U L2 B2 U R2 F2 U B2 U' L2 B2\n"
    result = runner.invoke(app, ["solve", "-r", "-v"])
    assert result.exit_code == 0
    assert "Solution:" in result.output
    assert "(" in result.output and ")" in result.output
    result = runner.invoke(app, ["solve", "-r", "-v", "-t", "0"])
    assert result.exit_code == 0
    assert "Solution: None" in result.output
    result = runner.invoke(app, ["solve", "-s", "L2 U R D' B2 D2 F B D", "-o", "-v"])
    assert result.exit_code == 0
    assert "Optimal:" in result.output
    result = runner.invoke(app, ["solve", "-r", "-o", "-vv", "-t", "1"])
    assert result.exit_code == 0
    assert "Suboptimal:" in result.output
