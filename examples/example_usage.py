from cube_solver import Cube, Maneuver, Kociemba, apply_maneuver

scramble = Maneuver.random()
print(f"Scramble: {scramble}")

cube = Cube(scramble)
print(cube)
print(f"Cube: {repr(cube)}")

solver = Kociemba()
solution = solver.solve(cube)
assert isinstance(solution, Maneuver)
assert apply_maneuver(cube, solution).is_solved
print(f"Solution: {solution} ({len(solution)})")
