from cube_solver import Cube, Maneuver, Kociemba

scramble = Maneuver.random()
print(f"Scramble: {scramble}")

cube = Cube(scramble)
print(cube)
print(f"Cube: {repr(cube)}")

solver = Kociemba()
solution = solver.solve(cube)
assert solution is not None
assert solution == scramble.inverse
print(f"Solution: {solution} ({len(solution)})")
