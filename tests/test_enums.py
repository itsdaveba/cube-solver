from cube_solver.cube.enums import Axis, Layer, Color, Face, Cubie, Move


def test_axis(response):
    assert hasattr(Axis, "NONE")
    assert len([*Axis.axes()]) == 5

    assert hasattr(Layer, "NONE")
    assert len([*Layer.layers()]) == 9

    assert hasattr(Color, "NONE")
    assert len([*Color.colors()]) == 6
    assert all(color == Color.from_char(color.char) for color in Color)

    assert hasattr(Face, "NONE")
    assert len([*Face.faces()]) == 6
    assert all(face == Face.from_char(face.char) for face in Face)

    assert hasattr(Cubie, "NONE")
    assert len([*Cubie.cubies()]) == 27
    assert len([*Cubie.corners()]) == 8
    assert all(cubie.is_corner for cubie in Cubie.corners())
    assert len([*Cubie.edges()]) == 12
    assert all(cubie.is_edge for cubie in Cubie.edges())
    assert len([*Cubie.centers()]) == 6
    assert all(cubie.is_center for cubie in Cubie.centers())
    assert all(cubie == Cubie.from_faces(cubie.faces) for cubie in Cubie)

    assert hasattr(Move, "NONE")
    assert len([*Move.moves()]) == 54
    assert len([*Move.face_moves()]) == 18
    assert all(move.is_face for move in Move.face_moves())
    assert len([*Move.slice_moves()]) == 9
    assert all(move.is_slice for move in Move.slice_moves())
    assert len([*Move.wide_moves()]) == 18
    assert all(move.is_wide for move in Move.wide_moves())
    assert len([*Move.rotations()]) == 9
    assert all(move.is_rotation for move in Move.rotations())
    assert all(move == Move.from_str(move.string) for move in Move)
    assert all(move == Move.from_str(move.string[0].lower() + move.string[2:]) for move in Move.wide_moves())
