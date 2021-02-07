from sympy import Matrix, S
from sympy.liealgebras.cartan_type import CartanType

def test_type_E():
    c = CartanType("E6")
    assert c.dimension() == 6
    assert c.roots() == 72
    assert c.basis() == 78
    diag = " "*8 + "6\n" + " "*8 + "0\n" + " "*8 + "|\n" + " "*8 + "|\n"
    diag += "---".join("0" for i in range(1, 6))+"\n"
    diag += "1   " + "   ".join(str(i) for i in range(2, 6))
    assert c.dynkin_diagram() == diag


def test_simpleroots():
    c = CartanType("E6")
    s = c.simple_roots()

    assert Matrix(s) == Matrix([
        [S.Half, -S.Half, -S.Half,
        -S.Half, -S.Half, -S.Half,
        -S.Half, S.Half],
        [-1, 1, 0, 0, 0, 0, 0, 0],
        [0, -1, 1, 0, 0, 0, 0, 0],
        [0, 0, -1, 1, 0, 0, 0, 0],
        [0, 0, 0, -1, 1, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0]])


    c = CartanType("E7")
    s = c.simple_roots()

    assert Matrix(s) == Matrix([
        [S.Half, -S.Half, -S.Half,
        -S.Half, -S.Half, -S.Half,
        -S.Half, S.Half],
       [-1, 1, 0, 0, 0, 0, 0, 0],
       [0, -1, 1, 0, 0, 0, 0, 0],
       [0, 0, -1, 1, 0, 0, 0, 0],
       [0, 0, 0, -1, 1, 0, 0, 0],
       [0, 0, 0, 0, -1, 1, 0, 0],
       [1, 1, 0, 0, 0, 0, 0, 0]])

    c = CartanType("E8")
    s = c.simple_roots()

    assert Matrix(s) == Matrix([
        [S.Half, -S.Half, -S.Half,
        -S.Half, -S.Half, -S.Half,
        -S.Half, S.Half],
        [-1, 1, 0, 0, 0, 0, 0, 0],
        [0, -1, 1, 0, 0, 0, 0, 0],
        [0, 0, -1, 1, 0, 0, 0, 0],
        [0, 0, 0, -1, 1, 0, 0, 0],
        [0, 0, 0, 0, -1, 1, 0, 0],
        [0, 0, 0, 0, 0, -1, 1, 0],
        [1, 1, 0, 0, 0, 0, 0, 0]])

def test_Cartan():
    c = CartanType("E6")
    s = c.cartan_matrix()

    assert Matrix(s) == Matrix([
        [ 2, -1,  0,  0,  0,  0],
        [-1,  2, -1,  0,  0,  0],
        [ 0, -1,  2, -1,  0, -1],
        [ 0,  0, -1,  2, -1,  0],
        [ 0,  0,  0, -1,  2,  0],
        [ 0,  0, -1,  0,  0,  2]])


    c = CartanType("E7")
    s = c.cartan_matrix()

    assert Matrix(s) == Matrix([
        [ 2, -1,  0,  0,  0,  0,  0],
        [-1,  2, -1,  0,  0,  0,  0],
        [ 0, -1,  2, -1,  0,  0, -1],
        [ 0,  0, -1,  2, -1,  0,  0],
        [ 0,  0,  0, -1,  2, -1,  0],
        [ 0,  0,  0,  0, -1,  2,  0],
        [ 0,  0, -1,  0,  0,  0,  2]])

    c = CartanType("E8")
    s = c.cartan_matrix()

    assert Matrix(s) == Matrix([
        [ 2, -1,  0,  0,  0,  0,  0,  0],
        [-1,  2, -1,  0,  0,  0,  0,  0],
        [ 0, -1,  2, -1,  0,  0,  0, -1],
        [ 0,  0, -1,  2, -1,  0,  0,  0],
        [ 0,  0,  0, -1,  2, -1,  0,  0],
        [ 0,  0,  0,  0, -1,  2, -1,  0],
        [ 0,  0,  0,  0,  0, -1,  2,  0],
        [ 0,  0, -1,  0,  0,  0,  0,  2]])

def test_orbit():
    """These results are compared to outputs of
    Mathematica's implementation of lie algebras"""
    # ignoring on large tests without numpy
    try:
        import numpy # noqa
    except ImportError:
        return

    c = CartanType("E6")
    s = c.simple_roots()
    orbit = c.orbit(s[0])

    assert len(orbit) == 72
    assert orbit[0] == Matrix([[1, 1, 0, 0, 0, 0, 0, 0]])
    assert orbit[10] == Matrix([
        [-S.Half, -S.Half, -S.Half, -S.Half,
         S.Half, -S.Half, -S.Half, S.Half]])
    assert orbit[30] == Matrix([[0, -1, 1, 0, 0, 0, 0, 0]])


    c = CartanType("E7")
    s = c.simple_roots()
    orbit = c.orbit(s[0])

    assert len(orbit) == 126
    assert orbit[0] == Matrix([[1, 1, 0, 0, 0, 0, 0, 0]])
    assert orbit[49] == Matrix([
        [0, -1, 0, 1, 0, 0, 0, 0]])
    assert orbit[99] == Matrix([
        [S.Half, -S.Half, S.Half, S.Half,
        S.Half, S.Half, -S.Half, S.Half]])

    c = CartanType("E8")
    s = c.simple_roots()
    orbit = c.orbit(s[0])

    assert len(orbit) == 240
    assert orbit[0] == Matrix([[1, 1, 0, 0, 0, 0, 0, 0]])
    assert orbit[49] == Matrix([
        [-S.Half, S.Half, -S.Half, -S.Half,
        -S.Half, S.Half, -S.Half, -S.Half]])
    assert orbit[-1] == Matrix([[1, 0, 1, 0, 0, 0, 0, 0]])

def test_positive_roots():
    """These results are compared to outputs of
    Mathematica's implementation of lie algebras"""

    # ignoring on large tests
    try:
        import numpy # noqa
    except ImportError:
        return
    c = CartanType("E6")
    p = c.positive_roots()

    assert Matrix(p) == Matrix([
        [S.Half, S.Half, S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, S.Half, -S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, S.Half, S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [S.Half, -S.Half, -S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [0, 0, 0, 1, 1, 0, 0, 0],
        [-S.Half, S.Half, S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [S.Half, -S.Half, S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [0, 0, 1, 0, 1, 0, 0, 0],
        [S.Half, -S.Half, S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [S.Half, S.Half, -S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [0, 0, 1, 1, 0, 0, 0, 0],
        [0, 1, 0, 0, 1, 0, 0, 0],
        [S.Half, S.Half, -S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, -S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [0, 1, 0, 1, 0, 0, 0, 0],
        [-1, 0, 0, 0, 1, 0, 0, 0],
        [1, 0, 0, 0, 1, 0, 0, 0],
        [S.Half, S.Half, S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, -S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [0, 1, 1, 0, 0, 0, 0, 0],
        [-1, 0, 0, 1, 0, 0, 0, 0],
        [1, 0, 0, 1, 0, 0, 0, 0],
        [0, -1, 0, 0, 1, 0, 0, 0],
        [-S.Half, -S.Half, S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [-1, 0, 1, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 0, 0, 0, 0],
        [0, -1, 0, 1, 0, 0, 0, 0],
        [0, 0, -1, 0, 1, 0, 0, 0],
        [-S.Half, S.Half, -S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [-1, 1, 0, 0, 0, 0, 0, 0],
        [0, -1, 1, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, -1, 1, 0, 0, 0, 0],
        [0, 0, 0, -1, 1, 0, 0, 0],
        [S.Half, -S.Half, -S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half],])

    c = CartanType("E7")
    p = c.positive_roots()

    assert Matrix(p) == Matrix([
        [0, 0, 0, 0, 0, 0, -1, 1],
        [-S.Half, S.Half, S.Half, S.Half, S.Half, S.Half, -S.Half, S.Half],
        [S.Half, -S.Half, S.Half, S.Half, S.Half, S.Half, -S.Half, S.Half],
        [S.Half, S.Half, -S.Half, S.Half, S.Half, S.Half, -S.Half, S.Half],
        [S.Half, S.Half, S.Half, -S.Half, S.Half, S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, -S.Half, S.Half, S.Half, S.Half, -S.Half, S.Half],
        [S.Half, S.Half, S.Half, S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, S.Half, -S.Half, S.Half, S.Half, -S.Half, S.Half],
        [S.Half, S.Half, S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, S.Half, S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [-S.Half, S.Half, -S.Half, -S.Half, S.Half, S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, S.Half, -S.Half, S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [S.Half, -S.Half, -S.Half, -S.Half, S.Half, S.Half, -S.Half, S.Half],
        [0, 0, 0, 0, 1, 1, 0, 0],
        [-S.Half, S.Half, -S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [S.Half, -S.Half, -S.Half, S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [0, 0, 0, 1, 0, 1, 0, 0],
        [-S.Half, S.Half, S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [S.Half, -S.Half, -S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [0, 0, 0, 1, 1, 0, 0, 0],
        [0, 0, 1, 0, 0, 1, 0, 0],
        [-S.Half, S.Half, S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [S.Half, -S.Half, S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [0, 0, 1, 0, 1, 0, 0, 0],
        [0, 1, 0, 0, 0, 1, 0, 0],
        [S.Half, -S.Half, S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [S.Half, S.Half, -S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [0, 0, 1, 1, 0, 0, 0, 0],
        [0, 1, 0, 0, 1, 0, 0, 0],
        [-1, 0, 0, 0, 0, 1, 0, 0],
        [1, 0, 0, 0, 0, 1, 0, 0],
        [S.Half, S.Half, -S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, -S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [0, 1, 0, 1, 0, 0, 0, 0],
        [-1, 0, 0, 0, 1, 0, 0, 0],
        [1, 0, 0, 0, 1, 0, 0, 0],
        [0, -1, 0, 0, 0, 1, 0, 0],
        [S.Half, S.Half, S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, -S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [0, 1, 1, 0, 0, 0, 0, 0],
        [-1, 0, 0, 1, 0, 0, 0, 0],
        [1, 0, 0, 1, 0, 0, 0, 0],
        [0, -1, 0, 0, 1, 0, 0, 0],
        [0, 0, -1, 0, 0, 1, 0, 0],
        [-S.Half, -S.Half, S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [-1, 0, 1, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 0, 0, 0, 0],
        [0, -1, 0, 1, 0, 0, 0, 0],
        [0, 0, -1, 0, 1, 0, 0, 0],
        [0, 0, 0, -1, 0, 1, 0, 0],
        [-S.Half, S.Half, -S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [-1, 1, 0, 0, 0, 0, 0, 0],
        [0, -1, 1, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, -1, 1, 0, 0, 0, 0],
        [0, 0, 0, -1, 1, 0, 0, 0],
        [0, 0, 0, 0, -1, 1, 0, 0],
        [S.Half, -S.Half, -S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half],
    ])

    c = CartanType("E8")
    p = c.positive_roots()

    assert Matrix(p) == Matrix([
        [0, 0, 0, 0, 0, 0, 1, 1],
        [0, 0, 0, 0, 0, 1, 0, 1],
        [0, 0, 0, 0, 1, 0, 0, 1],
        [0, 0, 0, 1, 0, 0, 0, 1],
        [0, 0, 1, 0, 0, 0, 0, 1],
        [0, 1, 0, 0, 0, 0, 0, 1],
        [-1, 0, 0, 0, 0, 0, 0, 1],
        [1, 0, 0, 0, 0, 0, 0, 1],
        [S.Half, S.Half, S.Half, S.Half, S.Half, S.Half, S.Half, S.Half],
        [0, -1, 0, 0, 0, 0, 0, 1],
        [-S.Half, -S.Half, S.Half, S.Half, S.Half, S.Half, S.Half, S.Half],
        [0, 0, -1, 0, 0, 0, 0, 1],
        [-S.Half, S.Half, -S.Half, S.Half, S.Half, S.Half, S.Half, S.Half],
        [0, 0, 0, -1, 0, 0, 0, 1],
        [-S.Half, S.Half, S.Half, -S.Half, S.Half, S.Half, S.Half, S.Half],
        [S.Half, -S.Half, -S.Half, S.Half, S.Half, S.Half, S.Half, S.Half],
        [0, 0, 0, 0, -1, 0, 0, 1],
        [-S.Half, S.Half, S.Half, S.Half, -S.Half, S.Half, S.Half, S.Half],
        [S.Half, -S.Half, S.Half, -S.Half, S.Half, S.Half, S.Half, S.Half],
        [0, 0, 0, 0, 0, -1, 0, 1],
        [-S.Half, S.Half, S.Half, S.Half, S.Half, -S.Half, S.Half, S.Half],
        [S.Half, -S.Half, S.Half, S.Half, -S.Half, S.Half, S.Half, S.Half],
        [S.Half, S.Half, -S.Half, -S.Half, S.Half, S.Half, S.Half, S.Half],
        [0, 0, 0, 0, 0, 0, -1, 1],
        [-S.Half, S.Half, S.Half, S.Half, S.Half, S.Half, -S.Half, S.Half],
        [S.Half, -S.Half, S.Half, S.Half, S.Half, -S.Half, S.Half, S.Half],
        [S.Half, S.Half, -S.Half, S.Half, -S.Half, S.Half, S.Half, S.Half],
        [-S.Half, -S.Half, -S.Half, -S.Half, S.Half, S.Half, S.Half, S.Half],
        [S.Half, -S.Half, S.Half, S.Half, S.Half, S.Half, -S.Half, S.Half],
        [S.Half, S.Half, -S.Half, S.Half, S.Half, -S.Half, S.Half, S.Half],
        [S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half, S.Half, S.Half],
        [-S.Half, -S.Half, -S.Half, S.Half, -S.Half, S.Half, S.Half, S.Half],
        [S.Half, S.Half, -S.Half, S.Half, S.Half, S.Half, -S.Half, S.Half],
        [S.Half, S.Half, S.Half, -S.Half, S.Half, -S.Half, S.Half, S.Half],
        [-S.Half, -S.Half, -S.Half, S.Half, S.Half, -S.Half, S.Half, S.Half],
        [-S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half, S.Half, S.Half],
        [S.Half, S.Half, S.Half, -S.Half, S.Half, S.Half, -S.Half, S.Half],
        [S.Half, S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half, S.Half],
        [-S.Half, -S.Half, -S.Half, S.Half, S.Half, S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, S.Half, -S.Half, S.Half, -S.Half, S.Half, S.Half],
        [-S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half, S.Half, S.Half],
        [S.Half, S.Half, S.Half, S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, S.Half, -S.Half, S.Half, S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half, S.Half],
        [-S.Half, S.Half, -S.Half, -S.Half, S.Half, -S.Half, S.Half, S.Half],
        [S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half, S.Half, S.Half],
        [0, 0, 0, 0, 0, 1, 1, 0],
        [S.Half, S.Half, S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, S.Half, S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [-S.Half, S.Half, -S.Half, -S.Half, S.Half, S.Half, -S.Half, S.Half],
        [-S.Half, S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half, S.Half],
        [S.Half, -S.Half, -S.Half, -S.Half, S.Half, -S.Half, S.Half, S.Half],
        [0, 0, 0, 0, 1, 0, 1, 0],
        [-S.Half, -S.Half, S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, S.Half, -S.Half, S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [-S.Half, S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half, S.Half],
        [S.Half, -S.Half, -S.Half, -S.Half, S.Half, S.Half, -S.Half, S.Half],
        [S.Half, -S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half, S.Half],
        [0, 0, 0, 0, 1, 1, 0, 0],
        [0, 0, 0, 1, 0, 0, 1, 0],
        [-S.Half, S.Half, -S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [S.Half, -S.Half, -S.Half, S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [S.Half, -S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half, S.Half],
        [0, 0, 0, 1, 0, 1, 0, 0],
        [0, 0, 1, 0, 0, 0, 1, 0],
        [-S.Half, S.Half, S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [S.Half, -S.Half, -S.Half, S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [S.Half, S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half, S.Half],
        [0, 0, 0, 1, 1, 0, 0, 0],
        [0, 0, 1, 0, 0, 1, 0, 0],
        [0, 1, 0, 0, 0, 0, 1, 0],
        [-S.Half, S.Half, S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [S.Half, -S.Half, S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half, S.Half],
        [0, 0, 1, 0, 1, 0, 0, 0],
        [0, 1, 0, 0, 0, 1, 0, 0],
        [-1, 0, 0, 0, 0, 0, 1, 0],
        [1, 0, 0, 0, 0, 0, 1, 0],
        [S.Half, -S.Half, S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [S.Half, S.Half, -S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half, -S.Half, S.Half],
        [0, 0, 1, 1, 0, 0, 0, 0],
        [0, 1, 0, 0, 1, 0, 0, 0],
        [-1, 0, 0, 0, 0, 1, 0, 0],
        [1, 0, 0, 0, 0, 1, 0, 0],
        [0, -1, 0, 0, 0, 0, 1, 0],
        [S.Half, S.Half, -S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, -S.Half, -S.Half, S.Half, -S.Half, -S.Half, S.Half],
        [0, 1, 0, 1, 0, 0, 0, 0],
        [-1, 0, 0, 0, 1, 0, 0, 0],
        [1, 0, 0, 0, 1, 0, 0, 0],
        [0, -1, 0, 0, 0, 1, 0, 0],
        [0, 0, -1, 0, 0, 0, 1, 0],
        [S.Half, S.Half, S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [-S.Half, -S.Half, -S.Half, S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [0, 1, 1, 0, 0, 0, 0, 0],
        [-1, 0, 0, 1, 0, 0, 0, 0],
        [1, 0, 0, 1, 0, 0, 0, 0],
        [0, -1, 0, 0, 1, 0, 0, 0],
        [0, 0, -1, 0, 0, 1, 0, 0],
        [0, 0, 0, -1, 0, 0, 1, 0],
        [-S.Half, -S.Half, S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [-1, 0, 1, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 0, 0, 0, 0],
        [0, -1, 0, 1, 0, 0, 0, 0],
        [0, 0, -1, 0, 1, 0, 0, 0],
        [0, 0, 0, -1, 0, 1, 0, 0],
        [0, 0, 0, 0, -1, 0, 1, 0],
        [-S.Half, S.Half, -S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half],
        [-1, 1, 0, 0, 0, 0, 0, 0],
        [0, -1, 1, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, -1, 1, 0, 0, 0, 0],
        [0, 0, 0, -1, 1, 0, 0, 0],
        [0, 0, 0, 0, -1, 1, 0, 0],
        [0, 0, 0, 0, 0, -1, 1, 0],
        [S.Half, -S.Half, -S.Half, -S.Half, -S.Half, -S.Half, -S.Half, S.Half],
    ])
