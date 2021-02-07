from sympy.liealgebras.cartan_type import CartanType
from sympy.matrices import Matrix

def test_type_B():
    c = CartanType("B3")
    m = Matrix(3, 3, [2, -1, 0, -1, 2, -2, 0, -1, 2])
    assert m == c.cartan_matrix()
    assert c.dimension() == 3
    assert c.roots() == 18
    assert c.simple_root(3) == Matrix([[0, 0, 1]])
    assert c.basis() == 3
    assert c.lie_algebra() == "so(6)"
    diag = "0---0=>=0\n1   2   3"
    assert c.dynkin_diagram() == diag
    assert c.positive_roots() ==  [
        Matrix([[1, 1, 0]]),
        Matrix([[1, 0, 1]]),
        Matrix([[0, 1, 1]]),
        Matrix([[1, 0, 0]]),
        Matrix([[0, 1, 0]]),
        Matrix([[1, 0, -1]]),
        Matrix([[0, 1, -1]]),
        Matrix([[0, 0, 1]]),
        Matrix([[1, -1, 0]])]

def test_type_B6():
    '''Testing numpy backend'''
    c = CartanType("B6")
    m = Matrix([
        [ 2, -1,  0,  0,  0,  0],
        [-1,  2, -1,  0,  0,  0],
        [ 0, -1,  2, -1,  0,  0],
        [ 0,  0, -1,  2, -1,  0],
        [ 0,  0,  0, -1,  2, -2],
        [ 0,  0,  0,  0, -1,  2]])
    assert m == c.cartan_matrix()

    simpleroot0 = Matrix([[1, -1, 0, 0, 0, 0]])

    assert simpleroot0 == c.simple_roots()[0]

    # take sample for brevity
    orbit = c.orbit(simpleroot0)
    assert len(orbit) == 60

    assert orbit[0] == Matrix([[1, 1, 0, 0, 0, 0]])
    assert orbit[6] == Matrix([[-1, 0, 0, 0, 0, 1]])

    assert Matrix(c.positive_roots()) == Matrix([
        [ 1,  1,  0,  0,  0,  0],
        [ 1,  0,  1,  0,  0,  0],
        [ 0,  1,  1,  0,  0,  0],
        [ 1,  0,  0,  1,  0,  0],
        [ 0,  1,  0,  1,  0,  0],
        [ 1,  0,  0,  0,  1,  0],
        [ 0,  1,  0,  0,  1,  0],
        [ 0,  0,  1,  1,  0,  0],
        [ 1,  0,  0,  0,  0,  1],
        [ 0,  1,  0,  0,  0,  1],
        [ 0,  0,  1,  0,  1,  0],
        [ 1,  0,  0,  0,  0,  0],
        [ 0,  1,  0,  0,  0,  0],
        [ 0,  0,  1,  0,  0,  1],
        [ 0,  0,  0,  1,  1,  0],
        [ 1,  0,  0,  0,  0, -1],
        [ 0,  1,  0,  0,  0, -1],
        [ 0,  0,  1,  0,  0,  0],
        [ 0,  0,  0,  1,  0,  1],
        [ 1,  0,  0,  0, -1,  0],
        [ 0,  1,  0,  0, -1,  0],
        [ 0,  0,  1,  0,  0, -1],
        [ 0,  0,  0,  1,  0,  0],
        [ 0,  0,  0,  0,  1,  1],
        [ 1,  0,  0, -1,  0,  0],
        [ 0,  1,  0, -1,  0,  0],
        [ 0,  0,  1,  0, -1,  0],
        [ 0,  0,  0,  1,  0, -1],
        [ 0,  0,  0,  0,  1,  0],
        [ 1,  0, -1,  0,  0,  0],
        [ 0,  1, -1,  0,  0,  0],
        [ 0,  0,  1, -1,  0,  0],
        [ 0,  0,  0,  1, -1,  0],
        [ 0,  0,  0,  0,  1, -1],
        [ 0,  0,  0,  0,  0,  1],
        [ 1, -1,  0,  0,  0,  0]
    ])
