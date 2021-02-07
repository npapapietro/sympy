from sympy.liealgebras.cartan_type import CartanType
from sympy.matrices import Matrix



def test_type_D():
    c = CartanType("D4")
    m = Matrix(4, 4, [2, -1, 0, 0, -1, 2, -1, -1, 0, -1, 2, 0, 0, -1, 0 , 2])
    assert c.cartan_matrix() == m
    assert c.basis() == 6
    assert c.lie_algebra() == "so(8)"
    assert c.roots() == 24
    assert c.simple_root(3) == Matrix([[0, 0, 1, -1]])
    diag = "    3\n    0\n    |\n    |\n0---0---0\n1   2   4"
    assert diag == c.dynkin_diagram()
    assert c.positive_roots() == [Matrix([[1, 1, 0, 0]]),
        Matrix([[1, 0, 1, 0]]),
        Matrix([[0, 1, 1, 0]]),
        Matrix([[1, 0, 0, 1]]),
        Matrix([[1, 0, 0, -1]]),
        Matrix([[0, 1, 0, 1]]),
        Matrix([[0, 1, 0, -1]]),
        Matrix([[1, 0, -1, 0]]),
        Matrix([[0, 1, -1, 0]]),
        Matrix([[0, 0, 1, 1]]),
        Matrix([[0, 0, 1, -1]]),
        Matrix([[1, -1, 0, 0]])]

def test_type_D6():
    '''Testing numpy backend'''
    c = CartanType("D6")
    m = Matrix([
        [ 2, -1,  0,  0,  0,  0],
        [-1,  2, -1,  0,  0,  0],
        [ 0, -1,  2, -1,  0,  0],
        [ 0,  0, -1,  2, -1, -1],
        [ 0,  0,  0, -1,  2,  0],
        [ 0,  0,  0, -1,  0,  2]])
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
        [ 1,  0,  0,  0,  0, -1],
        [ 0,  1,  0,  0,  0,  1],
        [ 0,  1,  0,  0,  0, -1],
        [ 0,  0,  1,  0,  1,  0],
        [ 1,  0,  0,  0, -1,  0],
        [ 0,  1,  0,  0, -1,  0],
        [ 0,  0,  1,  0,  0,  1],
        [ 0,  0,  1,  0,  0, -1],
        [ 0,  0,  0,  1,  1,  0],
        [ 1,  0,  0, -1,  0,  0],
        [ 0,  1,  0, -1,  0,  0],
        [ 0,  0,  1,  0, -1,  0],
        [ 0,  0,  0,  1,  0,  1],
        [ 0,  0,  0,  1,  0, -1],
        [ 1,  0, -1,  0,  0,  0],
        [ 0,  1, -1,  0,  0,  0],
        [ 0,  0,  1, -1,  0,  0],
        [ 0,  0,  0,  1, -1,  0],
        [ 0,  0,  0,  0,  1,  1],
        [ 0,  0,  0,  0,  1, -1],
        [ 1, -1,  0,  0,  0,  0]
    ])
