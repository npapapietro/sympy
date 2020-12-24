from sympy.core.numbers import Rational
from sympy.liealgebras.cartan_type import CartanType
from sympy.liealgebras.cartan_base import Standard_Cartan
from sympy import eye, Matrix

def test_Standard_Cartan():
    c = CartanType("A4")
    assert c.rank == 4
    assert c.series == "A"
    m = Standard_Cartan("A", 2)
    assert m.rank == 2
    assert m.series == "A"
    b = CartanType("B12")
    assert b.rank == 12
    assert b.series == "B"

def test_CartanType():
    c = CartanType("A4")

    omega = c.omega_matrix()
    cocar = c.cocartan_matrix()

    assert cocar * omega.T == eye(4)

    c2 = CartanType("A2")
    assert c2.omega_matrix() == Matrix([
                        [Rational(2,3), Rational(-1,3), Rational(-1,3)],
                        [Rational(1,3),  Rational(1,3), Rational(-2,3)]])

    assert c2.fundamental_weight(1) == Matrix([[Rational(2,3), Rational(-1,3), Rational(-1,3)]])

def test_basis_transforms():

    a2 = CartanType("A2")

    # orthogonal basis
    sr1 = Matrix([[1, -1, 0]])

    sr1_omega = a2.to_omega(sr1, "orthogonal")
    sr1_omegat = a2.to_omega(sr1.T, "orthogonal")

    assert sr1_omega == Matrix([[2, -1]])
    assert sr1_omegat == Matrix([[2, -1]]).T

    sr1_alpha = a2.to_alpha(sr1, "orthogonal")
    sr1_alphat = a2.to_alpha(sr1.T, "orthogonal")

    assert sr1_alpha == Matrix([[1,0]])
    assert sr1_alphat == Matrix([[1,0]]).T

    sr1_alpha = a2.to_alpha(sr1_omega, "omega")
    sr1_alphat = a2.to_alpha(sr1_omega.T, "omega")

    assert sr1_alpha == Matrix([[1,0]])
    assert sr1_alphat == Matrix([[1,0]]).T

    sr1_omega = a2.to_omega(sr1_alpha, "alpha")
    sr1_omegat = a2.to_omega(sr1_alpha.T, "alpha")

    assert sr1_omega == Matrix([[2, -1]])
    assert sr1_omegat == Matrix([[2, -1]]).T


    sr1_ortho = a2.to_orthogonal(Matrix([[2, -1]]), "omega")
    sr1_orthot = a2.to_orthogonal(Matrix([[2, -1]]).T, "omega")

    assert sr1_ortho == Matrix([[1, -1, 0]])
    assert sr1_orthot == Matrix([[1, -1, 0]]).T

    sr1_ortho = a2.to_orthogonal(Matrix([[1,0]]), "alpha")
    sr1_orthot = a2.to_orthogonal(Matrix([[1,0]]).T, "alpha")

    assert sr1_ortho == Matrix([[1, -1, 0]])
    assert sr1_orthot == Matrix([[1, -1, 0]]).T

def test_orbit():
    g2 = CartanType("G2")

    sr = g2.simple_roots()[0]

    orbit_sr1 = g2.orbit(sr)

    assert orbit_sr1 == [
        Matrix([[-1, 0, 1]]),
        Matrix([[-1, 1, 0]]),
        Matrix([[0, -1, 1]]),
        Matrix([[0, 1, -1]]),
        Matrix([[1, -1, 0]]),
        Matrix([[1, 0, -1]])]
    
    a2 = CartanType("A2")

    orbit_a = a2.orbit(Matrix([[1,0]]), basis="alpha")

    assert orbit_a == [
        Matrix([[-1, -1]]),
        Matrix([[-1, 0]]),
        Matrix([[0, -1]]),
        Matrix([[0, 1]]),
        Matrix([[1, 0]]),
        Matrix([[1, 1]])]

# def test_weight_tower():

#     a2 = CartanType("A2")

#     s = a2.simple_roots()[0]

#     a2.weight_tower(s)