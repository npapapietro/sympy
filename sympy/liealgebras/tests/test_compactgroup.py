from sympy.liealgebras.compactgroup_type import CompactLieGroup
from sympy.liealgebras.type_a import TypeA
from sympy.core.backend import Matrix
from sympy import I

def test_suN_properties():
    roots = TypeA(2)
    su3 = CompactLieGroup("SU3")

    assert roots.cartan_matrix() == su3.cartan_matrix
    assert roots.rank() == su3.rank

def test_suN_matricies():
    su3 = CompactLieGroup("SU3")

    omega_su3 = Matrix([
        ["2/3", "-1/3", "-1/3"],
        ["1/3", "1/3"," -2/3"]])

    cocartan_su3 = Matrix([
        ["1", "-1", "0"],
        ["0", "1"," -1"]])

    quadtratic = Matrix([
        ["2/3", "1/3"],
        ["1/3", "2/3"]])

    assert su3.omega_matrix == omega_su3
    assert su3.cocartan_matrix == cocartan_su3
    assert su3.quadratic_form == quadtratic


def test_suN_basis_rotations():
    su3 = CompactLieGroup("SU3")

    fundamental_omega = [1, 0]
    assert Matrix(fundamental_omega) == su3.to_omega(fundamental_omega, "omega")
    assert Matrix(["2/3", "-1/3", "-1/3"]) == su3.to_orthogonal(fundamental_omega, "omega")
    assert Matrix(["2/3", "1/3"]) == su3.to_alpha(fundamental_omega, "omega")

    fundamental_alpha = ["2/3", "1/3"]
    assert Matrix([1, 0]) == su3.to_omega(fundamental_alpha, "alpha")
    assert Matrix(["2/3", "-1/3", "-1/3"]) == su3.to_orthogonal(fundamental_alpha, "alpha")
    assert Matrix(fundamental_alpha) == su3.to_alpha(fundamental_alpha, "alpha")

    fundamental_ortho = ["2/3", "-1/3", "-1/3"]
    assert Matrix([1, 0]) == su3.to_omega(fundamental_ortho, "ortho")
    assert Matrix(fundamental_ortho) == su3.to_orthogonal(fundamental_ortho, "ortho")
    assert Matrix(["2/3", "1/3"]) == su3.to_alpha(fundamental_ortho, "ortho")

def test_suN_rep_dim():
    su3 = CompactLieGroup("SU3")

    f = [1,0]
    af = [0,1]
    ad = [1,1]

    assert su3.rep_dimension(f) == 3
    assert su3.rep_dimension(af) == 3
    assert su3.rep_dimension(ad) == 8


    su4 = CompactLieGroup("SU4")

    f = [1,0,0]
    af = [0,0,1]
    ad = [1,0,1]

    assert su4.rep_dimension(f) == 4
    assert su4.rep_dimension(af) == 4
    assert su4.rep_dimension(ad) == 15

def test_suN_gellmann():
    su2 = CompactLieGroup("SU2")
    pauli_matrices = [Matrix([
        [0, 1],
        [1, 0]]),
        Matrix([
        [0, -I],
        [I,  0]]),
        Matrix([
        [1,  0],
        [0, -1]])
    ]
    for a,b in zip(pauli_matrices, su2.gellmann_matrices):
        assert a == b
