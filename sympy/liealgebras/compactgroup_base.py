from sympy.core import Basic
from sympy.core.backend import Matrix
from .cartan_base import Standard_Cartan


class CompactLieGroupBase(Basic):
    """
    Concrete base type for all Compact Lie Groups.
    """
    def __new__(cls, grp, n, rootsystem):
        obj = Basic.__new__(cls, n)
        obj.n = n
        obj.grp = grp
        obj._rootsystem = rootsystem

        return obj

    @property
    def rootsystem(self) -> Standard_Cartan:
        """
        Returns the associated lie algebra root system.
        
        All roots and weights are defaulted to the orthogonal basis.
        """
        return self._rootsystem

    def dimension(self):
        """
        Returns the dimension of the group.
        """
        return self.n

    def rank(self):
        """
        Returns the rank of the group.
        """
        return self.rootsystem.rank()

    def cartan_matrix(self):
        """
        Returns the Cartan Matrix for the associated lie algebra.
        """
        return self.rootsystem.cartan_matrix()

    def omega_matrix(self):
        """
        Returns the matrix of fundamental weights.
        """
        return self.cocartan_matrix().pinv().T

    def cocartan_matrix(self):
        """
        Returns the matrix built out of the co-roots from the Cartan Matrix.
        """
        return Matrix([2 * x / x.dot(x) for x in self.rootsystem.simple_roots()])

    def rotate_basis(self, X, src_basis, tgt_basis):
        """
        Returns rotated src_basis X to tgt_basis. The possible choice of bases are
        {"ortho", "alpha", "omega"}.

        Explanation
        ===========

        In representing the lie algebra, there are generally 3 common
        bases used to express the roots and weights. We called them

        The ``α-basis``: where the roots and weights are in terms of the simple roots α_i.

        The ``ω-basis``: where the roots and weights in terms of the fundamental weights ω_i. 
        The coefficients in this basis are often called Dynkin labels.

        The ``orthogonal-basis``: where the root/weight-space are embedded into a bigger Euclidean space.
            
        """
        if tgt_basis == "ortho":
            return self.to_orthogonal(X, src_basis)
        if tgt_basis == "alpha":
            return self.to_alpha(X, src_basis)
        if tgt_basis == "omega":
            return self.to_omega(X, src_basis)
        
        raise ValueError("Unknown basis")

    def to_orthogonal(self, X, basis):
        """
        Returns the representation X under this group in the orthogonal basis.
        """

        if basis == "ortho":
            return X
        elif basis == "alpha":
            return self.omega_matrix.T * self.cartan_matrix().T * X
        elif basis == "omega":
            return self.omega_matrix.T * X
        else:
            raise ValueError("Unknown basis")

    def to_alpha(self, X, basis):
        """
        Returns the representation X under this group in the alpha basis.
        """
        
        if basis == "ortho":
            return self.cartan_matrix().pinv().T * self.omega_matrix().pinv().T * X
        elif basis == "alpha":
            return X
        elif basis == "omega":
            return self.cartan_matrix().pinv().T * X
        else:
            raise ValueError("Unknown basis")

    def to_omega(self, X, basis):
        """
        Returns the representation X under this group in the omega basis.
        """
        if basis == "ortho":
            return self.omega_matrix().pinv().T * X
        elif basis == "alpha":
            return self.cartan_matrix().T * X
        elif basis == "omega":
            return X
        else:
            raise ValueError("Unknown basis")





    