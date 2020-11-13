from sympy.matrices.dense import ones
from sympy.core import Basic
from sympy.core.backend import Matrix, zeros
from .cartan_base import Standard_Cartan

def _check_type(X):
    if not isinstance(X, Matrix):
        return Matrix(X)
    else:
        return X
    

class CompactLieGroupBase(Basic):
    """
    Concrete base type for all Compact Lie Groups.
    """
    def __new__(cls, grp, n, rootsystem):
        obj = Basic.__new__(cls, n)
        obj.n = n
        obj.grp = grp
        obj._rootsystem = rootsystem

        # Lazy evaluation objects
        obj._cartan_matrix = None
        obj._omega_matrix = None
        obj._cocartan_matrix = None
        obj._quadratic_form = None
        obj._positive_roots = []

        return obj

    @property
    def rootsystem(self) -> Standard_Cartan:
        """
        Returns the associated lie algebra root system.
        
        All roots and weights are defaulted to the orthogonal basis.
        """
        return self._rootsystem

    @property
    def dimension(self):
        """
        Returns the dimension of the group.
        """
        return self.n

    @property
    def rank(self):
        """
        Returns the rank of the group.
        """
        return self.rootsystem.rank()

    @property
    def cartan_matrix(self):
        """
        Returns the Cartan Matrix for the associated lie algebra.
        """
        if self._cartan_matrix is None:
            self._cartan_matrix = self.rootsystem.cartan_matrix()
        return self._cartan_matrix

    @property
    def omega_matrix(self):
        """
        Returns the matrix of fundamental weights.
        """
        if self._omega_matrix is None:
            self._omega_matrix = self.cocartan_matrix.pinv().T
        return self._omega_matrix

    @property
    def cocartan_matrix(self):
        """
        Returns the matrix built out of the co-roots from the Cartan Matrix.
        """
        if self._cocartan_matrix is None:
            self._cocartan_matrix = Matrix([
                2 * x / x.dot(x) for x in self.rootsystem.simple_roots()])
        return self._cocartan_matrix

    @property
    def quadratic_form(self):
        """
        Returns the quadratic form of the simple roots as a diagonal matrix
        """
        if self._quadratic_form is None:
            rank = self.rank
            temp = zeros(rank, rank)

            for i, w in enumerate(self.rootsystem.simple_roots()):
                temp[i,i] = w.dot(w) / 2
            self._quadratic_form = self.cartan_matrix.pinv() * temp
            
        return self._quadratic_form

    @property
    def positive_roots(self):
        """
        Returns the positive roots as a list of 
        """
        if len(self._positive_roots) == 0:
            for root in self.rootsystem.positive_roots().values():
                self._positive_roots.append(Matrix(root).T)
        return self._positive_roots

    def rotate_basis(self, X, src_basis, tgt_basis):
        """
        Returns rotated src_basis X to tgt_basis. The possible choice of 
        bases are {"ortho", "alpha", "omega"}.

        Explanation
        ===========

        In representing the lie algebra, there are generally 3 common
        bases used to express the roots and weights. We called them

        The ``α-basis``: where the roots and weights are in terms of 
        the simple roots α_i.

        The ``ω-basis``: where the roots and weights in terms of the 
        fundamental weights ω_i. The coefficients in this basis are 
        often called Dynkin labels.

        The ``orthogonal-basis``: where the root/weight-space are 
        embedded into a bigger Euclidean space.
            
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
        X = _check_type(X)
        if basis == "ortho":
            return X
        elif basis == "alpha":
            return self.omega_matrix.T * self.cartan_matrix.T * X
        elif basis == "omega":
            return self.omega_matrix.T * X
        else:
            raise ValueError("Unknown basis")

    def to_alpha(self, X, basis):
        """
        Returns the representation X under this group in the alpha basis.
        """
        X = _check_type(X)
        if basis == "ortho":
            return self.cartan_matrix.pinv().T * self.omega_matrix.pinv().T * X
        elif basis == "alpha":
            return X
        elif basis == "omega":
            return self.cartan_matrix.pinv().T * X
        else:
            raise ValueError("Unknown basis")

    def to_omega(self, X, basis):
        """
        Returns the representation X under this group in the omega basis.
        """
        X = _check_type(X)
        if basis == "ortho":
            return self.omega_matrix.pinv().T * X
        elif basis == "alpha":
            return self.cartan_matrix.T * X
        elif basis == "omega":
            return X
        else:
            raise ValueError("Unknown basis")

    def rep_dimension(self, X, basis="omega"):
        """
        Returns the dimension of the representation. The basis of 
        the representation should be in "omega" because that is
        a common basis in literature for writing down the 
        fundamental representations.

        Examples
        ========
        >>> from sympy.liealgebras.su import SU
        >>> su3 = SU(3)
        >>> fundamental = [1,0]
        >>> antifundamental = [0,1]
        >>> adjoint = [1,1]
        >>> print(su3.rep_dimension(fundamental))
        3
        >>> print(su3.rep_dimension(antifundamental))
        3
        >>> print(su3.rep_dimension(adjoint))
        8

        """
        X = self.to_omega(X, basis)
        r = ones(self.rank, 1)
        dim = 1
        for i in self.positive_roots:
            i = self.to_omega(i.T, "ortho")
            n = i.T * self.quadratic_form * (X + r)
            d = r.T * self.quadratic_form * i
            dim *= n[0,0] / d[0,0]
        return dim



    
