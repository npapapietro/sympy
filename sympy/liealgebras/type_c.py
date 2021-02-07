from .cartan_base import Standard_Cartan
from sympy import Matrix

class TypeC(Standard_Cartan):

    def __new__(cls, n):
        if n < 3:
            raise ValueError("n can not be less than 3")
        return Standard_Cartan.__new__(cls, "C", n)


    def dimension(self):
        """Dimension of the vector space V underlying the Lie algebra

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType("C3")
        >>> c.dimension()
        3
        """
        n = self.rank
        return n

    def basic_root(self, i, j):
        """Generate roots with 1 in ith position and a -1 in jth position
        """
        n = self.rank
        root = [0]*n
        root[i] = 1
        root[j] = -1
        return Matrix([root])

    def simple_root(self, i):
        """The ith simple root for the C series

        Every lie algebra has a unique root system.
        Given a root system Q, there is a subset of the
        roots such that an element of Q is called a
        simple root if it cannot be written as the sum
        of two elements in Q.  If we let D denote the
        set of simple roots, then it is clear that every
        element of Q can be written as a linear combination
        of elements of D with all coefficients non-negative.

        In C_n, the first n-1 simple roots are the same as
        the roots in A_(n-1) (a 1 in the ith position, a -1
        in the (i+1)th position, and zeroes elsewhere).  The
        nth simple root is the root in which there is a 2 in
        the nth position and zeroes elsewhere.

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType("C3")
        >>> c.simple_root(2)
        Matrix([[0, 1, -1]])

        """

        n = self.rank
        if i < n:
            return self.basic_root(i-1,i)
        else:
            root = [0]*self.rank
            root[n-1] = 2
            return Matrix([root])

    def roots(self):
        """
        Returns the total number of roots for C_n"
        """

        n = self.rank
        return 2*(n**2)

    def basis(self):
        """
        Returns the number of independent generators of C_n
        """

        n = self.rank
        return n*(2*n + 1)

    def lie_algebra(self):
        """
        Returns the Lie algebra associated with C_n"
        """

        n = self.rank
        return "sp(" + str(2*n) + ")"

    def dynkin_diagram(self):
        n = self.rank
        diag = "---".join("0" for i in range(1, n)) + "=<=0\n"
        diag += "   ".join(str(i) for i in range(1, n+1))
        return diag

    def orbit(self, weight, stabilizer=None):
        """Returns the weyl orbit of the weight. If
        rank of the algebra is >5, numpy backend is used"""
        if self.rank >= 5:
            backend = "numpy"
            dtype = float
        else:
            backend = "sympy"
            dtype = object
        return super().orbit(weight, stabilizer=stabilizer, dtype=dtype, backend=backend)
