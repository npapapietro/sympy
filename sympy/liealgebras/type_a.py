from .cartan_base import Standard_Cartan
from sympy import Matrix

class TypeA(Standard_Cartan):
    """
    This class contains the information about
    the A series of simple Lie algebras.
    ====
    """

    def __new__(cls, n):
        if n < 1:
            raise ValueError("n can not be less than 1")
        return Standard_Cartan.__new__(cls, "A", n)


    def dimension(self):
        """Dimension of the vector space V underlying the Lie algebra

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType("A4")
        >>> c.dimension()
        5
        """
        return self.rank+1


    def basic_root(self, i, j):
        """
        This is a method just to generate roots
        with a 1 iin the ith position and a -1
        in the jth position.

        """

        n = self.rank
        root = [0]*(n+1)
        root[i] = 1
        root[j] = -1
        return Matrix([root])

    def simple_root(self, i):
        """
        Every lie algebra has a unique root system.
        Given a root system Q, there is a subset of the
        roots such that an element of Q is called a
        simple root if it cannot be written as the sum
        of two elements in Q.  If we let D denote the
        set of simple roots, then it is clear that every
        element of Q can be written as a linear combination
        of elements of D with all coefficients non-negative.

        In A_n the ith simple root is the root which has a 1
        in the ith position, a -1 in the (i+1)th position,
        and zeroes elsewhere.

        This method returns the ith simple root for the A series.

        Examples
        ========

        >>> from sympy.liealgebras.cartan_type import CartanType
        >>> c = CartanType("A4")
        >>> c.simple_root(1)
        Matrix([[1, -1, 0, 0, 0]])

        """
        n = self.rank
        if i < 1 or i > n:
            raise ValueError("Simple root %s does not exist for A_%s." % (i, n))
        return self.basic_root(i-1, i)


    def highest_root(self):
        """
        Returns the highest weight root for A_n
        """

        return self.basic_root(0, self.rank)

    def roots(self):
        """
        Returns the total number of roots for A_n
        """
        n = self.rank
        return n*(n+1)

    def basis(self):
        """
        Returns the number of independent generators of A_n
        """
        n = self.rank
        return n**2 - 1

    def lie_algebra(self):
        """
        Returns the Lie algebra associated with A_n
        """
        n = self.rank
        return "su(" + str(n + 1) + ")"

    def dynkin_diagram(self):
        n = self.rank
        diag = "---".join("0" for i in range(1, n+1)) + "\n"
        diag += "   ".join(str(i) for i in range(1, n+1))
        return diag

    def orbit(self, weight, stabilizer=None):
        """Returns the weyl orbit of the weight. If
        rank of the algebra is >5, numpy backend is used"""
        if self.rank >= 5:
            backend = "numpy"
            dtype = int
        else:
            backend = "sympy"
            dtype=object
        return super().orbit(weight, stabilizer=stabilizer, dtype=dtype, backend=backend)
