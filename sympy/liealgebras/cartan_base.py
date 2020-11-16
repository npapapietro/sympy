from sympy.core import Basic
from sympy.core.backend import zeros, Matrix

class Standard_Cartan(Basic):
    """
    Semi-Concrete base class for Cartan types such as A4, etc
    """

    def __new__(cls, series, n):
        obj = Basic.__new__(cls, series, n)
        obj.n = n
        obj.series = series
        return obj

    def rank(self):
        """
        Returns the rank of the Lie algebra
        """
        return self.n

    def series(self):
        """
        Returns the type of the Lie algebra
        """
        return self.series

    def cartan_matrix(self):
        """
        Virtual method for generating the cartan matrix for this algebra.
        """
        r = self.rank()
        cartan_matrix = zeros(r,r)
        for i, sr_i in enumerate(self.simple_roots()):
            for j, sr_j in enumerate(self.simple_roots()):
                cartan_matrix[j,i] = 2 * sr_i.dot(sr_j) / sr_i.dot(sr_i)
        return cartan_matrix

    def simple_root(self, i):
        """
        Returns the i'th simple root in the orthogonal basis.
        """
        pass

    def simple_roots(self):
        """
        Returns the simple roots of the algebra.
        """
        return [Matrix(self.simple_root(i+1)).T for i in range(self.n)]

    
