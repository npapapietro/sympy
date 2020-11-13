from sympy.core import Basic

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
        pass

