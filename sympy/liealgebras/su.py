from sympy.core.backend import zeros, sqrt
from sympy import I, Sum

from .compactgroup_base import CompactLieGroupBase
from .type_a import TypeA

class SU(CompactLieGroupBase):

    def __new__(cls, n):
        if n < 1:
            raise ValueError("n cannot be less than 2")
        elif n < 2:
            raise ValueError("n cannot be less than 2, for U(1) use a phase.")
        
        obj = CompactLieGroupBase.__new__(cls, "SU", n, TypeA(n - 1))
        obj._gellmann = []
        return obj
    
    @property
    def gellmann_matrices(self):
        """
        Returns the generalized Gell-Mann matrices for the group. 
        
        The generalized Gell-Mann matrices are the n^2-1 matrices 
        generating the Lie algebra associated to the special 
        unitary group SU(n), n>=2. As their name suggests, these 
        matrices are intended to generalize both the standard 3×3 
        Gell-Mann matrices, which generate the Lie algebra associated 
        to SU(3), as well as the 2×2 Pauli matrices which generate 
        the Lie algebra associated to SU(2).

        Source: https://mathworld.wolfram.com/GeneralizedGell-MannMatrix.html

        Examples
        ========
        >>> from sympy.liealgebras import CompactLieGroup
        >>> su2 = CompactLieGroup("SU2")
        >>> su2.gellmann_matrices[2] # diagonal matrix
        Matrix([
        [1,  0],
        [0, -1]])
        """
        def Eij(i,j):
            _m = zeros(self.dimension)
            _m[i,j] = 1
            return _m

        if len(self._gellmann) == 0:
            n = self.dimension
            for i in range(0,n):
                for j in range(i+1, n):
                    symmetric = Eij(j,i) + Eij(i,j)
                    antisymmetric = - I * (Eij(i,j) - Eij(j,i))
                    self._gellmann += [symmetric, antisymmetric]

            for l in range(n-1):
                diag = zeros(n)
                for j in range(l+1):
                    diag += Eij(j,j) - (j+1) * Eij(l+1, l+1)
                self._gellmann += [sqrt(2) / sqrt((l+1)*(l+2)) * diag]
        return self._gellmann




    
