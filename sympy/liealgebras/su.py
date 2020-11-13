from sympy.core.backend import zeros

from .compactgroup_base import CompactLieGroupBase
from .type_a import TypeA

class SU(CompactLieGroupBase):

    def __new__(cls, n):
        if n < 1:
            raise ValueError("n cannot be less than 2")
        elif n < 2:
            raise ValueError("n cannot be less than 2, for U(1) use a phase.")
        
        return CompactLieGroupBase.__new__(cls, "SU", n, TypeA(n - 1))
    
    # def weight_dim(self, weight):
    #     if self._cocartan is None:
    #         self._generateMatrices()

    def _generateMatrices(self):
        rank = self._rootsystem.n
        cartan = self._rootsystem.cartan_matrix()

        cocartan = zeros(rank, rank + 1)
        for i in range(rank):
            cocartan
    