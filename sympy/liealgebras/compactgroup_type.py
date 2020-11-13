from sympy.core import Basic
from .su import SU

def group_parser(name: str):
    """Returns tuple that separates dimension from type of algebra/group

    Exampes
    =======
    >>> from sympy.liealgebras.compactgroup_type import group_parser
    >>> group_parser("SU2")
    ["SU", 2]
    """
    type_ = ""
    dim = ""
    for i in name:
        if i.isdigit():
            dim += i
        else:
            type_ += i
    return type_, int(dim)

class CompactLieGroup_generator(Basic):
    """Constructor for basic compact lie groups. Supported are
    SU(N), SO(N), SP(N) and any exceptional groups.
    """
    def __call__(self, *args):
        """Returns the Compact Lie group with associated roots and representations.
        Takes any of the compact lie groups (e.g 'Sp4' or 'SU3'). 

        Examples
        ========

        """
        c = args[0]
        if type(c) == list:
            letter, n = c[0], int(c[1])
        elif type(c) == str:
            letter, n = group_parser(c)
        else:
            raise TypeError("Argument must be a string (e.g. 'SU2') \
                or a list (e.g. ['A', 3])")

        if n < 0:
            raise ValueError("Lie algebra rank cannot be negative")

        if letter == "SU":
            return SU(n)       
        # elif letter == "SO":
        #     rootsystem = CartanType(["D", n / 2]) if n % 2 == 0 else CartanType(["B", (n - 1) / 2])  
        # elif letter == "Sp" and n % 2 == 0:
        #     rootsystem = CartanType(["C", n / 2])
        else:
            raise TypeError("Undefined Lie group")
        
        # return CompactLieGroupBase(letter, n, rootsystem)



CompactLieGroup = CompactLieGroup_generator()
