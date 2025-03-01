from scipy.sparse import csr_matrix as csr
from scipy.sparse import linalg

from ._block import Block

class BaseSolver(Block):
    """
    The Base Class initializing reservoir flow in Rectangular Cuboids;
    """

    def __init__(self,grids,rrock,fluid,tcomp=None):
        """
        Initialization of Base Class Solver.

        Properties:
        ----------
        grids  : It is a GridDelta instance.

        rrock  : It is a ResRock instance providing rock propertie
            at any given pressure (or time step).

        fluid  : It is a Fluid instance providing fluid properties
            at any given pressure (or time step).

        tcomp  : total compressibility, 1/psi

        """
        super().__init__(grids,rrock,fluid,tcomp)

    @property
    def iterator(self):
        return Iterator

    @property
    def residual(self):
        return Residual
    

class Iterator:

    @staticmethod
    def explicit(mat,Pprev):
        """Explicit pressure solution returning Pnext."""
        RHS = csr.dot(-(mat._T+mat._J),Pprev)+mat._Q+mat._G

        return Pprev+linalg.spsolve(mat._A,RHS)

    @staticmethod
    def mixed(mat,Pprev,theta:float=0.5):
        """Mixed pressure solution returning Pnext."""
        LHS = (1-theta)*(mat._T+mat._J)+mat._A
        RHS = csr.dot(mat._A-theta*(mat._T+mat._J),Pprev)+mat._Q+mat._G

        return linalg.spsolve(LHS,RHS)

    @staticmethod
    def implicit(mat,Pprev):
        """Implicit pressure solution returning Pnext."""
        LHS = mat._T+mat._J+mat._A
        RHS = csr.dot(mat._A,Pprev)+mat._Q+mat._G
        
        return linalg.spsolve(LHS,RHS)

class Residual:

    @staticmethod
    def explicit(mat,Pprev,Pnext):
        """Returns residual vector in SI units, (m**3)/(sec)"""
        RHS = csr.dot(mat._A-(mat._T+mat._J),Pprev)+mat._Q+mat._G

        return -csr.dot(mat._A,Pnext)+RHS

    @staticmethod
    def implicit(mat,Pprev,Pnext):
        """Returns residual vector in SI units, (m**3)/(sec)"""
        LHS = mat._T+mat._J+mat._A
        RHS = csr.dot(mat._A,Pprev)+mat._Q+mat._G

        return -csr.dot(LHS,Pnext)+RHS