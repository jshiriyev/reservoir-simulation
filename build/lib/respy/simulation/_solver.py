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
    def method(self):
        return Method

    @property
    def residual(self):
        return Residual
    

class Method:

    @staticmethod
    def explicit(Mprev,Pprev):
        """Explicit pressure solution returning Pnext."""
        RHS = csr.dot(-(Mprev._T+Mprev._J),Pprev)+Mprev._Q+Mprev._G

        return Pprev+linalg.spsolve(Mprev._A,RHS)

    @staticmethod
    def mixed(Mprev,Mnext,Pprev,theta:float=0.5):
        """Mixed pressure solution returning Pnext."""
        LHS = (1-theta)*(Mnext._T+Mnext._J)+Mnext._A
        RHS = csr.dot(Mprev._A-theta*(Mprev._T+Mprev._J),Pprev)+Mprev._Q+Mprev._G

        return linalg.spsolve(LHS,RHS)

    @staticmethod
    def implicit(Mnext,Pprev):
        """Implicit pressure solution returning Pnext."""
        LHS = Mnext._T+Mnext._J+Mnext._A
        RHS = csr.dot(Mnext._A,Pprev)+Mnext._Q+Mnext._G
        
        return linalg.spsolve(LHS,RHS)

class Residual:

    @staticmethod
    def explicit(Mprev,Pprev,Pnext):
        """Returns residual vector in SI units, (m**3)/(sec)"""
        RHS = csr.dot(Mprev._A-(Mprev._T+Mprev._J),Pprev)+Mprev._Q+Mprev._G

        return -csr.dot(Mprev._A,Pnext)+RHS

    @staticmethod
    def implicit(Mnext,Pprev,Pnext):
        """Returns residual vector in SI units, (m**3)/(sec)"""
        LHS = Mnext._T+Mnext._J+Mnext._A
        RHS = csr.dot(Mnext._A,Pprev)+Mnext._Q+Mnext._G

        return -csr.dot(LHS,Pnext)+RHS