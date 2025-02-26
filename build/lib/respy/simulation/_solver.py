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