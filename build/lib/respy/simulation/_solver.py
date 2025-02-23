from respy.rrock import GridDelta,ResRock
from respy.fluid import OnePhase,BlackOil
from respy.setup import Conds,Time,ResInit

class BaseSolver():
    """
    The Base Class initializing reservoir flow in Rectangular Cuboids;
    """

    def __init__(self,grid,depth=None,tcomp=None):
        """
        grid   : It is a GridDelta instance.
        
        depth  : reservoir depth, ft
        tcomp  : total compressibility, 1/psi

        """

        """

        grid   : GridDelta instance

        rrock  : It is a ResRock instance or any other rock class that
                 calculates rock properties at any given pressure.

        fluid  : It is a Fluid instance or any other fluid class that
                 calculates fluid properties at any given pressure. Multiple
                 fluid phases may exist with only one mopbile phase; the
                 rest of phases is at irreducible saturation, not mobile.

        wcond : tuple of WellCond instance.

        bcond : tuple of BoundCond instance.

        """

        self.grid = grid

        self._depth = depth*0.3048
        self._tcomp = tcomp/6894.76

    def set_rrock(self,rrock:ResRock):

        self.rrock = rrock

    def set_fluid(self,fluid):

        self.fluid = fluid

    def set_wcond(self,*args):

        self.wcond = Conds(*args)

    def set_bcond(self,*args):

        self.bcond = Conds(*args)        

    def set_time(self,*args,**kwargs):
        
        self.time = Time(*args,**kwargs)

    def init(self,*args,refp=None,DWOC=None,DGOC=None):
        """Calculates the initial pressure
        
        pzero   : initial pressure in psi; If not defined, it will be
                  calculated from reference point and reservoir rock depths.

        refp    : reference point (depth:ft,pressure:psi)
        grad    : fluid hydrostatic gradient, psi/ft
        """

        self._press = numpy.zeros(self.shape)

        if pzero is None:
            pzero = refp[1]+grad*(self.block.depth-refp[0])

        self._press[:,0] = pzero*6894.76

        self.set_block()
        self.set_build()

    def set_block(self):

        self.block = Block(grid,depth=depth,tcomp=tcomp)

    def set_build(self):

        pass

    def __call__(self,press=None,tcurr=0,tstep=None):

        if press is None:
            press = self._pzero

        rrock = self.rrock(press)
        fluid = self.fluid(press)
        wcond = self.wcond(tcurr)
        bcond = self.bcond(tcurr)

        if tstep is None:
            tstep = self.time._steps[0]

        return rrock,fluid,wcond,bcond,tstep

    def solve(self,*args,**kwargs):

        pass

    @property
    def depth(self):
        return self._depth/0.3048

    @property
    def tcomp(self):
        return self._tcomp*6894.76

    @property
    def press(self):
        return self._press/6894.76
    
    @property
    def pzero(self):
        return self._press[:,0]/6894.76

    @property
    def shape(self):
        return (self.grid.nums,self.time.nums+1)

    @property
    def isstatic(self):
        return all((self.rrock.isstatic,self.fluid.isstatic))

    @property
    def isdynamic(self):
        return all((self.rrock.isdynamic,self.fluid.isdynamic))

    @staticmethod
    def isstatic(prop):
        return not callable(prop)
