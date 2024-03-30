import sys

if __name__ == "__main__":
    # sys.path.append(r'C:\Users\javid.shiriyev\Documents\respy')
    sys.path.append(r'C:\Users\3876yl\Documents\respy')

from respy.rinit._time import Time

from respy.fluid import OnePhase, BlackOil

from respy.solver._onephase import OnePhase

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

        wconds : tuple of WellCond instance.

        bconds : tuple of BoundCond instance.

        """

        self.grid = grid

        self._depth = depth*0.3048
        self._tcomp = tcomp/6894.76

    def set_rrock(self,rrock:ResRock):

        self.rrock = rrock

    def set_fluid(self,fluid):

        self.fluid = fluid

    def set_wconds(self,wconds=None):

        self.__wconds = () if wconds is None else wconds

    def get_wconds(self,tcurr):

        return [cond for cond in self.__wconds if self.islive(cond,tcurr)]

    def set_bconds(self,bconds=None):

        self.__bconds = () if bconds is None else bconds

    def get_bconds(self,tcurr):

        return [cond for cond in self.__bconds if self.islive(cond,tcurr)]

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

        rrock  = self.rrock(press)
        fluid  = self.fluid(press)

        wconds = self.get_wconds(tcurr)
        bconds = self.get_bconds(tcurr)

        if tstep is None:
            tstep = self.time._steps[0]

        return rrock,fluid,wconds,bconds

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
    def tstat(self):
        return all((self.rstat,self.fstat))

    @property
    def rstat(self):
        return not callable(self.__rrock)

    @property
    def fstat(self):
        return not callable(self.__fluid)

    @staticmethod
    def isstat(prop):
        return not callable(prop)

    @staticmethod
    def islive(cond,tcurr):

        if cond._start>tcurr:
            return False

        if cond._stop is None:
            return True

        if cond._stop<=tcurr:
            return False

        return True
