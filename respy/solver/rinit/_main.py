import sys

if __name__ == "__main__":
    # sys.path.append(r'C:\Users\javid.shiriyev\Documents\respy')
    sys.path.append(r'C:\Users\3876yl\Documents\respy')

from respy.rinit._time import Time

from respy.solver._onephase import OnePhase

class Main():
    """
    The class initialized reservoir flow in Rectangular Cuboids;
    """

    def __init__(self,grid,depth=None,tcomp=None):
        """
        grid   : It is a GridDelta instance.
        
        depth  : reservoir depth, ft
        tcomp  : total compressibility, 1/psi

        """

        self.grid = grid

        self._depth = depth*0.3048
        self._tcomp = tcomp/6894.76

    def set_time(self,*args,**kwargs):
        
        self.time = Time(*args,**kwargs)

    def set_press(self,pzero=None,refp=None,grad=None):
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

    def set_satur(self):

        pass

    def __call__(self,rrock,fluid,wconds=None,bconds=None):

        return OnePhase(self.grid,rrock,fluid,wconds,bconds)

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
