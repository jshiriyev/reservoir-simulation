import sys

if __name__ == "__main__":
    sys.path.append(r'C:\Users\javid.shiriyev\Documents\respy')

from respy.solver._cube import RecCube

from respy.solver._time import Time

from respy.solver._block import Block
from respy.solver._build import Build

class OnePhase():
    """
    The class solves for single phase reservoir flow in Rectangular Cuboids;
    """

    def __init__(self,grid,rrock,fluid,wconds=None,bconds=None):
        """
        grid   : It is a GridDelta instance.

        rrock  : It is a ResRock instance or any other rock class that
                 calculates rock properties at any given pressure.

        fluid  : It is a Fluid instance or any other fluid class that
                 calculates fluid properties at any given pressure. Multiple
                 fluid phases may exist with only one mopbile phase; the
                 rest of phases is at irreducible saturation, not mobile.

        wconds : tuple of WellCond instance.

        bconds : tuple of BoundCond instance.

        """

        self.block = Block(
              cube = RecCube(grid.edge,grid.plat),
             rrock = rrock,
             fluid = fluid,
            wconds = () if wconds is None else wconds,
            bconds = () if bconds is None else bconds,
            )

        self.build = Build(
              cube = self.block.cube,
            wconds = self.block.wconds,
            bconds = self.block.bconds,
            )

    def set_pzero(self,pzero=None,refp=None,grad=None):
        """Calculates the initial pressure
        
        pzero   : initial pressure in psi; If not defined, it will be
                  calculated from reference point and reservoir rock depths.

        refp    : reference point (depth:ft,pressure:psi)
        grad    : fluid hydrostatic gradient, psi/ft
        """

        if pzero is None:
            pzero = refp[1]+grad*(self.block.rrock.depth-refp[0])

        self._pzero = pzero*6894.75729

    def set_time(self,time:Time):
        
        self.time = Time

    def solve(self,theta=0):
        """
        Solves the linear system of equation

        theta   : solution type determined in the mixed solver
                  0 means Implicit method
                  1 means Explicit method
                1/2 means Crank-Nicolson method
        """

        Pn = self._pzero.copy()

        self._press = numpy.zeros((self.block.cube.numtot,self.time.nums))
        
        for n,step in enumerate(self.time.steps):

            if n==0:
                mat = self.get_matrix(Pn,step)

            if not self.static:
                mat = self.get_matrix(Pn,step)

            Pn = mat.implicit_pressure(Pn)

            print(f"{n:10}",Pn.flatten())
            
            self._press[:,n] = Pn

            Pn = Pn.reshape((-1,1))

    def get_matrix(self,press,tstep,tcomp):

        vec = self.block(press)

        vec.set_A(tstep)
        vec.set_C(tcomp)

        return self.build(vec)

    @property
    def static(self):
        return all(self.block.static)

    @property
    def pzero(self):
        return self._pzero/6894.75729

    @property
    def comp(self):
        return self._comp*6894.75729

    @property
    def press(self):
        return self._press/6894.75729



    def pproc(self):

        Y = int((self.grid.numtot-1)/2)

        Pwf = self.press[Y,:]+self.Q[Y]/self.JR*self.Fluids.viscosity[0]

        return Pwf

def picard(self):

    for k in range(100):

        Fm = mat.implicit_residual(tstep,Pn,Pk)

        Pk = mat.implicit_pressure(tstep,Pn,Pk)

        error = np.linalg.norm(Fm,2)

        # print(f"{k:2}",f"{error:.5e}",Pk.flatten())

        if np.linalg.norm(Fm,2)<1e-6:
            break

def newton(self):

    for k in range(100):

        F = mat.implicit_residual(tstep,P,Pk)
        Z = jacobian(tstep,P,Pk)
        Pk += np.linalg.solve(Z,-F)
        error = np.linalg.norm(F,2)
        # print(f"{k:2}",f"{error:.5e}",Pk.flatten())
        if np.linalg.norm(F,2)<1e-6:
            break