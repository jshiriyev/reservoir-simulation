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
                 calculates fluid properties at any given pressure. In the
                 solution, there can be only one mobile phase. Multiple
                 fluid phases may exist with only one mopbile phase; the
                 rest is at irreducible saturation, not mobile.

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

    def init(self,pinit=None,refp=None,grad=None):
        """Calculates the initial pressure
        
        pinit   : initial pressure in psi; If not defined, will be
                  calculated from reference point and reservoir rock depths.

        refp    : reference point (depth:ft,pressure:psi)
        grad    : fluid hydrostatic gradient, psi/ft
        """

        if pinit is None:
            pinit = refp[1]+grad*(self.rrock.depth-refp[0])

        self._pinit = pinit*6894.75729

    def set_time(self,time:Time):
        
        self.time = Time

    def set_comp(self,comp=None):
        """
        comp  : total compressibility 1/psi
        """

        if comp is None:
            comp = self.fluid.comp+self.rrock.comp

        self._comp = comp/6894.75729

    def solve(self):
        """
        theta  : solution type determined in the mixed solver
            when theta = 0, it reduces to Implicit method
            when theta = 1, it reduces to Explicit method
            when theta = 1/2, it is Crank-Nicolson method
        """

        Pn = self._pinit

        vec = self.block()

        vec.set_A()
        vec.set_C()

        mat = self.build(vec)

        self._press = numpy.zeros((self.grid.numtot,self._time.size))
        
        for n in range(self.nstep):

            Pn = mat.implicit_pressure(Pn)

            print(f"{n:10}",Pn.flatten())
            
            self._press[:,n] = Pn

            Pn = Pn.reshape((-1,1))

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

    @property
    def press(self):
        return self._press/6894.75729

    def pproc(self):

        Y = int((self.grid.numtot-1)/2)

        Pwf = self.press[Y,:]+self.Q[Y]/self.JR*self.Fluids.viscosity[0]

        return Pwf

    @property
    def pinit(self):
        return self._pinit/6894.75729

    @property
    def comp(self):
        return self._comp*6894.75729
