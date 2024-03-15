import sys

if __name__ == "__main__":
    sys.path.append(r'C:\Users\javid.shiriyev\Documents\respy')

from respy.solver._cube import RecCube

from respy.solver._time import Time

from respy.solver._block import Block

from respy.solver._assemble import Assemble

class OnePhase():
    """
    The class solves for single phase reservoir flow in Rectangular Cuboids;
    """

    def __init__(self,grid,rrock,fluid,wconds=None,bconds=None):
        """
        grid   : It is a RecCuboid (rectangular cuboid) object.

        rrock  : It is a class with reservoir rock properties.

        fluid  : There is only one mobile phase in the system. There
                 can be two slightly compressible fluids where the
                 second one is at irreducible saturation, not mobile.

        wconds : tuple of WellCond item.

        bconds : tuple of BoundCond item.

        """

        self.grid = grid

        self.set_grid(rrock,fluid)

        self.wconds = () if wconds is None else wconds
        self.bconds = () if bconds is None else bconds

    def set_grid(self):

        self.cube = RecCube(
             edge = self.grid.edge,
             plat = self.grid.plat,
            rrock = rrock,
            fluid = fluid,
            )

    def set_rrock(self):

        props = ("poro","xperm","yperm","zperm","comp","depth")

        count = 0

        for prop in props:

            pp = getattr(self.rrock,prop)

            if callable(pp):
                count += 1

            setattr(self.cube,prop,pp)

        if count>0:
            return True

        return False

    def set_fluid(self):

        props = ("visc","rho","comp","fvf","rel0")

        count = 0

        for prop in props:

            pp = getattr(self.fluid,prop)

            if callable(pp):
                count += 1

            setattr(self.cube,prop,pp)

        if count>0:
            return True

        return False

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

    def solve(self,**kwargs):
        """
        theta  : solution type determined in the mixed solver
            when theta = 0, it reduces to Implicit method
            when theta = 1, it reduces to Explicit method
            when theta = 1/2, it is Crank-Nicolson method
        """

        pass

    def static(self):

        Pn = self._pinit

        vec = Block(cube,Pn,tstep,comp)
        mat = Shape(cube,Pn,vec)

        self._pressure = numpy.zeros((self.grid.numtot,self._time.size))
        
        for n in range(self.nstep):

            Pn = mat.implicit_pressure(Pn)

            print(f"{n:10}",Pn.flatten())
            
            self._pressure[:,n] = Pn

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
    def pressure(self):
        return self._pressure/6894.75729

    def postprocess(self):

        Y = int((self.grid.numtot-1)/2)

        Pwf = self.pressure[Y,:]+self.Q[Y]/self.JR*self.Fluids.viscosity[0]

        return Pwf

    @property
    def pinit(self):
        return self._pinit/6894.75729

    @property
    def comp(self):
        return self._comp*6894.75729
