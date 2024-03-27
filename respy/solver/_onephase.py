import sys

if __name__ == "__main__":
    # sys.path.append(r'C:\Users\javid.shiriyev\Documents\respy')
    sys.path.append(r'C:\Users\3876yl\Documents\respy')


from respy.solver._time import Time

from respy.solver._block import Block
from respy.solver._build import Build

from respy.solver._vector import Vector
from respy.solver._matrix import Matrix

class OnePhase():
    """
    The class solves for single phase reservoir flow in Rectangular Cuboids;
    """

    def __init__(self,grid,rrock,fluid,wconds=None,bconds=None,depth=None,tcomp=None):
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

        depth  : depth of reservoir rock, ft

        """

        self.rrock = rrock
        self.fluid = fluid

        wconds = () if wconds is None else wconds
        bconds = () if bconds is None else bconds

        self.block = Block(grid)

        self.build = Build(
              grid = grid,
            wconds = wconds,
            bconds = bconds,
            )

        self._depth = self.set_prop(depth,0.3048)

    def __call__(self,press,tstep,tcomp=None):

        # updating rock and fluid properties with input pressure values
        if callable(self.rrock):
            rrock = self.rrock(press)
        else:
            rrock = self.rrock
            rrock._press = press

        if callable(self.fluid):
            fluid = self.fluid(press)
        else:
            fluid = self.fluid
            fluid._press = press

        vec = self.block(rrock,fluid)

        vec.set_A(tstep)
        vec.set_C(tcomp)

        mat = self.build(rrock,fluid,vec)

        return mat

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

        Pn = numpy.copy(self._pzero)

        self._press = numpy.zeros(self.shape)
        
        for n,step in enumerate(self.time.steps):

            if n==0:
                mat = self.get_matrix(Pn,step)

            if not self.static:
                mat = self.get_matrix(Pn,step)

            Pn = mat.implicit_pressure(Pn)

            print(f"{n:10}",Pn.flatten())
            
            self._press[:,n] = Pn

            Pn = Pn.reshape((-1,1))

    def iterate(self,mat:Matrix,jacobian=None,maxiter=100,tol=1e-6):

        Pn = numpy.copy(mat._P)

        Pk = numpy.copy(Pn)

        for k in range(maxiter):

            Rv = mat.implicit_residual(Pn)

            if jacobian is None:
                Pk = mat.implicit_pressure(Pn)
            else:
                Jm = jacobian(Pk,tstep,tcomp)
                Pk = Pk+np.linalg.solve(Jm,-Rv)
                
            error = np.linalg.norm(Rv,2)

            print(f"{k:2}",f"{error:.5e}",Pk.flatten())

            if np.linalg.norm(Rv,2)<tol:
                break

            mat = self.get_matrix(Pk,tstep,tcomp)

        else:

            print(f"It could not converge after {maxiter} iterations.")

        return mat

    @property
    def static(self):
        return all((callable(self.rrock),callable(self.fluid)))

    @property
    def depth(self):
        if self._depth is not None:
            return self._depth/0.3048

    @property
    def pzero(self):
        return self._pzero/6894.75729

    @property
    def tcomp(self):
        return self._tcomp*6894.75729

    @property
    def press(self):
        return self._press/6894.75729

    @property
    def shape(self):
        return (self.block.cube.numtot,self.time.nums)

    @staticmethod
    def set_prop(prop,conv=1.):
        if prop is not None:
            return numpy.asarray(prop).astype(numpy.float_)*conv