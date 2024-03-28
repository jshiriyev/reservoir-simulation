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
        tcomp  : total compressibility of rock and fluid system, 1/psi

        """

        kwargs = {}

        if depth is not None:
            kwargs["depth"] = depth*0.3048

        if tcomp is not None:
            kwargs["tcomp"] = tcomp/6894.76

        self.block  = Block(grid,**kwargs)
        self.build  = Build(grid)

        self.rrock  = rrock
        self.fluid  = fluid

        self.wconds = () if wconds is None else wconds
        self.bconds = () if bconds is None else bconds

    def set_pzero(self,pzero=None,refp=None,grad=None):
        """Calculates the initial pressure
        
        pzero   : initial pressure in psi; If not defined, it will be
                  calculated from reference point and reservoir rock depths.

        refp    : reference point (depth:ft,pressure:psi)
        grad    : fluid hydrostatic gradient, psi/ft
        """

        if pzero is None:
            pzero = refp[1]+grad*(self.block.depth-refp[0])

        self._pzero = pzero*6894.76

    def set_time(self,time:Time):
        
        self.time = Time

    def __call__(self,press,tcurr,tstep):

        rrock  = self.get_rrock(press)
        fluid  = self.get_fluid(press)

        wconds = self.get_wconds(tcurr)
        bconds = self.get_bconds(tcurr)

        vec = self.block(
            rrock,fluid,wconds,bconds,tstep)

        mat = self.build(vec)

        return mat

    def get_rrock(self,press):

        if self.rstat:
            return self.rrock(press)
        
        self.rrock._press = press

        return self.rrock

    def get_fluid(self,press):

        if self.fstat:
            return self.fluid(press)

        self.fluid._press = press

        return self.fluid

    def get_wconds(self,tcurr):

        return self.wconds

    def get_bconds(self,tcurr):

        return self.bconds

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
    def pzero(self):
        return self._pzero/6894.75729

    @property
    def press(self):
        return self._press/6894.75729

    @property
    def shape(self):
        return (self.block.nums,self.time.nums)

    @property
    def rstat(self):
        return callable(self.rrock)

    @property
    def fstat(self):
        return callable(self.fluid)