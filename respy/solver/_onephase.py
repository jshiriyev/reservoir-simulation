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

    def __init__(self,grid,rrock,fluid,wconds=None,bconds=None):
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

        self.set_block(grid)
        self.set_build(grid)

        self.__rrock  = rrock
        self.__fluid  = fluid

        self.__wconds = () if wconds is None else wconds
        self.__bconds = () if bconds is None else bconds

    def set_block(self,grid,**kwargs):

        self.block = Block(grid,depth=depth,tcomp=tcomp)

    def set_build(self,grid):

        self.build = Build(grid)

    def __call__(self,press=None,tcurr=0,tstep=None):

        if press is None:
            press = self._pzero

        rrock  = self.rrock(press)
        fluid  = self.fluid(press)

        wconds = self.wconds(tcurr)
        bconds = self.bconds(tcurr)

        if tstep is None:
            tstep = self.time._steps[0]

        vec = self.block(
            rrock,fluid,wconds,bconds,tstep)

        mat = self.build(vec)

        return mat

    def rrock(self,press):

        if self.rstat and press is None:
            return self.__rrock

        if self.rstat:
            self.__rrock._press = press
            return self.__rrock

        return self.__rrock(press)

    def fluid(self,press):

        if self.fstat and press is None:
            return self.__fluid

        if self.fstat:
            self.__fluid._press = press
            return self.__fluid
        
        return self.__fluid(press)

    def wconds(self,tcurr):

        return [cond for cond in self.__wconds if self.islive(cond,tcurr)]

    def bconds(self,tcurr):

        return [cond for cond in self.__bconds if self.islive(cond,tcurr)]

    def solve(self,**kwargs):
        """Solves the linear system of equation"""

        Pn = numpy.copy(self._press[:,0])
        
        for index,tcurr,tstep in self.time:

            if index==0 or not self.tstat:
                mat = self(Pn,tcurr,tstep)

            if not self.tstat:
                mat = self.iterate(mat,**kwargs)

            Pn = mat.imppress(Pn)

            print(f"{index:10}",Pn.flatten())
            
            self._press[:,index+1] = Pn

            Pn = Pn.reshape((-1,1))

    def iterate(self,mat:Matrix,jacobian=None,maxiter=100,tol=1e-6):

        Pn = numpy.copy(mat._P)

        Pk = numpy.copy(Pn)

        for k in range(maxiter):

            Rv = mat.impresid(Pn)

            if jacobian is None:
                Pk = mat.imppress(Pn)
            elif jacobian is True:
                Jm = self.jacobian(Pk,tcurr,tstep,Pn)
                Pk = Pk+np.linalg.solve(Jm,-Rv)
            else:
                Jm = jacobian(Pk,tstep,tcomp)
                Pk = Pk+np.linalg.solve(Jm,-Rv)
                
            error = np.linalg.norm(Rv,2)

            print(f"{k:2}",f"{error:.5e}",Pk.flatten())

            if np.linalg.norm(Rv,2)<tol:
                break

            mat = self(Pk,tcurr,tstep)

        else:

            print(f"It could not converge after {maxiter} iterations.")

        return mat

    def jacobian(self,press,tcurr,tstep,pprev,delta=10):

        mat0 = self(press,tcurr,tstep)

        resid0 = mat.impresid(pprev)

        press1 = numpy.zeros(press.shape)
        press2 = numpy.zeros(press.shape)

        press1[0::2] = press[0::2]+delta
        press2[1::2] = press[1::2]+delta

        mat1 = self(press1,tcurr,tstep)
        mat2 = self(press2,tcurr,tstep)

        resid1 = mat1.impresid(pprev)
        resid2 = mat2.impresid(pprev)

        vector = numpy.zeros(press.shape)

        vector[0::2] = (resid1-resid0)/delta
        vector[1::2] = (resid2-resid0)/delta

        return self.build.get_diag(vector)

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
