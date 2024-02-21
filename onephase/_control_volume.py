from dataclasses import dataclass,field

import matplotlib.pyplot as plt

import numpy

from scipy.sparse import linalg

from scipy.sparse import identity
from scipy.sparse import diags
from scipy.sparse import csr_matrix as csr

# It should include Slightly Compressible and Compressible Flows

class Tclass():

    def __init__(self,grid,*pconst):
        """
        grid    : RecCuboid instance, rectangular cuboid grids
        pconst  : constant pressure boundary condition tuple (face,pressure)
        
        face    : string indicating the direction where the constant pressure boundary
                  condition is implemented: (xmin,xmax,ymin,ymax,zmin,zmax)
        """

        self.grid = grid

        self.set_static(*pconst)

    def set_static(self,*pconst):
        """Self assigns static transmissibility values."""

        self._staticx = self.set_stat_innface('x')
        self._staticy = self.set_stat_innface('y')
        self._staticz = self.set_stat_innface('z')

        self._staticb = []

        for face, _ in pconst:

            statics = self.set_stat_surface(
                getattr(self.grid,face)._dims,
                getattr(self.grid,face)._area,
                getattr(self.grid,face)._perm)

            self._staticb.append(statics)

    def Tmatrix(self,fluid):
        """Returns transmissibility matrix with dynamic transmissibility values."""

        tmatrix = csr(self.matrix)

        tmatrix = self.tcharge('x',fluid,tmatrix)

        if self.grid.flodim>1:
            tmatrix = self.tcharge('y',fluid,tmatrix)

        if self.grid.flodim>2:
            tmatrix = self.tcharge('z',fluid,tmatrix)

        return tmatrix

    def Jmatrix(self,fluid,*pconst):
        """
        pconst   : constant pressure boundary condition tuple (face,pressure)
        
        face     : string indicating the direction where the constant pressure boundary
                   condition is implemented: (xmin,xmax,ymin,ymax,zmin,zmax)

        Return J matrix filled with constant pressure boundary indices.

        """

        jmatrix = csr(self.matrix)

        for index,(face, _) in enumerate(pconst):

            diag = getattr(self.grid,face)

            vals = 2*self._staticb[index]/fluid._viscosity

            jmatrix += csr((vals,(diag,diag)),shape=self.matrix)

        return jmatrix

    def Amatrix(self,tstep):

        pore_volume = self.grid._volume*self.grid._poro

        return diags(pore_volume.flatten()/tstep)

    def Qvector(self,fluid,*pconst):
        """
        pconst   : constant pressure boundary condition tuple (face,pressure)

        face     : string indicating the direction where the constant pressure boundary
                   condition is implemented: (xmin,xmax,ymin,ymax,zmin,zmax)

        pressure : the value of the constant pressure in psi.

        Return Q vector filled with constant pressure values.

        """

        qvector = csr(self.vector)

        for index,(face,pressure) in enumerate(pconst):

            diag = getattr(self.grid,face)

            vals = 2*self._staticb[index]/fluid._viscosity*pressure*6894.76

            qvector += csr((vals,(diag,numpy.zeros(diag.sum()))),shape=self.vector)

        return qvector

    def Gvector(self,fluid,Tmatrix):

        return fluid._density*self._gravity*Tmatrix.dot(self.grid._depth)

    @property
    def _gravity(self):
        """Gravitational acceleration in SI units"""
        return 9.807

    @property
    def vector(self):
        """Shape of the vectors of transmissibility calculations"""
        return (self.grid.numtot,1)

    @property
    def matrix(self):
        """Shape of the matrices of transmissibility calculations"""
        return (self.grid.numtot,self.grid.numtot)
    
    @property
    def field2si(self):
        """Conversion factor for transmissibility value,
        from Oil Field Units to SI Units."""
        return 1/(3.28084**3*(24*60*60)*6894.76)

    @property
    def si2field(self):
        """Conversion factor for transmissibility value,
        from SI Units to Oil Field Units."""
        return (3.28084**3*(24*60*60)*6894.76)

    def tcharge(self,axis:str,fluid,tmatrix:csr):
        """
        Returns updated transmissibility matrix:

        axis    : x, y, or z

        tmatrix : csr_matrix object defined for transmissibility matrix

        The transmissibility matrix is updated for diagonal and offset entries.
        """

        stat = getattr(self,f"_static{axis}")

        vals = stat/fluid._viscosity

        dneg = getattr(self.grid,f"{axis}neg")
        dpos = getattr(self.grid,f"{axis}pos")

        tmatrix += csr((vals,(dneg,dneg)),shape=tmatrix.shape)
        tmatrix -= csr((vals,(dneg,dpos)),shape=tmatrix.shape)

        tmatrix += csr((vals,(dpos,dpos)),shape=tmatrix.shape)
        tmatrix -= csr((vals,(dpos,dneg)),shape=tmatrix.shape)

        return tmatrix

    def set_stat_innface(self,axis):
        """Returns static transmissibility values for the given
        direction and inner faces (interfaces)."""
        dims = self.stat_diff(
            getattr(self.grid,f"{axis}neg")._dims,
            getattr(self.grid,f"{axis}pos")._dims)

        area = self.stat_area(
            getattr(self.grid,f"{axis}neg")._area,
            getattr(self.grid,f"{axis}pos")._area)
        
        perm = self.stat_perm(
            getattr(self.grid,f"{axis}neg")._dims,
            getattr(self.grid,f"{axis}pos")._dims,
            getattr(self.grid,f"{axis}neg")._perm,
            getattr(self.grid,f"{axis}pos")._perm)

        return self.stat_stat(dims,area,perm)

    def set_stat_surface(self,face):
        """Returns static transmissibility values for the given
        surface on the exterior boundary."""
        
        dims = getattr(self.grid,face)._dims
        area = getattr(self.grid,face)._area
        perm = getattr(self.grid,face)._perm

        return self.stat_stat(dims,area,perm)

    @staticmethod
    def stat_diff(dims_neg,dims_pos):
        return (dims_neg+dims_pos)/2

    @staticmethod
    def stat_area(area_neg,area_pos):
        return (area_neg+area_pos)/2

    @staticmethod
    def stat_perm(dims_neg,dims_pos,perm_neg,perm_pos):
        return (dims_neg+dims_pos)/(dims_neg/perm_neg+dims_pos/perm_pos)

    @staticmethod
    def stat_stat(dims,area,perm):
        return (perm*area)/dims

@dataclass(frozen=True)
class WellCondition:        # bottom hole conditions
    time    : float         # the time for implementing the condition, days
    status  : str           # status of the well, open or shut
    bhp     : float = None  # constant bottom hole pressure if that is the control
    orate   : float = None  # constant oil rate if that is the control
    wrate   : float = None  # constant water rate if that is the control
    grate   : float = None  # constant gas rate if that is the control

@dataclass(frozen=True)
class Well:                 # It is a well dictionary used in the simulator
    name    : str           # name of the well
    block   : tuple         # block indices containing the well 
    sort    : str           # vertical or horizontal

    radius  : float         # well radius, ft
    skin    : float = 0     # skin factor of the well, dimensionless

    conds   : list  = field(default_factory=list)

    def __post_init__(self):
        object.__setattr__(self,'_radius',self.radius*0.3048)

    def add(self,*args,**kwargs):
        """Adds a bottom hole condition to the conds property"""
        self.conds.append(WellCondition(*args,**kwargs)) 

class MixedSolver(Tclass):
    """
    This class solves for single phase reservoir flow in Rectangular Cuboids.
    """

    def __init__(self,grid,fluid,*pconst,well=None,theta=0):
        """
        grid  : It is a RecCuboid (rectangular cuboid) object.

        fluid : There is only one mobile phase in the system. There can be
            two slightly compressible fluids where the second one is at irreducible
            saturation, not mobile.
        
        well  : well schedule object.

        theta : solution type determined in the mixed solver
            when theta = 0, it reduces to Implicit method
            when theta = 1, it reduces to Explicit method
            when theta = 1/2, it is Crank-Nicolson method

        """

        super().__init__(grid,*pconst)

        self.fluid = fluid

        self.pconst = pconst

        self.well  = well

        self.theta = theta

    def set_time(self,tstep:float,total:float=None,nstep:int=1):
        """
        tstep   : time step defined in days
        total   : total simulation time defined in days
        nstep   : number of time steps
        """

        self._tstep = tstep*24*60*60

        if total is None:
            self._ttime = self._tstep*nstep
            self._nstep = nstep
        else:
            self._ttime = total*24*60*60
            self._nstep = int(self._ttime/self._tstep)

        self._time = numpy.arange(
            self._tstep,self._ttime+self._tstep/2,self._tstep)

    def initialize(self,pref,ctotal):
        """Initializing the reservoir pressure
        
        pref    : reference depth and pressure,
                  tuple of (depth in ft, pressure in psi)

        ctotal  : total compressibility 1/psi
        """

        depth,pressure = pref

        depth_diff = (self.grid._depth-depth*0.3048)

        hydro_diff = self.fluid._density*self._gravity*depth_diff

        self._pinit = pressure*6894.76+hydro_diff

        self._ctotal = ctotal/6894.76

        self._pressure = numpy.zeros((self.grid.numtot,self._time.size))

    def solve(self):

        P = self._pinit

        T = self.Tmatrix(self.fluid)
        J = self.Jmatrix(self.fluid,*self.pconst)
        A = self.Amatrix(self._tstep)
        Q = self.Qvector(self.fluid,*self.pconst)
        G = self.Gvector(self.fluid,T)
        
        for k in range(self.nstep):

            if self.theta==0:
                P = self.implicit(P,T,J,A,Q,G)
            elif self.theta==1:
                P = self.explicit(P,T,J,A,Q,G)
            else:
                P = self.mixed(P,T,J,A,Q,G)

            print(f"{k} time step is complete...")
            
            self._pressure[:,k] = P

            P = P.reshape((-1,1))

    def implicit(self,P,T,J,A,Q,G):

        Act = A*self._ctotal

        return linalg.spsolve(T+J+Act,csr.dot(Act,P)+Q+G)

    def mixed(self,P,T,J,A,Q,G):

        Act = A*self._ctotal

        LHS = (1-self.theta)(T+J)+Act
        RHS = csr.dot(Act-self.theta*(T+J),P)+Q+G

        return linalg.spsolve(LHS,RHS)

    def explicit(self,P,T,J,A,Q,G):

        Act = A*self._ctotal

        return P+linalg.spsolve(Act,csr.dot(-(T+J),P)+Q+G)

    def postprocess(self):

        Y = int((self.grid.numtot-1)/2)

        Pwf = self.pressure[Y,:]+self.Q[Y]/self.JR*self.Fluids.viscosity[0]

        return Pwf

    @property
    def pinit(self):
        return self._pinit/6894.76

    @property
    def ctotal(self):
        return self._ctotal*6894.76

    @property
    def tstep(self):
        return self._tstep/(24*60*60)

    @property
    def ttime(self):
        return self._ttime/(24*60*60)

    @property
    def nstep(self):
        return self._nstep
    
    @property
    def time(self):
        return self._time/(24*60*60)

    @property
    def pressure(self):
        return self._pressure/6894.76

    @property
    def pa2psi(self):
        return 1/6894.76

    @property
    def psi2pa(self):
        return 6894.76
    
def newton_solver(grid,timestep,timesteps,T,J,Q):

    array = np.zeros((grid.numtot,timesteps))

    P = grid.pressure_initial

    for j in range(timesteps):

        Pk = P.copy()

        error = 1

        k = 1

        while error>1e-6:

            A = np.diag(100/Pk.flatten())

            F = -np.matmul(T+J+A,Pk)+np.matmul(A,P)+Q

            error = np.linalg.norm(F.flatten(),2)

            JACOB = -(T+J)+np.matmul(-A,np.diag(P.flatten()/Pk.flatten()))

            Pk += np.linalg.solve(JACOB,-F)

            print(f"iteration #{k}: {error = }")
            print(f"{Pk}\n")#{F}\n

            k += 1

        P = Pk.copy()

        array[:,j] = P.flatten()

    return array

def picard_solver(grid,timestep,timesteps,T,J,Q):

    array = np.zeros((grid.numtot,timesteps))

    P = grid.pressure_initial

    for j in range(timesteps):

        Pk = P.copy()

        firstIteration = True

        k = 1

        while firstIteration or error>1e-6:

            firstIteration = False

            A = np.eye(grid.numtot)*100/Pk

            D = T+J+A

            V = np.matmul(A,P)+Q

            F = -np.matmul(D,Pk)+V

            error = np.linalg.norm(F,2)

            Pk = np.linalg.solve(D,V)

            print(f"iteration #{k}: {error=}")
            print(f"{Pk}\n")

            k += 1

        P = Pk.copy()

        array[:,j] = P.flatten()

    return array

if __name__ == "__main__":

    # well = Well('RS',5,'vertical',0.3)

    # well.add(5,'open',orate=127)

    # print(well.conds)

    import unittest

    from tests import test_porous_media

    unittest.main(test_porous_media)