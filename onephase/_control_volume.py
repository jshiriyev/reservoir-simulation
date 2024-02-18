import matplotlib.pyplot as plt

import numpy

from scipy.sparse import linalg

from scipy.sparse import identity
from scipy.sparse import diags
from scipy.sparse import csr_matrix as csr

# It should include Slightly Compressible and Compressible Flows

class Tclass():

    gravity = 9.81

    def __init__(self,grid):

        self.grid = grid

        self.set_delta()
        self.set_perm()
        self.set_static()

    def set_delta(self):
        """delta has the shape of (number of grid, flow dimension x 2),
        and the columns are dx_m, dx_p, dy_m, dy_p, dz_m, dz_p."""

        self._delta = numpy.zeros((self.grid.numtot,self.grid.dimension*2))

        # self._delta = (self.grid._size+self.grid._size[self.grid.index[:,1:],:])/2

        self._delta[:,0] = (self.grid._size[:,0]+self.grid._size[self.grid.index[:,1],0])/2
        self._delta[:,1] = (self.grid._size[:,0]+self.grid._size[self.grid.index[:,2],0])/2

        if self.grid.dimension>1:
            self._delta[:,2] = (self.grid._size[:,1]+self.grid._size[self.grid.index[:,3],1])/2
            self._delta[:,3] = (self.grid._size[:,1]+self.grid._size[self.grid.index[:,4],1])/2

        if self.grid.dimension>2:
            self._delta[:,4] = (self.grid._size[:,2]+self.grid._size[self.grid.index[:,5],2])/2
            self._delta[:,5] = (self.grid._size[:,2]+self.grid._size[self.grid.index[:,6],2])/2

    def set_perm(self):
        """perm has the shape of (number of grid, flow dimension x 2),
        and the columns are kx_m, kx_p, ky_m, ky_p, kz_m, kz_p."""

        self._perm = numpy.zeros((self.grid.numtot,self.grid.dimension*2))

        self._perm[:,0] = (2*self._delta[:,0])/(self.grid._size[:,0]/self.grid._perm[:,0]+
            self.grid._size[self.grid.index[:,1],0]/self.grid._perm[self.grid.index[:,1],0])
        self._perm[:,1] = (2*self._delta[:,1])/(self.grid._size[:,0]/self.grid._perm[:,0]+
            self.grid._size[self.grid.index[:,2],0]/self.grid._perm[self.grid.index[:,2],0])

        if self.grid.dimension>1:
            self._perm[:,2] = (2*self._delta[:,2])/(self.grid._size[:,1]/self.grid._perm[:,1]+
                self.grid._size[self.grid.index[:,3],1]/self.grid._perm[self.grid.index[:,3],1])
            self._perm[:,3] = (2*self._delta[:,3])/(self.grid._size[:,1]/self.grid._perm[:,1]+
                self.grid._size[self.grid.index[:,4],1]/self.grid._perm[self.grid.index[:,4],1])

        if self.grid.dimension>2:
            self._perm[:,4] = (2*self._delta[:,4])/(self.grid._size[:,2]/self.grid._perm[:,2]+
                self.grid._size[self.grid.index[:,5],2]/self.grid._perm[self.grid.index[:,5],2])
            self._perm[:,5] = (2*self._delta[:,5])/(self.grid._size[:,2]/self.grid._perm[:,2]+
                self.grid._size[self.grid.index[:,6],2]/self.grid._perm[self.grid.index[:,6],2])

    def self_static(self):
        """static part of the transmissibility values,
        has the shape of (number of grids, flow dimension x 2)."""
        self._static = (self._perm*self.grid._area)/(self._delta)

    def array(self,fluid):
        """transmissibility arrays of the size (number of grids, flow dimension x 2)"""
        return self._static/fluid._viscosity

    def Tmatrix(self,array):

        tmatrix = csr(self.shape)

        tmatrix = self.utility(tmatrix,0,array)
        tmatrix = self.utility(tmatrix,1,array)

        if self.grid.dimension>1:
            tmatrix = self.utility(tmatrix,2,array)
            tmatrix = self.utility(Tmatrix,3,array)

        if self.grid.dimension>2:
            tmatrix = self.utility(tmatrix,4,array)
            tmatrix = self.utility(tmatrix,5,array)

        return tmatrix

    def Jmatrix(self,array,*columns):
        """
        array   : transmissibility array of size (number of grids, flow dimension x 2)
        columns : integer indicating the direction where the
            constant pressure boundary condition is implemented:

            0   : xmin direction
            1   : xmax direction
            2   : ymin direction
            3   : ymax direction
            4   : zmin direction
            5   : zmax direction

        Return J matrix filled with 2 x transmissibility values.

        """

        jmatrix = csr(self.shape)

        for column in columns:

            fringes = (self.grid.index[:,0]==self.grid.index[:,column+1])

            indices = (self.grid.index[fringes,0],self.grid.index[fringes,0])

            jmatrix += csr((2*array[fringes,column],indices),shape=self.shape)

        return jmatrix

    def Amatrix(self,tstep):

        pore_volume = self.grid._volume*self.grid._poro

        return diags(pore_volume.flatten()/tstep)

    def Qvector(self,array,*pconst):
        """
        array    : transmissibility array of size (number of grids, flow dimension x 2)
        pconst   : constant pressure boundary condition tuple (column,pressure)
        columns  : integer indicating the direction where the
            constant pressure boundary condition is implemented:

            0    : xmin direction
            1    : xmax direction
            2    : ymin direction
            3    : ymax direction
            4    : zmin direction
            5    : zmax direction

        pressure : the value of the constant pressure in psi.

        Return Q vector filled with 2 x transmissibility x pressure values.

        """

        shape = (self.grid.numtot,1)

        qvector = csr(shape)

        for column,pressure in pconst:

            fringes = (self.grid.index[:,0]==self.grid.index[:,column+1])

            indices = (self.grid.index[fringes,0],numpy.zeros(fringes.size))

            qvector += csr((2*array[fringes,column]*pressure*6894.76,indices),shape=shape)

        return qvector

    def Gvector(self,array,fluid):

        return fluid.density*self.gravity*self.Tmatrix(array).dot(self.grid._depth)

    @property
    def shape(self):
        return (self.grid.numtot,self.grid.numtot)
    
    def utility(self,tmatrix:csr,column:int,array:numpy.ndarray):
        """
        Returns updated transmissibility matrix:

        tmatrix : csr_matrix object defined for transmissibility matrix
 
        column  : integer indicating the direction:
            0   : xmin direction
            1   : xmax direction
            2   : ymin direction
            3   : ymax direction
            4   : zmin direction
            5   : zmax direction

        array   : transmissibility array with the size of
                  (number of grids, 2*dimension) and float type

        The transmissibility matrix is updated for:

        1. transmissibility matrix diagonal entries 
        2. transmissibility matrix offset entries

        """
        notOnBorder = ~(self.grid.index[:,0]==self.grid.index[:,column+1])

        # indices of grids not located at the border for the given direction
        diag_indices = (self.grid.index[notOnBorder,0],self.grid.index[notOnBorder,0])

        # indices of neighbor grids in that direction
        offs_indices = (self.grid.index[notOnBorder,0],self.grid.index[notOnBorder,column+1])

        tmatrix += csr((array[notOnBorder,column],diag_indices),shape=tmatrix.shape)
        tmatrix -= csr((array[notOnBorder,column],offs_indices),shape=tmatrix.shape)

        return tmatrix

class MixedSolver():
    """
    This class solves for single phase reservoir flow in Rectangular Grid.
    """

    def __init__(self,grid,fluid,well=None,theta=0):
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

        self.grid  = grid
        self.fluid = fluid
        self.well  = well

        self.theta = theta

        self.tclass = Tclass(self.grid)

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

    def initialize(self,pressure,Swirr=0,ctotal=None):
        """Initializing the reservoir pressure
        
        pressure    : initial pressure in psi
        ctotal      : total compressibility 1/psi
        """
        self._pinit = numpy.ones((self.grid.numtot,1))*pressure*6894.76

        # self.Sw = Swirr
        # self.So = 1-self.Sw

        # coSo = self.Fluids.compressibility[0]*self.So
        # cwSw = self.Fluids.compressibility[1]*self.Sw
        # ctotal = coSo+cwSw+self.PorRock.compressibility

        self._ctotal = ctotal/6894.76

        self._pressure = numpy.zeros((self.grid.numtot,self._time.size))

        self._pressure[:,0] = self._pinit.flatten()

    def solve(self,timesteps):

        P = self._pinit

        LHS = self._Act+self._T+self._J
        
        Psol = np.zeros((self.grid.numtot,timesteps))
        
        for k in range(timesteps):

            # array = self.tclass.array(self.fluid)
            
            RHS = csr.dot(self._Act,P)+self._Q
            
            P = linalg.spsolve(LHS,RHS)
            
            Psol[:,k] = P/6894.76

            P = P.reshape((-1,1))

        # mshape = (grid.numtot,grid.numtot)
        # vshape = (grid.numtot,1)

        # vzeros = np.zeros(self.Wells.number,dtype=int)

        # indices = self.PorRock.grid_indices

        # shape = (grid.numtot,self.grid.numtot)

        # oneperdtime = np.ones(grid.numtot)/self.time_step
        
        # t_correction = csr((oneperdtime,(indices[:,0],indices[:,0])),shape=shape)

        # J_v = (self.JR)/(self.Fluids.viscosity[0]) #*self.Fluids.fvf[0]

        # self.J = csr((J_v[self.well_bhpflags],(self.well_grid[self.well_bhpflags],self.well_grid[self.well_bhpflags])),shape=mshape)
        
        # self.Q = csr(vshape)

        # q_cp = self.well_limits[self.well_bhpflags]*J_v[self.well_bhpflags]

        # q_cr = self.well_limits[~self.well_bhpflags]

        # self.Q += csr((q_cp,(self.well_grid[self.well_bhpflags],vzeros[self.well_bhpflags])),shape=vshape)

        # self.Q += csr((q_cr,(self.well_grid[~self.well_bhpflags],vzeros[~self.well_bhpflags])),shape=vshape)

        # self.Q = self.Q.toarray().flatten()

        # G = self.Gmatrix-t_correction

        # CCC = self.Q/grid.volume/self.PorRock.porosity/self.compressibility
        
        # for i in range(1,self.times.size):
            
        #     b = -self.pressure[:,i-1]/self.time_step+self.b_correction-CCC
            
        #     self.pressure[:,i] = linalg.spsolve(G,b)

        return Psol

    def postprocess(self):

        Y = int((self.grid.numtot-1)/2)

        Pwf = self.pressure[Y,:]+self.Q[Y]/self.JR*self.Fluids.viscosity[0]

        return Pwf

    @property
    def shape(self):
        return (self.grid.numtot,self.grid.numtot)

    @property
    def pinit(self):
        return self._pinit/6894.76

    @property
    def ctotal(self):
        return self._ctotal*6894.76

    @property
    def dtime(self):
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
    def T(self):
        return self._T*24*60*60*6894.76*3.28084**3
    
    @property
    def J(self):
        return self._J*24*60*60*6894.76*3.28084**3

    @property
    def A(self):
        return self._A*24*60*60*3.28084**3

    @property
    def Q(self):
        return self._Q*24*60*60*3.28084**3

    @property
    def G(self):
        return self._G
    
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

    import unittest

    from tests import test_porous_media

    unittest.main(test_porous_media)