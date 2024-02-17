import matplotlib.pyplot as plt

import numpy

from scipy.sparse import linalg

from scipy.sparse import identity
from scipy.sparse import diags
from scipy.sparse import csr_matrix as csr

# It should include Slightly Compressible and Compressible Flows

class Tclass():

    def __init__(self,grid):

        self.grid = grid

    def delta(self):
        """Returned delta has the shape of (number of grid, flow dimension x 2),
        and the columns are dx_m, dx_p, dy_m, dy_p, dz_m, dz_p."""

        grid = self.grid

        _delta = numpy.zeros((grid.numtot,grid.dimension*2))

        _delta[:,0] = (grid._size[:,0]+grid._size[grid.index[:,1],0])/2
        _delta[:,1] = (grid._size[:,0]+grid._size[grid.index[:,2],0])/2

        if grid.dimension>1:
            _delta[:,2] = (grid._size[:,1]+grid._size[grid.index[:,3],1])/2
            _delta[:,3] = (grid._size[:,1]+grid._size[grid.index[:,4],1])/2

        if grid.dimension>2:
            _delta[:,4] = (grid._size[:,2]+grid._size[grid.index[:,5],2])/2
            _delta[:,5] = (grid._size[:,2]+grid._size[grid.index[:,6],2])/2

        return _delta

    def perm(self,delta=None):
        """Returned perm has the shape of (number of grid, flow dimension x 2),
        and the columns are kx_m, kx_p, ky_m, ky_p, kz_m, kz_p."""

        grid = self.grid

        delta = self.delta() if delta is None else delta

        _perm = numpy.zeros((grid.numtot,grid.dimension*2))

        _perm[:,0] = (2*delta[:,0])/(grid._size[:,0]/grid._perm[:,0]+
            grid._size[grid.index[:,1],0]/grid._perm[grid.index[:,1],0])
        _perm[:,1] = (2*delta[:,1])/(grid._size[:,0]/grid._perm[:,0]+
            grid._size[grid.index[:,2],0]/grid._perm[grid.index[:,2],0])

        if grid.dimension>1:
            _perm[:,2] = (2*delta[:,2])/(grid._size[:,1]/grid._perm[:,1]+
                grid._size[grid.index[:,3],1]/grid._perm[grid.index[:,3],1])
            _perm[:,3] = (2*delta[:,3])/(grid._size[:,1]/grid._perm[:,1]+
                grid._size[grid.index[:,4],1]/grid._perm[grid.index[:,4],1])

        if grid.dimension>2:
            _perm[:,4] = (2*delta[:,4])/(grid._size[:,2]/grid._perm[:,2]+
                grid._size[grid.index[:,5],2]/grid._perm[grid.index[:,5],2])
            _perm[:,5] = (2*delta[:,5])/(grid._size[:,2]/grid._perm[:,2]+
                grid._size[grid.index[:,6],2]/grid._perm[grid.index[:,6],2])

        return _perm

    def array(self,fluid,delta=None,perm=None):

        grid = self.grid

        delta = self.delta() if delta is None else delta

        perm = self.perm(delta) if perm is None else perm

        # transmissibility[:,0] = (perm[:,0]*grid._area[:,0])/(fluid._viscosity*dx_m)
        # transmissibility[:,1] = (perm[:,1]*grid._area[:,1])/(fluid._viscosity*dx_p)
        # transmissibility[:,2] = (perm[:,2]*grid._area[:,2])/(fluid._viscosity*dy_m)
        # transmissibility[:,3] = (perm[:,3]*grid._area[:,3])/(fluid._viscosity*dy_p)
        # transmissibility[:,4] = (perm[:,4]*grid._area[:,4])/(fluid._viscosity*dz_m)
        # transmissibility[:,5] = (perm[:,5]*grid._area[:,5])/(fluid._viscosity*dz_p)

        return (perm*grid._area)/(fluid._viscosity*delta) # transmissibility

    def matrix(self,fluid,delta=None,perm=None):

        grid = self.grid

        delta = self.delta() if delta is None else delta
        
        perm = self.perm(delta) if perm is None else perm

        array = self.array(fluid,delta,perm)

        indices = grid.index

        noxmin = ~(indices[:,0]==indices[:,1])
        noxmax = ~(indices[:,0]==indices[:,2])

        if grid.dimension>1:
            noymin = ~(indices[:,0]==indices[:,3])
            noymax = ~(indices[:,0]==indices[:,4])

        if grid.dimension>2:
            nozmin = ~(indices[:,0]==indices[:,5])
            nozmax = ~(indices[:,0]==indices[:,6])

        id_noxmin = indices[noxmin,0] #index of grid not at xmin boundary
        id_noxmax = indices[noxmax,0] #index of grid not at xmax boundary

        if grid.dimension>1:
            id_noymin = indices[noymin,0] #index of grid not at ymin boundary
            id_noymax = indices[noymax,0] #index of grid not at ymax boundary

        if grid.dimension>2:
            id_nozmin = indices[nozmin,0] #index of grid not at zmin boundary
            id_nozmax = indices[nozmax,0] #index of grid not at zmax boundary

        idNnoxmin = indices[noxmin,1] #index of xmin neighbors for id_noxmin grid
        idNnoxmax = indices[noxmax,2] #index of xmax neighbors for id_noxmax grid

        if grid.dimension>1:
            idNnoymin = indices[noymin,3] #index of ymin neighbors for id_noymin grid
            idNnoymax = indices[noymax,4] #index of ymax neighbors for id_noymax grid

        if grid.dimension>2:
            idNnozmin = indices[nozmin,5] #index of zmin neighbors for id_nozmin grid
            idNnozmax = indices[nozmax,6] #index of zmax neighbors for id_nozmax grid

        shape = (grid.numtot,grid.numtot)

        Tmatrix = csr(shape)

        Tmatrix -= csr((array[noxmin,0],(id_noxmin,idNnoxmin)),shape=shape)
        Tmatrix += csr((array[noxmin,0],(id_noxmin,id_noxmin)),shape=shape)
        Tmatrix -= csr((array[noxmax,1],(id_noxmax,idNnoxmax)),shape=shape)
        Tmatrix += csr((array[noxmax,1],(id_noxmax,id_noxmax)),shape=shape)

        if grid.dimension>1:
            Tmatrix -= csr((array[noymin,2],(id_noymin,idNnoymin)),shape=shape)
            Tmatrix += csr((array[noymin,2],(id_noymin,id_noymin)),shape=shape)
            Tmatrix -= csr((array[noymax,3],(id_noymax,idNnoymax)),shape=shape)
            Tmatrix += csr((array[noymax,3],(id_noymax,id_noymax)),shape=shape)
        
        if grid.dimension>2:
            Tmatrix -= csr((array[nozmin,4],(id_nozmin,idNnozmin)),shape=shape)
            Tmatrix += csr((array[nozmin,4],(id_nozmin,id_nozmin)),shape=shape)
            Tmatrix -= csr((array[nozmax,5],(id_nozmax,idNnozmax)),shape=shape)
            Tmatrix += csr((array[nozmax,5],(id_nozmax,id_nozmax)),shape=shape)

        return Tmatrix

class MixedSolver():
    """
    This class solves for single phase reservoir flow in Rectangular Grid.
    """

    def __init__(self,grid,fluid,wells=None,theta=0):
        """
        grid : It is a RecCuboid (rectangular cuboid) object.

        fluid : There is only one mobile phase in the system. There can be
            two slightly compressible fluids where the second one is at irreducible
            saturation, not mobile.
        
        wells : well schedule object.

        theta : solution type determined in the mixed solver
            when theta = 0, it reduces to Implicit method
            when theta = 1, it reduces to Explicit method
            when theta = 1/2, it is Crank-Nicolson method

        """

        self.grid  = grid
        self.fluid = fluid
        self.wells = wells
        self.theta = theta

        self.tclass = Tclass(self.grid)

    def initialize(self,pressure,Swirr=0,ctotal=None):
        """Initializing the reservoir pressure
        
        pressure    : initial pressure in psi
        ctotal      : total compressibility 1/psi
        """
        self._pinit = np.ones((self.grid.numtot,1))*pressure*6894.76

        # self.Sw = Swirr
        # self.So = 1-self.Sw

        # coSo = self.Fluids.compressibility[0]*self.So
        # cwSw = self.Fluids.compressibility[1]*self.Sw
        # ctotal = coSo+cwSw+self.PorRock.compressibility

        self._ctotal = ctotal/6894.76

    def set_time(self,step:float,total:float):
        """
        step    : time step defined in days
        total   : total simulation time defined in days
        """

        self._dtime = step*24*60*60
        self._ttime = total*24*60*60

        self._times = np.arange(self._dtime,self._ttime+self._dtime,self._dtime)

        self._pressure = np.zeros((self.grid.numtot,self._times.size))

        self._pressure[:,0] = self._pinit

    def set_Tmatrix(self):

        self._T = self.tclass.matrix(self.fluid)

    def set_Jmatrix(self):

        coeff = np.asarray(2*self._trans[-1,0]).flatten()
        index = np.asarray(self.grid.numtot-1).flatten()

        J = csr((coeff,(index,index)),shape=self.shape)

        self._J = J

    def set_Amatrix(self):
        
        self._dt = dt*(24*60*60)

        coeff = (self.grid._area[:,0].reshape((-1,1))*self.grid._size*self.grid._poro.reshape((-1,1)))/(self._dt)

        A = csr((coeff[:,0],(self.grid.index[:,0],self.grid.index[:,0])),shape=self.shape)

        self._A = A

        self._Act = self._A*self._ct

    def set_Qvector(self,Pbound):
        """Pbound in psi"""

        shape = (self.grid.numtot,1)

        coeff = np.asarray(2*self._trans[-1,0]*Pbound*6894.76).flatten()
        index = np.asarray(self.grid.numtot-1).flatten()

        Q = csr((coeff,(index,np.asarray(0).flatten())),shape=shape)

        self._Q = Q

    def set_Gvector(self):

        pass

    def solve(self,timesteps):

        P = self._pinit

        LHS = self._Act+self._T+self._J
        
        Psol = np.zeros((self.grid.numtot,timesteps))
        
        for k in range(timesteps):
            
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
        return self._dtime/(24*60*60)

    @property
    def ttime(self):
        return self._ttime/(24*60*60)

    @property
    def times(self):
        return self._times/(24*60*60)

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
    def Act(self):
        return self._Act*24*60*60*6894.76*3.28084**3

    @property
    def Q(self):
        return self._Q*24*60*60*3.28084**3

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