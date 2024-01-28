import matplotlib.pyplot as plt
import numpy as np

from scipy.sparse import csr_matrix as csr
from scipy.sparse import diags

from scipy.sparse.linalg import spsolve as sps

# import fluids

# from borepy.items import PorRock
# from borepy.items import Wells

# It should include Slightly Compressible and Compressible Flows

def TransMatrix1D(grids,TransVector1D):

    T = np.zeros((grids.numtot,grids.numtot))

    for i in range(grids.numtot):

        if i>0:
            T[i,i-1] = -TransVector1D[i,0]
            T[i,i] += TransVector1D[i,0]

        if i<grids.numtot-1:
            T[i,i+1] = -TransVector1D[i,1]
            T[i,i] += TransVector1D[i,1]

    return T

class Implicit1D():

    @staticmethod
    def solve(length,time,pinit,pright,eta,ngrids,nsteps):
        """
        lenght  : length of core sample
        time    : time at which to calculate pressure
        pinit   : initial pressure
        pright  : constant pressure boundary implemented at right hand side
        eta     : hydraulic diffusivity
        ngrids  : number of pressure calculation points
        nsteps  : number of time steps to be used in finite difference solution
        """

        x = np.arange(length/ngrids/2,length,length/ngrids)

        T = CoreFlow.tondimtime(time,eta,length)

        P = CoreFlow.ndimpressure(ngrids,T/nsteps,nsteps)
        
        P = CoreFlow.todimpressure(P,pinit,pright)

        return x,P
    
    @staticmethod
    def tondimx(x,length):
        """Converts dimensional x to non-dimensional x."""

        return x/length

    @staticmethod
    def todimx(x,length):
        """Converts non-dimensional x to dimensional x."""

        return x*length

    @staticmethod
    def tondimtime(time,eta,length):
        """Converts dimensional time to non-dimensional time."""

        return eta*time/length**2

    @staticmethod
    def todimtime(time,eta,length):
        """Converts non-dimensional time to dimensional time."""
        
        return time*length**2/eta

    @staticmethod
    def tondimpressure(pressure,pinit,pright):
        """Converts dimensional pressure to non-dimensional pressure."""

        return (pressure-pright)/(pinit-pright)

    @staticmethod
    def todimpressure(pressure,pinit,pright):
        """Converts non-dimensional pressure to dimensional pressure."""

        return pressure*(pinit-pright)+pright

    @staticmethod
    def ndimpressure(ngrids:int,deltat:int,nsteps:int):
        """Calculates the non-dimensional pressure with finite difference method.

        """

        shape = (ngrids,ngrids)

        deltax = 1/ngrids

        Gmatrix = csr(shape)

        indices = np.array([i for i in range(ngrids-1)],dtype="int16")

        coeffis = [deltat/deltax**2 for _ in range(ngrids-1)]

        Gmatrix -= csr((coeffis,(indices,indices+1)),shape=shape)
        Gmatrix += csr((coeffis,(indices,indices)),shape=shape)
        
        Gmatrix -= csr((coeffis,(indices+1,indices)),shape=shape)
        Gmatrix += csr((coeffis,(indices+1,indices+1)),shape=shape)

        Gmatrix += csr(([2*deltat/deltax**2],([ngrids-1],[ngrids-1])),shape=shape)

        Gmatrix += diags([1 for _ in range(ngrids)],shape=shape)

        pinit = np.ones((ngrids,))

        pn = pinit

        for n in range(nsteps):

            pn = sps(Gmatrix,pn)

        return pn

def TransVector3D(grids):
    
    dx_m = (grids.size[:,0]+grids.size[grids.index[:,1],0])/2
    dx_p = (grids.size[:,0]+grids.size[grids.index[:,2],0])/2
    dy_m = (grids.size[:,1]+grids.size[grids.index[:,3],1])/2
    dy_p = (grids.size[:,1]+grids.size[grids.index[:,4],1])/2
    dz_m = (grids.size[:,2]+grids.size[grids.index[:,5],2])/2
    dz_p = (grids.size[:,2]+grids.size[grids.index[:,6],2])/2

    kx_m = (2*dx_m)/(grids.size[:,0]/grids.permeability[:,0]+
                     grids.size[grids.index[:,1],0]/grids.permeability[grids.index[:,1],0])
    kx_p = (2*dx_p)/(grids.size[:,0]/grids.permeability[:,0]+
                     grids.size[grids.index[:,2],0]/grids.permeability[grids.index[:,2],0])
    ky_m = (2*dy_m)/(grids.size[:,1]/grids.permeability[:,1]+
                     grids.size[grids.index[:,3],1]/grids.permeability[grids.index[:,3],1])
    ky_p = (2*dy_p)/(grids.size[:,1]/grids.permeability[:,1]+
                     grids.size[grids.index[:,4],1]/grids.permeability[grids.index[:,4],1])
    kz_m = (2*dz_m)/(grids.size[:,2]/grids.permeability[:,2]+
                     grids.size[grids.index[:,5],2]/grids.permeability[grids.index[:,5],2])
    kz_p = (2*dz_p)/(grids.size[:,2]/grids.permeability[:,2]+
                     grids.size[grids.index[:,6],2]/grids.permeability[grids.index[:,6],2])

    transmissibility = np.zeros((grids.numtot,6))

    transmissibility[:,0] = (kx_m*grids.area[:,0])/(grids.viscosity*dx_m)
    transmissibility[:,1] = (kx_p*grids.area[:,1])/(grids.viscosity*dx_p)
    transmissibility[:,2] = (ky_m*grids.area[:,2])/(grids.viscosity*dy_m)
    transmissibility[:,3] = (ky_p*grids.area[:,3])/(grids.viscosity*dy_p)
    transmissibility[:,4] = (kz_m*grids.area[:,4])/(grids.viscosity*dz_m)
    transmissibility[:,5] = (kz_p*grids.area[:,5])/(grids.viscosity*dz_p)

    return transmissibility

class Implicit3D():
    
    """
    This class is supposed to generate regular grids in cartesian
    coordinates and perform reservoir simulation.
    """
    
    def __init__(self):
        
        self.PorRock = PorRock(geo="rectangle")(window=None)

        # There can be two slightly compressible fluids where the
        # second one is at irreducible saturation, not mobile

        self.Fluids = Fluids()(number=2)

        self.Wells = Wells()()

    def initialize(self,pressure0,Swirr=0,ctotal=None):

        self.pressure0 = pressure0

        self.Sw = Swirr
        self.So = 1-self.Sw

        if ctotal is None:
            coSo = self.Fluids.compressibility[0]*self.So
            cwSw = self.Fluids.compressibility[1]*self.Sw
            ctotal = coSo+cwSw+self.PorRock.compressibility

        self.compressibility = ctotal

        cons = self.PorRock.porosity*self.Fluids.viscosity[0]*self.compressibility

        self.diffusivity = (self.PorRock.permeability)/(np.reshape(cons,(-1,1)))

    def set_times(self,step,total):

        self.time_step = float(step)
        self.time_total = float(total)

        self.times = np.arange(step,total+step,step)

        self.pressure = np.zeros((self.PorRock.grid_numtot,self.times.size))

        self.pressure[:,0] = self.pressure0

    def set_transmissibility(self):

        sizes = self.PorRock.grid_sizes

        indices = self.PorRock.grid_indices

        permeability = self.PorRock.permeability
        
        dx_m = (sizes[:,0]+sizes[indices[:,1],0])/2
        dx_p = (sizes[:,0]+sizes[indices[:,2],0])/2
        dy_m = (sizes[:,1]+sizes[indices[:,3],1])/2
        dy_p = (sizes[:,1]+sizes[indices[:,4],1])/2
        dz_m = (sizes[:,2]+sizes[indices[:,5],2])/2
        dz_p = (sizes[:,2]+sizes[indices[:,6],2])/2

        kx_m = (2*dx_m)/(sizes[:,0]/permeability[:,0]+
                         sizes[indices[:,1],0]/permeability[indices[:,1],0])
        kx_p = (2*dx_p)/(sizes[:,0]/permeability[:,0]+
                         sizes[indices[:,2],0]/permeability[indices[:,2],0])
        ky_m = (2*dy_m)/(sizes[:,1]/permeability[:,1]+
                         sizes[indices[:,3],1]/permeability[indices[:,3],1])
        ky_p = (2*dy_p)/(sizes[:,1]/permeability[:,1]+
                         sizes[indices[:,4],1]/permeability[indices[:,4],1])
        kz_m = (2*dz_m)/(sizes[:,2]/permeability[:,2]+
                         sizes[indices[:,5],2]/permeability[indices[:,5],2])
        kz_p = (2*dz_p)/(sizes[:,2]/permeability[:,2]+
                         sizes[indices[:,6],2]/permeability[indices[:,6],2])

        porosity = self.PorRock.porosity

        viscosity = self.Fluids.viscosity[0]

        etax_m = kx_m/(porosity*viscosity*self.compressibility)
        etax_p = kx_p/(porosity*viscosity*self.compressibility)
        etay_m = ky_m/(porosity*viscosity*self.compressibility)
        etay_p = ky_p/(porosity*viscosity*self.compressibility)
        etaz_m = kz_m/(porosity*viscosity*self.compressibility)
        etaz_p = kz_p/(porosity*viscosity*self.compressibility)

        self.transmissibility = np.zeros((self.PorRock.grid_numtot,6))

        self.transmissibility[:,0] = (2*etax_m)/(dx_m*(dx_m+dx_p))
        self.transmissibility[:,1] = (2*etax_p)/(dx_p*(dx_m+dx_p))
        self.transmissibility[:,2] = (2*etay_m)/(dy_m*(dy_m+dy_p))
        self.transmissibility[:,3] = (2*etay_p)/(dy_p*(dy_m+dy_p))
        self.transmissibility[:,4] = (2*etaz_m)/(dz_m*(dz_m+dz_p))
        self.transmissibility[:,5] = (2*etaz_p)/(dz_p*(dz_m+dz_p))

    def set_matrix(self,order=2):

        indices = self.PorRock.grid_indices

        noxmin = ~(indices[:,0]==indices[:,1])
        noxmax = ~(indices[:,0]==indices[:,2])
        noymin = ~(indices[:,0]==indices[:,3])
        noymax = ~(indices[:,0]==indices[:,4])
        nozmin = ~(indices[:,0]==indices[:,5])
        nozmax = ~(indices[:,0]==indices[:,6])

        id_noxmin = indices[noxmin,0] #id of grids not at xmin boundary
        id_noxmax = indices[noxmax,0] #id of grids not at xmax boundary
        id_noymin = indices[noymin,0] #id of grids not at ymin boundary
        id_noymax = indices[noymax,0] #id of grids not at ymax boundary
        id_nozmin = indices[nozmin,0] #id of grids not at zmin boundary
        id_nozmax = indices[nozmax,0] #id of grids not at zmax boundary

        idNnoxmin = indices[noxmin,1] #id of xmin neighbors for id_noxmin grids
        idNnoxmax = indices[noxmax,2] #id of xmax neighbors for id_noxmax grids
        idNnoymin = indices[noymin,3] #id of ymin neighbors for id_noymin grids
        idNnoymax = indices[noymax,4] #id of ymax neighbors for id_noymax grids
        idNnozmin = indices[nozmin,5] #id of zmin neighbors for id_nozmin grids
        idNnozmax = indices[nozmax,6] #id of zmax neighbors for id_nozmax grids

        cx_m = self.transmissibility[noxmin,0]
        cx_p = self.transmissibility[noxmax,1]
        cy_m = self.transmissibility[noymin,2]
        cy_p = self.transmissibility[noymax,3]
        cz_m = self.transmissibility[nozmin,4]
        cz_p = self.transmissibility[nozmax,5]

        shape = (self.PorRock.grid_numtot,self.PorRock.grid_numtot)

        self.Gmatrix = csr(shape)

        self.Gmatrix += csr((cx_m,(id_noxmin,idNnoxmin)),shape=shape)
        self.Gmatrix -= csr((cx_m,(id_noxmin,id_noxmin)),shape=shape)
        
        self.Gmatrix += csr((cx_p,(id_noxmax,idNnoxmax)),shape=shape)
        self.Gmatrix -= csr((cx_p,(id_noxmax,id_noxmax)),shape=shape)
        
        self.Gmatrix += csr((cy_m,(id_noymin,idNnoymin)),shape=shape)
        self.Gmatrix -= csr((cy_m,(id_noymin,id_noymin)),shape=shape)
        
        self.Gmatrix += csr((cy_p,(id_noymax,idNnoymax)),shape=shape)
        self.Gmatrix -= csr((cy_p,(id_noymax,id_noymax)),shape=shape)
        
        self.Gmatrix += csr((cz_m,(id_nozmin,idNnozmin)),shape=shape)
        self.Gmatrix -= csr((cz_m,(id_nozmin,id_nozmin)),shape=shape)

        self.Gmatrix += csr((cz_p,(id_nozmax,idNnozmax)),shape=shape)
        self.Gmatrix -= csr((cz_p,(id_nozmax,id_nozmax)),shape=shape)
        
    def set_externalBC(self,b_xmin=(0,1,0),b_xmax=(0,1,0),b_ymin=(0,1,0),b_ymax=(0,1,0),b_zmin=(0,1,0),b_zmax=(0,1,0)):

        """
        b_xmin,b_xmax,b_ymin,b_ymax,b_zmin and b_zmax have three entries:
        - dirichlet boundary condition coefficient,
        - neumann boundary condition coefficient,
        - function value of boundary condition
        Default is no flow boundary conditions at the exterior boundaries
        """

        sizes = self.PorRock.grid_sizes

        indices = self.PorRock.grid_indices

        self.b_correction = np.zeros(self.PorRock.grid_numtot)

        questbound = lambda x: True if x>1 else False

        if questbound(self.PorRock.grid_num[0]):

            xmin = indices[:,0]==indices[:,1]
            xmax = indices[:,0]==indices[:,2]

            id_xmin = indices[xmin,0]
            id_xmax = indices[xmax,0]

            dx_xmin = sizes[id_xmin,0]
            dx_xmax = sizes[id_xmax,0]

            tx_xmin = self.transmissibility[id_xmin,0]
            tx_xmax = self.transmissibility[id_xmax,1]

            bc_xmin = (2*tx_xmin*dx_xmin)/(b_xmin[0]*dx_xmin-2*b_xmin[1])
            bc_xmax = (2*tx_xmax*dx_xmax)/(b_xmax[0]*dx_xmax+2*b_xmax[1])
            
            self.Gmatrix[id_xmin,id_xmin] -= bc_xmin*b_xmin[0]
            self.Gmatrix[id_xmax,id_xmax] -= bc_xmax*b_xmax[0]
            
            self.b_correction[id_xmin] -= bc_xmin*b_xmin[2]
            self.b_correction[id_xmax] -= bc_xmax*b_xmax[2]

        if questbound(self.PorRock.grid_num[1]):

            ymin = indices[:,0]==indices[:,3]
            ymax = indices[:,0]==indices[:,4]

            id_ymin = indices[ymin,0]
            id_ymax = indices[ymax,0]

            dy_ymin = sizes[id_ymin,1]
            dy_ymax = sizes[id_ymax,1]

            ty_ymin = self.transmissibility[id_ymin,2]
            ty_ymax = self.transmissibility[id_ymax,3]

            bc_ymin = (2*ty_ymin*dy_ymin)/(b_ymin[0]*dy_ymin-2*b_ymin[1])
            bc_ymax = (2*ty_ymax*dy_ymax)/(b_ymax[0]*dy_ymax+2*b_ymax[1])
            
            self.Gmatrix[id_ymin,id_ymin] -= bc_ymin*b_ymin[0]
            self.Gmatrix[id_ymax,id_ymax] -= bc_ymax*b_ymax[0]
            
            self.b_correction[id_ymin] -= bc_ymin*b_ymin[2]
            self.b_correction[id_ymax] -= bc_ymax*b_ymax[2]

        try:
            grid_num_x = self.PorRock.grid_num[2]
        except IndexError:
            grid_num_x = 1
            
        if questbound(grid_num_x):

            zmin = indices[:,0]==indices[:,5]
            zmax = indices[:,0]==indices[:,6]

            id_zmin = indices[zmin,0]
            id_zmax = indices[zmax,0]

            dz_zmin = sizes[id_zmin,2]
            dz_zmax = sizes[id_zmax,2]

            tz_zmin = self.transmissibility[id_zmin,4]
            tz_zmax = self.transmissibility[id_zmax,5]

            bc_zmin = (2*tz_zmin*dz_zmin)/(b_zmin[0]*dz_zmin-2*b_zmin[1])
            bc_zmax = (2*tz_zmax*dz_zmax)/(b_zmax[0]*dz_zmax+2*b_zmax[1])

            self.Gmatrix[id_zmin,id_zmin] -= bc_zmin*b_zmin[0]
            self.Gmatrix[id_zmax,id_zmax] -= bc_zmax*b_zmax[0]
            
            self.b_correction[id_zmin] -= bc_zmin*b_zmin[2]
            self.b_correction[id_zmax] -= bc_zmax*b_zmax[2]

    def set_wells(self):

        self.well_grids      = np.array([],dtype=int)
        self.well_indices    = np.array([],dtype=int)
        self.well_bhpflags   = np.array([],dtype=bool)
        self.well_limits     = np.array([],dtype=float)

        wells = zip(
            self.Wells.Trajectory.tracks,
            self.Wells.consbhp,
            self.Wells.limits,
            )

        for index,(track,flag,limit) in enumerate(wells):

            ttrack = np.transpose(track[:,:,np.newaxis],(2,1,0))

            vector = self.PorRock.grid_centers[:,:,np.newaxis]-ttrack

            distance = np.sqrt(np.sum(vector**2,axis=1))

            well_grid = np.unique(np.argmin(distance,axis=0))

            well_index = np.full(well_grid.size,index,dtype=int)
    
            well_bhpflag = np.full(well_grid.size,flag,dtype=bool)

            well_limit = np.full(well_grid.size,limit,dtype=float)

            self.well_grids = np.append(self.well_grids,well_grid)

            self.well_indices = np.append(self.well_indices,well_index)

            self.well_bhpflags = np.append(self.well_bhpflags,well_bhpflag)

            self.well_limits = np.append(self.well_limits,well_limit)

        dx = self.PorRock.grid_sizes[self.well_grids,0]
        dy = self.PorRock.grid_sizes[self.well_grids,1]
        dz = self.PorRock.grid_sizes[self.well_grids,2]

        req = 0.14*np.sqrt(dx**2+dy**2)

        rw = self.Wells.radii[self.well_indices]

        skin = self.Wells.skinfactors[self.well_indices]

        kx = self.PorRock.permeability[:,0][self.well_grids]
        ky = self.PorRock.permeability[:,1][self.well_grids]

        dz = self.PorRock.grid_sizes[:,2][self.well_grids]

        self.JR = (2*np.pi*dz*np.sqrt(kx*ky))/(np.log(req/rw)+skin)
            
    def solve(self):

        mshape = (self.PorRock.grid_numtot,self.PorRock.grid_numtot)
        vshape = (self.PorRock.grid_numtot,1)

        vzeros = np.zeros(self.Wells.number,dtype=int)

        indices = self.PorRock.grid_indices

        shape = (self.PorRock.grid_numtot,self.PorRock.grid_numtot)

        oneperdtime = np.ones(self.PorRock.grid_numtot)/self.time_step
        
        t_correction = csr((oneperdtime,(indices[:,0],indices[:,0])),shape=shape)

        J_v = (self.JR)/(self.Fluids.viscosity[0]) #*self.Fluids.fvf[0]

        self.J = csr((J_v[self.well_bhpflags],(self.well_grids[self.well_bhpflags],self.well_grids[self.well_bhpflags])),shape=mshape)
        
        self.Q = csr(vshape)

        q_cp = self.well_limits[self.well_bhpflags]*J_v[self.well_bhpflags]

        q_cr = self.well_limits[~self.well_bhpflags]

        self.Q += csr((q_cp,(self.well_grids[self.well_bhpflags],vzeros[self.well_bhpflags])),shape=vshape)

        self.Q += csr((q_cr,(self.well_grids[~self.well_bhpflags],vzeros[~self.well_bhpflags])),shape=vshape)

        self.Q = self.Q.toarray().flatten()

        G = self.Gmatrix-t_correction

        CCC = self.Q/self.PorRock.grid_volumes/self.PorRock.porosity/self.compressibility
        
        for i in range(1,self.times.size):
            
            b = -self.pressure[:,i-1]/self.time_step+self.b_correction-CCC
            
            self.pressure[:,i] = sps(G,b)

    def postprocess(self):

        N = self.PorRock.grid_numtot

        Y = int((N-1)/2)

        Pwf = self.pressure[Y,:]+self.Q[Y]/self.JR*self.Fluids.viscosity[0]

        return Pwf

    @staticmethod
    def transmissibility(grids):
        
        dx_m = (grids.size[:,0]+grids.size[grids.index[:,1],0])/2
        dx_p = (grids.size[:,0]+grids.size[grids.index[:,2],0])/2
        dy_m = (grids.size[:,1]+grids.size[grids.index[:,3],1])/2
        dy_p = (grids.size[:,1]+grids.size[grids.index[:,4],1])/2
        dz_m = (grids.size[:,2]+grids.size[grids.index[:,5],2])/2
        dz_p = (grids.size[:,2]+grids.size[grids.index[:,6],2])/2

        kx_m = (2*dx_m)/(grids.size[:,0]/grids.permeability[:,0]+
                         grids.size[grids.index[:,1],0]/grids.permeability[grids.index[:,1],0])
        kx_p = (2*dx_p)/(grids.size[:,0]/grids.permeability[:,0]+
                         grids.size[grids.index[:,2],0]/grids.permeability[grids.index[:,2],0])
        ky_m = (2*dy_m)/(grids.size[:,1]/grids.permeability[:,1]+
                         grids.size[grids.index[:,3],1]/grids.permeability[grids.index[:,3],1])
        ky_p = (2*dy_p)/(grids.size[:,1]/grids.permeability[:,1]+
                         grids.size[grids.index[:,4],1]/grids.permeability[grids.index[:,4],1])
        kz_m = (2*dz_m)/(grids.size[:,2]/grids.permeability[:,2]+
                         grids.size[grids.index[:,5],2]/grids.permeability[grids.index[:,5],2])
        kz_p = (2*dz_p)/(grids.size[:,2]/grids.permeability[:,2]+
                         grids.size[grids.index[:,6],2]/grids.permeability[grids.index[:,6],2])

        transmissibility = np.zeros((grids.numtot,6))

        transmissibility[:,0] = (kx_m*grids.area[:,0])/(grids.viscosity*dx_m)
        transmissibility[:,1] = (kx_p*grids.area[:,1])/(grids.viscosity*dx_p)
        transmissibility[:,2] = (ky_m*grids.area[:,2])/(grids.viscosity*dy_m)
        transmissibility[:,3] = (ky_p*grids.area[:,3])/(grids.viscosity*dy_p)
        transmissibility[:,4] = (kz_m*grids.area[:,4])/(grids.viscosity*dz_m)
        transmissibility[:,5] = (kz_p*grids.area[:,5])/(grids.viscosity*dz_p)

        return transmissibility

    @staticmethod
    def transmatrix(grids,transmissibility):

        T = np.zeros((grids.numtot,grids.numtot))

        for i in range(grids.numtot):

            if i>0:
                T[i,i-1] = -transmissibility[i,0]
                T[i,i] += transmissibility[i,0]

            if i<grids.numtot-1:
                T[i,i+1] = -transmissibility[i,1]
                T[i,i] += transmissibility[i,1]

        return T

    @staticmethod
    def picardsolver(grids,timestep,timesteps,T,J,Q):

        array = np.zeros((grids.numtot,timesteps))

        P = grids.pressure_initial

        for j in range(timesteps):

            Pk = P.copy()

            firstIteration = True

            k = 1

            while firstIteration or error>1e-6:

                firstIteration = False

                A = np.eye(grids.numtot)*100/Pk

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

def NewtonSolver(grids,timestep,timesteps,T,J,Q):

    array = np.zeros((grids.numtot,timesteps))

    P = grids.pressure_initial

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

def PicardSolver(grids,timestep,timesteps,T,J,Q):

    pass

if __name__ == "__main__":

    import unittest

    from tests import test_porous_media

    unittest.main(test_porous_media)