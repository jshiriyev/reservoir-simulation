import matplotlib.pyplot as plt

import numpy

from scipy.sparse import linalg

from scipy.sparse import identity
from scipy.sparse import diags
from scipy.sparse import csr_matrix as csr

# from ._relperm import RelPerm

# from ._cappres import BrooksCorey
# from ._cappres import VanGenuchten
# from ._cappres import JFunction
# from ._cappres import ScanCurves

"""
The solver must include:

- One-Phase Slightly Compressible Flow
- One-Phase Compressible Flow
- Multi-Phase Flow

"""

class WellCond:
    """
    It is a well condition object used in the simulator.
    """

    def __init__(self,block:tuple,axis:str,radius:float,skin:float,start:float,end:float,**kwargs):
        """
        block   : block indices containing the well
        axis    : (z) vertical or (x,y) horizontal well

        radius  : well radius, ft
        skin    : skin factor of the well, dimensionless

        start   : start time for implementing the condition, days
        end     : end time for implementing the condition, days

        Assign only one of the following conditions:

        bhp     : constant bottom hole pressure
        whp     : constant well head pressure

        orate   : constant oil rate
        wrate   : constant water rate
        grate   : constant gas rate
        """

        self._block     = block
        self._axis      = axis
        self._radius    = radius*0.3048
        self._skin      = skin
        self._start     = start*(24*60*60)
        self._end       = end*(24*60*60)

        for key,value in kwargs.items():
            if value is not None:
                self._sort = key
                if key == "press":
                    self._cond = value*6894.76
                elif key in ("orate","grate"):
                    self._cond = value*1.84013e-6
                elif key == "grate":
                    self._cond = value*3.27741e-7
                break

    @property
    def block(self):
        return self._block

    @property
    def axis(self):
        return self._axis

    @property
    def radius(self):
        return self._radius/0.3048

    @property
    def skin(self):
        return self._skin

    @property
    def start(self):
        return self._start/(24*60*60)

    @property
    def end(self):
        return self._end/(24*60*60)

    @property
    def sort(self):
        return self._sort
    
    @property
    def cond(self):
        if self._sort in ("bhp","whp"):
            return self._cond/6894.76
        elif self._sort in ("orate","wrate"):
            return self._cond/1.84013e-6
        elif self._sort == "grate":
            return self._cond/3.27741e-7

class BoundCond():

    def __init__(self,face,**kwargs):
        """
        face    : boundary, xmin, xmax, ymin, ymax, zmin, or zmax

        Assign only one of the following conditions:

        press   : constant pressure values
        orate   : constant flow boundary condition, 0 = no flow
        wrate   : constant flow boundary condition, 0 = no flow
        grate   : constant flow boundary condition, 0 = no flow
        """
        self._face = face

        for key,value in kwargs.items():
            if value is not None:
                self._sort = key
                if key == "press":
                    self._cond = value*6894.76
                elif key in ("orate","wrate"):
                    self._cond = value*1.84013e-6
                elif key == "grate":
                    self._cond = value*3.27741e-7
                break

    @property
    def face(self):
        return self._face

    @property
    def sort(self):
        return self._sort
    
    @property
    def cond(self):
        if self._sort == "press":
            return self._cond/6894.76
        elif self._sort in ("orate","wrate"):
            return self._cond/1.84013e-6
        elif self._sort == "grate":
            return self._cond/3.27741e-7

class Matrix():

    def __init__(self,grid,fluid,wconds=None,bconds=None):
        """
        grid    : RecCuboid instance, rectangular cuboid grids

        fluid   : fluid item defining properties and methods for calculations
        
        wcond   : well condition, tuple of WellCond instance
        bcond   : boundary condition, tuple of BoundCond instance
        """

        self.grid = grid

        self.fluid = fluid

        self.wconds = () if wconds is None else wconds
        self.bconds = () if bconds is None else bconds

        self.set_static()

    def set_static(self):
        """Self assigns static transmissibility values."""

        self._staticx = self.set_stat_axis('x')

        if self.grid.flodim>1:
            self._staticy = self.set_stat_axis('y')

        if self.grid.flodim>2:
            self._staticz = self.set_stat_axis('z')

        self._staticw = [self.set_stat_well(wcond) for wcond in self.wconds]
        self._staticb = [self.set_stat_face(bcond) for bcond in self.bconds]

        for bcond in self.bconds:
            bcond.block = getattr(self.grid,bcond.face)

    def Tmatrix(self):
        """Returns transmissibility matrix with dynamic transmissibility values."""

        tmatrix = csr(self.matrix)

        tmatrix = self.tcharge('x',tmatrix)

        if self.grid.flodim>1:
            tmatrix = self.tcharge('y',tmatrix)

        if self.grid.flodim>2:
            tmatrix = self.tcharge('z',tmatrix)

        return tmatrix

    def Jmatrix(self):
        """
        Returns J matrix filled with constant pressure boundary indices.
        """

        jmatrix = csr(self.matrix)

        for index,wcond in enumerate(self.wconds):
            jmatrix += self.jcharge(wcond,self._staticw[index],jmatrix)

        for index,bcond in enumerate(self.bconds):
            jmatrix += self.jcharge(bcond,self._staticb[index],jmatrix)

        return jmatrix

    def Amatrix(self,tstep):

        pore_volume = self.grid._volume*self.grid._poro

        return diags(pore_volume.flatten()/tstep)

    def Qvector(self):
        """
        Returns Q vector filled with constant pressure values.
        """

        qvector = csr(self.vector)

        for index,wcond in enumerate(self.wconds):
            qvector += self.qcharge(wcond,self._staticw[index],qvector)

        for index,bcond in enumerate(self.bconds):
            qvector += self.qcharge(bcond,self._staticb[index],qvector)

        return qvector

    def Gvector(self,tmatrix):

        return self.fluid._density*self._gravity*tmatrix.dot(self.grid._depth)

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

    def tcharge(self,axis:str,tmatrix:csr):
        """
        Returns updated transmissibility matrix:

        axis    : x, y, or z

        tmatrix : csr_matrix object defined for transmissibility matrix

        The transmissibility matrix is updated for diagonal and offset entries.
        """

        stat = getattr(self,f"_static{axis}")

        vals = stat/self.fluid._viscosity

        dneg = getattr(self.grid,f"{axis}neg")
        dpos = getattr(self.grid,f"{axis}pos")

        tmatrix += csr((vals,(dneg,dneg)),shape=tmatrix.shape)
        tmatrix -= csr((vals,(dneg,dpos)),shape=tmatrix.shape)

        tmatrix += csr((vals,(dpos,dpos)),shape=tmatrix.shape)
        tmatrix -= csr((vals,(dpos,dneg)),shape=tmatrix.shape)

        return tmatrix

    def jcharge(self,cond,static,jmatrix:csr):

        if cond.sort=="press":

            jvalues = 2*static/self.fluid._viscosity
            jmatrix += csr((jvalues,(cond.block,cond.block)),shape=self.matrix)

        return jmatrix

    def qcharge(self,cond,static,qvector:csr):

        if cond.sort=="press":
            qvalues = 2*static/self.fluid._viscosity*cond._cond
        else:
            qvalues = cond._cond

        indices = numpy.zeros(cond.block.size)

        qvector += csr((qvalues,(cond.block,indices)),shape=self.vector)

        return qvector

    def set_stat_axis(self,axis):
        """Returns static transmissibility values for the given
        direction and inner faces (interfaces)."""
        
        dims = self.stat_dims(
            getattr(getattr(self.grid,f"{axis}neg"),f"_{axis}dims"),
            getattr(getattr(self.grid,f"{axis}pos"),f"_{axis}dims"))

        area = self.stat_area(
            getattr(getattr(self.grid,f"{axis}neg"),f"_{axis}area"),
            getattr(getattr(self.grid,f"{axis}pos"),f"_{axis}area"))
        
        perm = self.stat_perm(
            getattr(getattr(self.grid,f"{axis}neg"),f"_{axis}dims"),
            getattr(getattr(self.grid,f"{axis}pos"),f"_{axis}dims"),
            getattr(getattr(self.grid,f"{axis}neg"),f"_{axis}perm"),
            getattr(getattr(self.grid,f"{axis}pos"),f"_{axis}perm"))

        return self.stat_stat(dims,area,perm)

    def set_stat_well(self,wcond):
        """Returns static transmissibility values for the given
        vertical well."""

        dx = self.grid._xdims[wcond.block]
        dy = self.grid._ydims[wcond.block]
        dz = self.grid._zdims[wcond.block]

        kx = self.grid._xperm[wcond.block]
        ky = self.grid._yperm[wcond.block]
        kz = self.grid._zperm[wcond.block]

        if wcond.axis == "x":
            dkh = dx*numpy.sqrt(ky*kz)
            req = self.stat_well(dy,dz,ky,kz)
        elif wcond.axis == "y":
            dkh = dy*numpy.sqrt(kz*kx)
            req = self.stat_well(dz,dx,kz,kx)
        elif wcond.axis == "z":
            dkh = dz*numpy.sqrt(kx*ky)
            req = self.stat_well(dx,dy,kx,ky)

        return (2*numpy.pi*dkh)/(numpy.log(req/wcond.radius)+wcond.skin)

    def set_stat_face(self,bcond):
        """Returns static transmissibility values for the given
        surface on the exterior boundary."""
        
        dims = getattr(getattr(self.grid,bcond),f"_{bcond[0]}dims")
        area = getattr(getattr(self.grid,bcond),f"_{bcond[0]}area")
        perm = getattr(getattr(self.grid,bcond),f"_{bcond[0]}perm")

        return self.stat_stat(dims,area,perm)

    @staticmethod
    def stat_dims(dims_neg,dims_pos):
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

    @staticmethod
    def stat_well(dims1,dims2,perm1,perm2):

        sqrt21 = numpy.power(perm2/perm1,1/2)*numpy.power(dims1,2)
        sqrt12 = numpy.power(perm1/perm2,1/2)*numpy.power(dims2,2)

        quar21 = numpy.power(perm2/perm1,1/4)
        quar12 = numpy.power(perm1/perm2,1/4)

        return 0.28*numpy.sqrt(sqrt21+sqrt12)/(quar21+quar12)

class OnePhase(Matrix):
    """
    This class solves for single phase reservoir flow in Rectangular Cuboids.
    """

    def __init__(self,grid,fluid,wconds=None,bconds=None,theta=0):
        """
        grid  : It is a RecCuboid (rectangular cuboid) object.

        fluid : There is only one mobile phase in the system. There can be
            two slightly compressible fluids where the second one is at irreducible
            saturation, not mobile.

        wcond : tuple of WellCond item.

        bcond : tuple of BoundCond item.

        theta : solution type determined in the mixed solver
            when theta = 0, it reduces to Implicit method
            when theta = 1, it reduces to Explicit method
            when theta = 1/2, it is Crank-Nicolson method

        """

        super().__init__(grid,fluid,wconds,bconds)

        self._theta = theta

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
        J = self.Jmatrix(self.fluid)
        A = self.Amatrix(self._tstep)
        Q = self.Qvector(self.fluid)
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
    def theta(self):
        return self._theta

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

class BlackOil():
    """IMPES Solution"""
    def __init__(self,res,fluids,relperm,wells):

        self.PorRock = PorRock("rectangle")()

        # There can be two slightly compressible fluids where the
        # second one is at irreducible saturation, not mobile

        self.Fluids = Fluids(number=2)

        self.Wells = Wells()

        self.rp = relperm

    def set_transmissibility(self):

        dx0 = self.res.grid_sizes[:,0]
        dy0 = self.res.grid_sizes[:,1]
        dz0 = self.res.grid_sizes[:,2]

        dxm = self.res.grid_sizes[self.res.grid_indices[:,1],0]
        dxp = self.res.grid_sizes[self.res.grid_indices[:,2],0]
        dym = self.res.grid_sizes[self.res.grid_indices[:,3],1]
        dyp = self.res.grid_sizes[self.res.grid_indices[:,4],1]
        dzm = self.res.grid_sizes[self.res.grid_indices[:,5],2]
        dzp = self.res.grid_sizes[self.res.grid_indices[:,6],2]

        Ax = self.res.grid_areas[:,0]
        Ay = self.res.grid_areas[:,1]
        Az = self.res.grid_areas[:,2]

        kx0 = self.res.permeability[:,0]
        ky0 = self.res.permeability[:,1]
        kz0 = self.res.permeability[:,2]

        kxm = self.res.permeability[self.res.grid_indices[:,1],0]
        kxp = self.res.permeability[self.res.grid_indices[:,2],0]
        kym = self.res.permeability[self.res.grid_indices[:,3],1]
        kyp = self.res.permeability[self.res.grid_indices[:,4],1]
        kzm = self.res.permeability[self.res.grid_indices[:,5],2]
        kzp = self.res.permeability[self.res.grid_indices[:,6],2]

        dxm_mean = (dx0+dxm)/2
        dxp_mean = (dx0+dxp)/2
        dym_mean = (dy0+dym)/2
        dyp_mean = (dy0+dyp)/2
        dzm_mean = (dz0+dzm)/2
        dzp_mean = (dz0+dzp)/2

        kxm_mean = (2*dxm_mean)/(dx0/kx0+dxm/kxm)
        kxp_mean = (2*dxp_mean)/(dx0/kx0+dxp/kxp)
        kym_mean = (2*dym_mean)/(dy0/ky0+dym/kym)
        kyp_mean = (2*dyp_mean)/(dy0/ky0+dyp/kyp)
        kzm_mean = (2*dzm_mean)/(dz0/kz0+dzm/kzm)
        kzp_mean = (2*dzp_mean)/(dz0/kz0+dzp/kzp)

        self.TXM = (Ax*kxm_mean)/(dxm_mean)
        self.TYM = (Ay*kym_mean)/(dym_mean)
        self.TZM = (Az*kzm_mean)/(dzm_mean)

        self.TXP = (Ax*kxp_mean)/(dxp_mean)
        self.TYP = (Ay*kyp_mean)/(dyp_mean)
        self.TZP = (Az*kzp_mean)/(dzp_mean)

    def set_wells(self):

        self.well_grids      = np.array([],dtype=int)
        self.well_indices    = np.array([],dtype=int)
        self.well_bhpflags   = np.array([],dtype=bool)
        self.well_limits     = np.array([],dtype=float)
        self.well_waterflags = np.array([],dtype=bool)
        self.well_oilflags   = np.array([],dtype=bool)

        wells = zip(self.wells.tracks,self.wells.consbhp,self.wells.limits,self.wells.water,self.wells.oil)

        for index,(track,flag,limit,water,oil) in enumerate(wells):

            ttrack = np.transpose(track[:,:,np.newaxis],(2,1,0))

            vector = self.res.grid_centers[:,:,np.newaxis]-ttrack

            distance = np.sqrt(np.sum(vector**2,axis=1))

            well_grid = np.unique(np.argmin(distance,axis=0))

            well_index = np.full(well_grid.size,index,dtype=int)
    
            well_bhpflag = np.full(well_grid.size,flag,dtype=bool)

            well_limit = np.full(well_grid.size,limit,dtype=float)

            well_waterflag = np.full(well_grid.size,water,dtype=bool)
            
            well_oilflag = np.full(well_grid.size,oil,dtype=bool)

            self.well_grids = np.append(self.well_grids,well_grid)

            self.well_indices = np.append(self.well_indices,well_index)

            self.well_bhpflags = np.append(self.well_bhpflags,well_bhpflag)

            self.well_limits = np.append(self.well_limits,well_limit)

            self.well_waterflags = np.append(self.well_waterflags,well_waterflag)

            self.well_oilflags = np.append(self.well_oilflags,well_oilflag)

        dx = self.res.grid_sizes[self.well_grids,0]
        dy = self.res.grid_sizes[self.well_grids,1]
        dz = self.res.grid_sizes[self.well_grids,2]

        req = 0.14*np.sqrt(dx**2+dy**2)

        rw = self.wells.radii[self.well_indices]

        skin = self.wells.skinfactors[self.well_indices]

        kx = self.res.permeability[:,0][self.well_grids]
        ky = self.res.permeability[:,1][self.well_grids]

        dz = self.res.grid_sizes[:,2][self.well_grids]

        self.JR = (2*np.pi*dz*np.sqrt(kx*ky))/(np.log(req/rw)+skin)

    def set_time(self,step,total):

        self.time_step = step
        self.time_total = total

        self.time_array = np.arange(self.time_step,self.time_total+self.time_step,self.time_step)

    def initialize(self,pressure,saturation):

        N = self.time_array.size

        self.pressure = np.empty((self.res.grid_numtot,N+1))
        self.pressure[:,0] = pressure

        self.Sw = np.empty((self.res.grid_numtot,N+1))
        self.Sw[:,0] = saturation

        # for index,(p,sw) in enumerate(zip(self.pressure[:,-1],self.Sw[:,-1])):
        #   print("{:d}\tP\t{:3.0f}\tSw\t{:.4f}".format(index,p,sw))

    def solve(self):

        Vp = self.res.grid_volumes*self.res.porosity

        muw = self.fluids.viscosity[0]
        muo = self.fluids.viscosity[1]
        
        fvfw = self.fluids.fvf[0]
        fvfo = self.fluids.fvf[1]
        
        cr = self.res.compressibility

        cw = self.fluids.compressibility[0]
        co = self.fluids.compressibility[1]

        wflow_1phase = ~np.logical_and(self.well_waterflags,self.well_oilflags)

        sp_wf = np.logical_and(self.well_waterflags,wflow_1phase)
        sp_of = np.logical_and(self.well_oilflags,wflow_1phase)

        cpw = np.logical_and(self.well_bhpflags,self.well_waterflags)
        cpo = np.logical_and(self.well_bhpflags,self.well_oilflags)

        cqw = np.logical_and(~self.well_bhpflags,self.well_waterflags)
        cqo = np.logical_and(~self.well_bhpflags,self.well_oilflags)

        cxm_0 = self.res.grid_indices[self.res.grid_hasxmin,0]
        cxm_m = self.res.grid_indices[self.res.grid_hasxmin,1]

        cxp_0 = self.res.grid_indices[self.res.grid_hasxmax,0]
        cxp_p = self.res.grid_indices[self.res.grid_hasxmax,2]

        cym_0 = self.res.grid_indices[self.res.grid_hasymin,0]
        cym_m = self.res.grid_indices[self.res.grid_hasymin,3]

        cyp_0 = self.res.grid_indices[self.res.grid_hasymax,0]
        cyp_p = self.res.grid_indices[self.res.grid_hasymax,4]

        czm_0 = self.res.grid_indices[self.res.grid_haszmin,0]
        czm_m = self.res.grid_indices[self.res.grid_haszmin,5]

        czp_0 = self.res.grid_indices[self.res.grid_haszmax,0]
        czp_p = self.res.grid_indices[self.res.grid_haszmax,6]

        mshape = (self.res.grid_numtot,self.res.grid_numtot)
        vshape = (self.res.grid_numtot,1)

        vzeros = np.zeros(len(self.wells.itemnames),dtype=int)

        krw,kro = self.rp.water_oil(Sw=self.Sw[:,0])

        krw_xm = krw[cxm_0]
        krw_xp = krw[cxp_0]
        krw_ym = krw[cym_0]
        krw_yp = krw[cyp_0]
        krw_zm = krw[czm_0]
        krw_zp = krw[czp_0]

        kro_xm = kro[cxm_0]
        kro_xp = kro[cxp_0]
        kro_ym = kro[cym_0]
        kro_yp = kro[cyp_0]
        kro_zm = kro[czm_0]
        kro_zp = kro[czp_0]

        for index,time in enumerate(self.time_array):

            print("@{:5.1f}th time step".format(time))

            d11 = (Vp*self.Sw[:,index])/(self.time_step*fvfw)*(cr+cw)
            d12 = (Vp)/(self.time_step*fvfw)
            d21 = (Vp*(1-self.Sw[:,index]))/(self.time_step*fvfo)*(cr+co)
            d22 = (Vp)/(self.time_step*fvfo)*-1

            self.D = diags(-d22/d12*d11+d21)

            Uxm = self.pressure[:,index][cxm_0]>self.pressure[:,index][cxm_m]
            Uxp = self.pressure[:,index][cxp_0]>self.pressure[:,index][cxp_p]
            Uym = self.pressure[:,index][cym_0]>self.pressure[:,index][cym_m]
            Uyp = self.pressure[:,index][cyp_0]>self.pressure[:,index][cyp_p]
            Uzm = self.pressure[:,index][czm_0]>self.pressure[:,index][czm_m]
            Uzp = self.pressure[:,index][czp_0]>self.pressure[:,index][czp_p]

            krw,kro = self.rp.water_oil(Sw=self.Sw[:,index])

            krw_xm[Uxm] = krw[cxm_0][Uxm]
            krw_xp[Uxp] = krw[cxp_0][Uxp]
            krw_ym[Uym] = krw[cym_0][Uym]
            krw_yp[Uyp] = krw[cyp_0][Uyp]
            krw_zm[Uzm] = krw[czm_0][Uzm]
            krw_zp[Uzp] = krw[czp_0][Uzp]

            krw_xm[~Uxm] = krw[cxm_m][~Uxm]
            krw_xp[~Uxp] = krw[cxp_p][~Uxp]
            krw_ym[~Uym] = krw[cym_m][~Uym]
            krw_yp[~Uyp] = krw[cyp_p][~Uyp]
            krw_zm[~Uzm] = krw[czm_m][~Uzm]
            krw_zp[~Uzp] = krw[czp_p][~Uzp]

            kro_xm[Uxm] = kro[cxm_0][Uxm]
            kro_xp[Uxp] = kro[cxp_0][Uxp]
            kro_ym[Uym] = kro[cym_0][Uym]
            kro_yp[Uyp] = kro[cyp_0][Uyp]
            kro_zm[Uzm] = kro[czm_0][Uzm]
            kro_zp[Uzp] = kro[czp_0][Uzp]

            kro_xm[~Uxm] = kro[cxm_m][~Uxm]
            kro_xp[~Uxp] = kro[cxp_p][~Uxp]
            kro_ym[~Uym] = kro[cym_m][~Uym]
            kro_yp[~Uyp] = kro[cyp_p][~Uyp]
            kro_zm[~Uzm] = kro[czm_m][~Uzm]
            kro_zp[~Uzp] = kro[czp_p][~Uzp]

            TXMw = (self.TXM[self.res.grid_hasxmin]*krw_xm)/(muw*fvfw)*6.33e-3 # unit conversion
            TYMw = (self.TYM[self.res.grid_hasymin]*krw_ym)/(muw*fvfw)*6.33e-3 # unit conversion
            TZMw = (self.TZM[self.res.grid_haszmin]*krw_zm)/(muw*fvfw)*6.33e-3 # unit conversion

            TXPw = (self.TXP[self.res.grid_hasxmax]*krw_xp)/(muw*fvfw)*6.33e-3 # unit conversion
            TYPw = (self.TYP[self.res.grid_hasymax]*krw_yp)/(muw*fvfw)*6.33e-3 # unit conversion
            TZPw = (self.TZP[self.res.grid_haszmax]*krw_zp)/(muw*fvfw)*6.33e-3 # unit conversion

            TXMn = (self.TXM[self.res.grid_hasxmin]*kro_xm)/(muo*fvfo)*6.33e-3 # unit conversion
            TYMn = (self.TYM[self.res.grid_hasymin]*kro_ym)/(muo*fvfo)*6.33e-3 # unit conversion
            TZMn = (self.TZM[self.res.grid_haszmin]*kro_zm)/(muo*fvfo)*6.33e-3 # unit conversion

            TXPn = (self.TXP[self.res.grid_hasxmax]*kro_xp)/(muo*fvfo)*6.33e-3 # unit conversion
            TYPn = (self.TYP[self.res.grid_hasymax]*kro_yp)/(muo*fvfo)*6.33e-3 # unit conversion
            TZPn = (self.TZP[self.res.grid_haszmax]*kro_zp)/(muo*fvfo)*6.33e-3 # unit conversion

            self.Tw = csr(mshape)

            self.Tw -= csr((TXMw,(cxm_0,cxm_m)),shape=mshape)
            self.Tw += csr((TXMw,(cxm_0,cxm_0)),shape=mshape)

            self.Tw -= csr((TXPw,(cxp_0,cxp_p)),shape=mshape)
            self.Tw += csr((TXPw,(cxp_0,cxp_0)),shape=mshape)

            self.Tw -= csr((TYMw,(cym_0,cym_m)),shape=mshape)
            self.Tw += csr((TYMw,(cym_0,cym_0)),shape=mshape)

            self.Tw -= csr((TYPw,(cyp_0,cyp_p)),shape=mshape)
            self.Tw += csr((TYPw,(cyp_0,cyp_0)),shape=mshape)

            self.Tw -= csr((TZMw,(czm_0,czm_m)),shape=mshape)
            self.Tw += csr((TZMw,(czm_0,czm_0)),shape=mshape)

            self.Tw -= csr((TZPw,(czp_0,czp_p)),shape=mshape)
            self.Tw += csr((TZPw,(czp_0,czp_0)),shape=mshape)

            self.Tn = csr(mshape)

            self.Tn -= csr((TXMn,(cxm_0,cxm_m)),shape=mshape)
            self.Tn += csr((TXMn,(cxm_0,cxm_0)),shape=mshape)

            self.Tn -= csr((TXPn,(cxp_0,cxp_p)),shape=mshape)
            self.Tn += csr((TXPn,(cxp_0,cxp_0)),shape=mshape)

            self.Tn -= csr((TYMn,(cym_0,cym_m)),shape=mshape)
            self.Tn += csr((TYMn,(cym_0,cym_0)),shape=mshape)

            self.Tn -= csr((TYPn,(cyp_0,cyp_p)),shape=mshape)
            self.Tn += csr((TYPn,(cyp_0,cyp_0)),shape=mshape)

            self.Tn -= csr((TZMn,(czm_0,czm_m)),shape=mshape)
            self.Tn += csr((TZMn,(czm_0,czm_0)),shape=mshape)

            self.Tn -= csr((TZPn,(czp_0,czp_p)),shape=mshape)
            self.Tn += csr((TZPn,(czp_0,czp_0)),shape=mshape)

            self.T = diags(-d22/d12)*self.Tw+self.Tn

            krw,kro = self.rp.water_oil(Sw=self.Sw[:,index][self.well_grids])

            Jw_v = (self.JR*krw)/(muw*fvfw)*6.33e-3 # unit conversion
            Jn_v = (self.JR*kro)/(muo*fvfo)*6.33e-3 # unit conversion

            self.Jw = csr((Jw_v[cpw],(self.well_grids[cpw],self.well_grids[cpw])),shape=mshape)
            self.Jn = csr((Jn_v[cpo],(self.well_grids[cpo],self.well_grids[cpo])),shape=mshape)

            self.J = diags(-d22/d12)*self.Jw+self.Jn
            
            self.Qw = csr(vshape)
            self.Qn = csr(vshape)

            qw_cp = self.well_limits[cpw]*Jw_v[cpw]
            qo_cp = self.well_limits[cpo]*Jn_v[cpo]

            qw_cr = self.well_limits[cqw]*(krw[cqw]*muo)/(krw[cqw]*muo+kro[cqw]*muw)
            qo_cr = self.well_limits[cqo]*(kro[cqo]*muw)/(krw[cqo]*muo+kro[cqo]*muw)

            qw_cr[sp_wf[cqw]] = self.well_limits[sp_wf]
            qo_cr[sp_of[cqo]] = self.well_limits[sp_of]

            self.Qw += csr((qw_cp,(self.well_grids[cpw],vzeros[cpw])),shape=vshape)
            self.Qn += csr((qo_cp,(self.well_grids[cpo],vzeros[cpo])),shape=vshape)

            self.Qw += csr((qw_cr*5.61,(self.well_grids[cqw],vzeros[cqw])),shape=vshape) # unit conversion
            self.Qn += csr((qo_cr*5.61,(self.well_grids[cqo],vzeros[cqo])),shape=vshape) # unit conversion

            self.Qw = self.Qw.toarray().flatten()
            self.Qn = self.Qn.toarray().flatten()

            self.Q = -d22/d12*self.Qw+self.Qn

            self.pressure[:,index+1] = sps(self.T+self.J+self.D,self.D.dot(self.pressure[:,index])+self.Q)

            delta_p = (self.pressure[:,index+1]-self.pressure[:,index])
            
            tjp = csr.dot(self.Tw+self.Jw,self.pressure[:,index+1])

            self.Sw[:,index+1] = self.Sw[:,index]-(d11*delta_p-self.Qw+tjp)/d12

        for index,(p,sw) in enumerate(zip(self.pressure[:,-1],self.Sw[:,-1])):
            print("{:d}\tP\t{:4.1f}\tSw\t{:.5f}".format(index,p,sw))

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