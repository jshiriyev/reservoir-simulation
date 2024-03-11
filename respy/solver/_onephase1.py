class OnePhase1():
    """
    This class solves for single phase reservoir flow in Rectangular Cuboids;
    Static Rock and Fluid Properties.
    """

    _gravity = 9.807  # Gravitational acceleration in SI units

    def __init__(self,grid,rock,fluid,wconds=None,bconds=None,theta=0):
        """
        grid  : It is a RecCuboid (rectangular cuboid) object.

        rock  : It is a 

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

        self.grid = grid

        self.rock = rock

        self.set_static()

        self.fluid = fluid

        self.wconds = () if wconds is None else wconds
        self.bconds = () if bconds is None else bconds

        self._theta = theta

    @property
    def theta(self):
        return self._theta

    def get_property(self,quality,coeff=1.,dtype=None):
        """coeff    : conversion factor"""

        quality = numpy.asarray(quality)

        if dtype is not None:
            quality = quality.astype(dtype)

        quality = quality.flatten()*coeff

        if quality.size==1:
            quality = quality.repeat(self.numtot)

        return quality.reshape((-1,1))

    def set_depth(self,depth):
        """Assigns the depth values in ft to the grids."""

        self.grid.cube.cvol.set_prop(
            depth=self.get_property(depth,coeff=0.3048,dtype=numpy.float_))

    def set_poro(self,poro):
        """Assigns the porosity values in fractions to the grids."""

        self.grid.cube.cvol.set_prop(
            poro=self.get_property(poro,dtype=numpy.float_))

    def set_perm(self,xperm,yperm=None,zperm=None,yreduce=1.,zreduce=1.):
        """Assigns the permeability values in mD to the grids."""

        self.grid.cube.cvol.set_prop(
            xperm=self.get_property(xperm,coeff=9.869233e-16,dtype=numpy.float_))

        self.grid.cube.cvol.set_prop(
            yperm=self.xperm*yreduce if yperm is None else self.get_property(
            yperm,coeff=9.869233e-16,dtype=numpy.float_))

        self.grid.cube.cvol.set_prop(
            zperm=self.xperm*zreduce if zperm is None else self.get_property(
                zperm,coeff=9.869233e-16,dtype=numpy.float_))

    def set_comp(self,P):
        return self.__comp(P)

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

    def init(self,dref,pref,tcomp):
        """Initializing the reservoir pressure
        
        dref    : reference depth in feet
        pref    : reference pressure in psi

        tcomp  : total compressibility 1/psi
        """

        self._pinit = self.pzero(
            dref*0.3048,pref*6894.75729,self.grid._depth,self.fluid._density)

        self._tcomp = tcomp/6894.75729

    @property
    def pinit(self):
        return self._pinit/6894.75729

    @property
    def tcomp(self):
        return self._tcomp*6894.75729

    def __call__(self,P):

        fluid, = self.fluids

        phase = fluid(P,None)

        T = self.get_Tmatrix(phase)
        S = self.get_Smatrix()
        G = self.get_Gvector(phase)
        J = self.get_Jmatrix(phase)
        Q = self.get_Qvector(phase)

        return Matrix(T,V,G,J,Q,self.tcomp,self.tstep)

    def solve(self):

        P = self._pinit

        self._pressure = numpy.zeros((self.grid.numtot,self._time.size))
        
        for k in range(self.nstep):

            mat = self(P,self.tstep)

            if self.theta==0:
                P = self.implicit(P,T,J,Act,Q,G)
            elif self.theta==1:
                P = self.explicit(P,T,J,Act,Q,G)
            else:
                P = self.mixed(P,T,J,Act,Q,G,self.theta)

            print(f"{k} time step is complete...")
            
            self._pressure[:,k] = P

            P = P.reshape((-1,1))

    @property
    def pressure(self):
        return self._pressure/6894.75729

    def postprocess(self):

        Y = int((self.grid.numtot-1)/2)

        Pwf = self.pressure[Y,:]+self.Q[Y]/self.JR*self.Fluids.viscosity[0]

        return Pwf

    @staticmethod
    def pzero(dref,pref,depths,density):
        """Calculates the initial pressure
        
        dref    : reference depth in m
        pref    : reference pressure in Pa

        depths  : depths where to calculate the pressure in m
        density : fluid density in kg/m3
        """
        return pref+density*OnePhase._gravity*(depths-dref)
