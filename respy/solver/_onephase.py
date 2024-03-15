class OnePhase():
    """
    This class solves for single phase reservoir flow in Rectangular Cuboids;
        
        - Static Rock and Static Fluid Properties.
        - Static Rock and Dynamic Fluid Properties
        - Dynamic Rock and Dynamic Fluid Properties

    also it includes:

        - Implicit Pressure Solver
        - Explicit Pressure Solver
        - Mixed Pressure Solver
        
    """

    def __init__(self,grid,rrock,fluid,wconds=None,bconds=None,theta=0):
        """
        grid   : It is a RecCuboid (rectangular cuboid) object.

        rrock  : It is a class with reservoir rock properties.

        fluid  : There is only one mobile phase in the system. There can be
            two slightly compressible fluids where the second one is at irreducible
            saturation, not mobile.

        wconds : tuple of WellCond item.

        bconds : tuple of BoundCond item.

        theta  : solution type determined in the mixed solver
            when theta = 0, it reduces to Implicit method
            when theta = 1, it reduces to Explicit method
            when theta = 1/2, it is Crank-Nicolson method

        """

        self.grid = grid

        self.rrock = rrock
        self.fluid = fluid

        self.wconds = () if wconds is None else wconds
        self.bconds = () if bconds is None else bconds

        self.theta = theta

    def init(self,dref,pref,depth,gradf):
        """Calculates the initial pressure
        
        dref    : reference depth, ft
        pref    : reference pressure, psi

        depth   : depths where to calculate the pressure, ft
        gradf   : fluid hydrostatic gradient, psi/ft
        """
        self._pinit = (pref+gradf*(depth-dref))*6894.75729

    def set_time(self,tstep:float,total:float=None,nstep:int=1):
        """Setting the numerical parameters

        tstep   : time step defined in days
        total   : total simulation time defined in days
        nstep   : number of total time steps
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

    def set_tcomp(self,tcomp):
        """
        tcomp  : total compressibility 1/psi
        """

        self._tcomp = tcomp/6894.75729

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
    def pinit(self):
        return self._pinit/6894.75729

    @property
    def tcomp(self):
        return self._tcomp*6894.75729

    def __call__(self,press):

        fluid, = self.fluids

        phase = fluid(press,None)

        for wcond in self.wconds:
            wcond.cube = self.cube[wcond.block]

        for bcond in self.bconds:
            bcond.cube = getattr(getattr(self.cube,bcond.face),"rows")

        T = self.get_Tmatrix(phase)
        S = self.get_Smatrix()
        G = self.get_Gvector(phase)
        J = self.get_Jmatrix(phase)
        Q = self.get_Qvector(phase)

        return Matrix(T,V,G,J,Q,self.tcomp,self.tstep)

    def solve(self,**kwargs):

        pass

    def static(self):

        P = self._pinit

        vec = Block(cube,P,tstep,comp)
        mat = Shape(cube,P,vec)

        self._pressure = numpy.zeros((self.grid.numtot,self._time.size))
        
        for k in range(self.nstep):

            P = mat.implicit_pressure(P)

            print(f"{k} time step is complete...")
            
            self._pressure[:,k] = P

            P = P.reshape((-1,1))

    def picard(self):

        pass

    def newton(self):

        pass

    @property
    def pressure(self):
        return self._pressure/6894.75729

    def postprocess(self):

        Y = int((self.grid.numtot-1)/2)

        Pwf = self.pressure[Y,:]+self.Q[Y]/self.JR*self.Fluids.viscosity[0]

        return Pwf
