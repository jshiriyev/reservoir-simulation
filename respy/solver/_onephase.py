class OnePhase(ResRock):
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

    def initialize(self,dref,pref,ctotal):
        """Initializing the reservoir pressure
        
        dref    : reference depth in feet
        pref    : reference pressure in psi

        ctotal  : total compressibility 1/psi
        """

        self._pinit = self.pzero(
            dref*0.3048,pref*6894.76,self.grid._depth,self.fluid._density)

        self._ctotal = ctotal/6894.76

    def solve(self):

        P = self._pinit

        self.update()

        T = self._Tmatrix
        J = self._Jmatrix
        A = self._Amatrix
        Q = self._Qvector
        G = self._Gvector

        Act = A*(self._ctotal/self._tstep)

        self._pressure = numpy.zeros((self.grid.numtot,self._time.size))
        
        for k in range(self.nstep):

            if self.theta==0:
                P = self.implicit(P,T,J,Act,Q,G)
            elif self.theta==1:
                P = self.explicit(P,T,J,Act,Q,G)
            else:
                P = self.mixed(P,T,J,Act,Q,G,self.theta)

            print(f"{k} time step is complete...")
            
            self._pressure[:,k] = P

            P = P.reshape((-1,1))

    @staticmethod
    def implicit(P,T,J,Act,Q,G):
        """Implicit solution of one-phase flow."""
        return linalg.spsolve(T+J+Act,csr.dot(Act,P)+Q+G)

    @staticmethod
    def mixed(P,T,J,Act,Q,G,theta):
        """Mixed solution of one-phase flow."""

        LHS = (1-theta)(T+J)+Act
        RHS = csr.dot(Act-theta*(T+J),P)+Q+G

        return linalg.spsolve(LHS,RHS)

    @staticmethod
    def explicit(P,T,J,Act,Q,G):
        """Explicit solution of one-phase flow."""
        return P+linalg.spsolve(Act,csr.dot(-(T+J),P)+Q+G)

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
