class Transient():

    # Line source solution based on exponential integral

    def __init__(self,flow_rate):

        super().__init__(flow_rate)

    def set_tmin(self):

        # setting mimimum time limit because of the wellbore size

        self.tmin = 100*self.Well.radii[0]**2/self.diffusivity

        return self.tmin

    def set_tmax(self,tmax=None):

        # setting maximum time limit because of the external flow radius
        # if the tmax is provided, new external radius will be calculated 

        if tmax is None:
            tmax = 0.25*self.PorRock.radii[0]**2/self.diffusivity
        else:
            tmax_radius = float(np.sqrt(tmax*self.diffusivity/0.25))
            self.PorRock.set_radii(tmax_radius)

        self.tmax = tmax

        return self.tmax

    def set_times(self,times=None,timespace="linear"):

        if not hasattr(self,"tmin"):
            self.set_tmin()

        if not hasattr(self,"tmax"):
            self.set_tmax()

        if times is None:
            if timespace=="linear":
                times = np.linspace(self.tmin,self.tmax,1000)
            elif timespace=="log":
                times = np.logspace(np.log10(self.tmin),np.log10(self.tmax),1000)
        else:
            bound_int = times>=self.tmin
            bound_ext = times<=self.tmax

            if np.any(~bound_int):
                raise Warning("Not all times satisfy the early time limits!")

            if np.any(~bound_ext):
                raise Warning("Not all times satisfy the late time limits!")

            validtimes = np.logical_and(bound_int,bound_ext)

            times = times[validtimes]

        self.times = times.reshape((1,-1))

    def set_observers(self,observers=None,number=50):

        if observers is not None:
            self.observers = observers
        else:
            inner = self.Well.radii[0]
            if self.shape == "circular":
                outer = self.PorRock.radii[0]
            elif self.shape == "square":
                outer = self.PorRock.lengths[0]
            self.observers = np.linspace(inner,outer,number)

        self.observers = self.observers.reshape((-1,1))

    def solve(self):

        if not hasattr(self,"times"):
            self.set_times()

        rateNumer = self.Well.limits*self.Fluids.viscosity[0]
        rateDenom = self.PorRock.permeability*self.PorRock.thickness

        constRate = (rateNumer)/(2*np.pi*rateDenom)

        Ei = expi(-(self.observers**2)/(4*self.diffusivity*self.times))

        self.deltap = -1/2*constRate*Ei

        self.pressure = self.pressure0-self.deltap

