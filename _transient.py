import numpy as np

from scipy.special import expi

class Linear():

    def __init__(self,length:float,pinit:float):
        """
        length  : length of the porous media, ft
        pinit   : initial pressure, psi
        """

        self.length,self.pinit = length*0.3048,pinit*6894.76

    def set_eta(self,k,phi,mu,ct):
        """
        k   : permeability, milli-Darcy
        phi : porosity, non-dimensional
        mu  : viscosity, centi-poise
        ct  : total compressibility, 1/psi
        """
        k *= 9.869233e-16
        mu *= 1e-3
        ct /= 6894.76

        self.eta = self._eta(k,phi,mu,ct)

    def solve(self,time,ngrids,pright=None,noflux=False):
        """
        time    : in hours
        ngrids  : number of pressure calculation points
        pright  : pressure on the right boundary, psi
        noflux  : if True, no-flux boundary conditions are implementted on RHS.
        """

        time *= 3600

        if not noflux:
            pright = 101325 if pright is None else pright/6894.76

        xaxis,Pressure = self._solve(self.length,time,self.pinit,self.eta,ngrids,pright,noflux)

        return xaxis,Pressure/6894.76

    @staticmethod
    def _solve(length,time,pinit,eta,ngrids,pright,noflux):
        """
        lenght  : length of core sample, meter
        time    : time at which to calculate pressure, second
        pinit   : initial pressure, Pascal
        eta     : hydraulic diffusivity, m2/second
        ngrids  : number of pressure calculation points
        pright  : constant pressure boundary implemented at right hand side
        """

        x = np.arange(length/ngrids/2,length,length/ngrids)
        
        X = SinglePhaseLinear._tondimx(x,length)

        T = SinglePhaseLinear._tondimtime(time,eta,length)
        
        if noflux:
            P = SinglePhaseLinear._ndimpressure_noflux(X,T)
        else:
            P = SinglePhaseLinear._ndimpressure_pconst(X,T)
        
        P = SinglePhaseLinear._todimpressure(P,pinit,pright)

        return x,P

    @staticmethod
    def _eta(k,phi,ct,mu):
        """
        k   : permeability, meter-square
        phi : porosity, non-dimensional
        mu  : viscosity, Pa.second
        ct  : total compressibility, 1/Pascal
        """
        return k/(phi*ct*mu)
    
    @staticmethod
    def _tondimx(dimx,length):
        """Converts dimensional x to non-dimensional x."""
        return dimx/length

    @staticmethod
    def _todimx(ndimx,length):
        """Converts non-dimensional x to dimensional x."""
        return ndimx*length

    @staticmethod
    def _tondimtime(dimtime,eta,length):
        """Converts dimensional time to non-dimensional time."""
        return eta*dimtime/length**2

    @staticmethod
    def _todimtime(ndimtime,eta,length):
        """Converts non-dimensional time to dimensional time."""
        return ndimtime*length**2/eta

    @staticmethod
    def _tondimpressure(dimpressure,pinit,pright):
        """Converts dimensional pressure to non-dimensional pressure."""
        return (dimpressure-pright)/(pinit-pright)

    @staticmethod
    def _todimpressure(ndimpressure,pinit,pright):
        """Converts non-dimensional pressure to dimensional pressure."""
        return ndimpressure*(pinit-pright)+pright

    @staticmethod
    def _ndimpressure_pconst(ndimx:np.ndarray,ndimtime:float,terms:int=100):
        """Calculates the non-dimensional pressure for constant pressure
        boundary conditions on both sides.

        ndimx      : non-dimensional x
        ndimtime   : non-dimensional time
        terms      : number of terms after which to terminate the summation

        """

        psum1 = np.zeros(ndimx.shape)
        psum2 = np.zeros(ndimx.shape)

        for n in range(1,terms):

            fstterm = (-1)**n/n

            sinterm1 = np.sin(n*np.pi*ndimx)
            sinterm2 = np.sin(n*np.pi*(ndimx-1))

            expterm = np.exp(-n**2*np.pi**2*ndimtime)
            
            psum1 += fstterm*sinterm1*expterm
            psum2 += fstterm*sinterm2*expterm

        psum1 = psum1*2/np.pi
        psum2 = psum2*2/np.pi

        pdi = 0.2

        return 1-ndimx-pdi*psum1-(1-pdi)*psum2

    @staticmethod
    def _ndimpressure_noflux(ndimx:np.ndarray,ndimtime:float,terms:int=100):
        """Calculates the non-dimensional pressure for constant pressure on left side
        and no-flux on right side of sample.

        ndimx      : non-dimensional x
        ndimtime   : non-dimensional time
        terms      : number of terms after which to terminate the summation

        """

        psum = np.zeros(ndimx.shape)

        for n in range(terms):

            fstterm = (-1)**n/(2*n+1)
            
            expterm = np.exp(-(2*n+1)**2*np.pi**2*ndimtime/4)
            costerm = np.cos((2*n+1)*np.pi/2*ndimx)
            
            psum += fstterm*expterm*costerm

        return psum*4/np.pi

class Radial():

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