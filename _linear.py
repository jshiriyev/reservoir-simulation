import numpy

from scipy.special import erfc

class Steady():

	def __init__(self,flowRate):

        self.flowRate = flowRate

    def incompressible(self,xvals,permeability,width,height,viscosity):

        coeff = (permeability*width*height)/viscosity

        return self.flowRate/coeff*xvals

class CoreFlow():

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

class CoreFlow_v2(): # will be DEPRECIATED

    @staticmethod
    def solve(length,time,pinit,pright,eta,ngrids):
        """
        lenght  : length of core sample
        time    : time at which to calculate pressure
        pinit   : initial pressure
        pright  : constant pressure boundary implemented at right hand side
        eta     : hydraulic diffusivity
        ngrids  : number of pressure calculation points
        """

        x = np.arange(length/ngrids/2,length,length/ngrids)
        
        X = CoreFlow.tondimx(x,length)

        T = CoreFlow.tondimtime(time,eta,length)
        
        P = CoreFlow.ndimpressure(X,T)
        
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
    def ndimpressure(ndimx:np.ndarray,ndimtime:float,terms:int=100):
        """Calculates the non-dimensional pressure given with equation 5 on the word document.

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


class BuckleyLeverett():
    
    """
    based on Buckley-Leverett solution
    """

    def __init__(self,Sor,Swr,muo,muw):

        self.Sor = Sor
        self.Swr = Swr
        self.muo = muo
        self.muw = muw

    def k_model(self):

        N = 1000

        self.Sw = numpy.linspace(self.Swr,1-self.Sor,N)

        self.kro = 2*(1-self.Sw-self.Sor)**2
        self.krw = (self.Sw-self.Swr)**3

    def coreymodel(self,koro,korw,m,n):

        N = 1000

        self.Sw = numpy.linspace(self.Swr,1-self.Sor,N)

        S = (self.Sw-self.Swr)/(1-self.Swr-self.Sor)

        self.kro = koro*(1-S)**m
        self.krw = korw*S**n

        ## end-point mobility ratio calculation
        self.Mo = (korw/self.muw)/(koro/self.muo)

    def fractionalflow(self):

        self.fw = (self.krw*self.muo)/(self.krw*self.muo+self.kro*self.muw)
        
        N = self.fw.size

        one = numpy.ones(N-1)

        idx = numpy.array(list(range(N)))

        row = numpy.concatenate(((idx[0],idx[-1]),idx[:-1],idx[1:]))
        col = numpy.concatenate(((idx[0],idx[-1]),idx[1:],idx[:-1]))

        val = numpy.concatenate(((-1,1),one,-one))

        G = csr((val,(row,col)),shape=(N,N))

        fw_diff = G*self.fw
        Sw_diff = G*self.Sw

        self.fw_der = fw_diff/Sw_diff
        
    def shockfront(self,Swi):

        self.Swi = Swi

        IC = self.Sw>=self.Swi

        ## loosing some data in fw, Sw and fw_der for the saturations below the initial value
        self.fw_IC = self.fw[IC]
        self.Sw_IC = self.Sw[IC]
        ## fw_der_IC is the fw_der corrected for the shock front as well
        self.fw_der_IC = self.fw_der[IC]
        
        self.fwi = self.fw_IC[0]
        
        idx = numpy.argmax(self.fw_der_IC)

        fw_dh = self.fw_IC[idx:]
        Sw_dh = self.Sw_IC[idx:]

        pseudo_fw_der = (fw_dh-self.fwi)/(Sw_dh-self.Swi)

        self.fwf = fw_dh[numpy.argmin(numpy.abs(self.fw_der_IC[idx:]-pseudo_fw_der))]
        self.Swf = Sw_dh[numpy.argmin(numpy.abs(self.fw_der_IC[idx:]-pseudo_fw_der))]

        self.Sw_avg = self.Swi+(1-self.fwi)/(self.fwf-self.fwi)*(self.Swf-self.Swi)

        fw_der_corrected = numpy.empty_like(self.fw_der_IC)

        fw_der_corrected[self.Sw_IC>=self.Swf] = self.fw_der_IC[self.Sw_IC>=self.Swf]
        fw_der_corrected[self.Sw_IC<self.Swf] = self.fw_der_IC[self.Sw_IC==self.Swf]

        self.fw_der_IC = fw_der_corrected

        self.fwf_der = self.fw_der_IC[0]

    def production(self,q,A,L,phi):

        """
        given that L is "ft", A is "ft2", and q is "ft3/day", tp will be in days.
        phi is porosity and is a dimensionless value.
        calculated volumes are normalized with respect to pore volume Vp,
        both produced oil Np and injected water Wi
        """

        v = q/(phi*A)

        self.v = v

        self.Vp = A*L*phi
        
        self.tbt = L/(v*self.fwf_der)

        self.Nbt = v/L*(1-self.fwi)*self.tbt

        self.tp = L/(v*self.fw_der_IC)
        self.Np = (1-self.fw_IC)/(self.fw_der_IC)+self.Sw_IC-self.Swi

        idx = self.tp<=self.tbt

        N = numpy.count_nonzero(idx)

        self.tp[idx] = numpy.linspace(0,self.tbt,N)
        self.Np[idx] = v/L*(1-self.fwi)*self.tp[idx]
        
        self.Wi = v/L*self.tp

    def profile(self,q,A,phi,t):

        """
        time is in days assuming that L is "ft", A is "ft2", and q is "ft3/day".
        phi is porosity and is a dimensionless value.
        calculated x_Sw for a saturation profile will be in "ft".
        """

        v = q/(phi*A)        

        self.x_Sw = v*t*self.fw_der_IC

class MiscibleDisplacement():

    def obj(XX,v,L):

        data = np.loadtxt('Lecture_08_dispersion.txt',skiprows=1)

        time = data[:,0]*3600   # seconds

        cc0_c = concentration(XX,v,L,time)

        cc0 = data[:,1]         # dimensionless

        return np.sum((cc0-cc0_c)**2)

    def concentration(XX,v,L,time):

        DL = 10**(-XX)

        term1 = erfc((L-v*time)/(2*np.sqrt(DL*time)))
        term2 = erfc((L+v*time)/(2*np.sqrt(DL*time)))

        if not np.any(term2):
            return 1/2*(term1)
        else:
            term3 = np.exp(v*L/DL)
            return 1/2*(term1+term2*term3)