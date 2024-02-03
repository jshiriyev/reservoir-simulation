import numpy

from scipy.special import erfc

class OnePhase():

    def __init__(self,length:float,k:float,phi:float,mu:float,ct:float):
        """
        length  : length of the porous media, ft

        k       : permeability, milli-Darcy
        phi     : porosity, non-dimensional
        mu      : viscosity, centi-poise
        ct      : total compressibility, 1/psi

        """

        self.length = length*0.3048

        self.eta = (k*9.869233e-16)/(phi*(ct/6894.76)*(mu*1e-3))

    def solve(self,time:float,ngrids:int,pinit:float,pleft:float,pright=None,noflux=False):
        """
        time    : when to calculate pressures, in hours
        ngrids  : number of pressure calculation points
        pinit   : initial pressure, psi
        pleft   : upstream pressure, psi
        pright  : downstream pressure, psi
        noflux  : downstream pressure boundary conditions
                  True, no-flux boundary condition is implementted
                  False, constant pressure boundary condition is implemented

        returns x (ft) locations where pressures are calculated and pressure values (psi).

        """

        time *= 3600

        pinit *= 6894.76

        pleft *= 6894.76

        xaxis = numpy.arange(self.length/ngrids/2,self.length,self.length/ngrids)
        
        X = self.tondimx(xaxis,self.length)

        T = self.tondimtime(time,self.eta,self.length)
        
        if not noflux: # constant downstream pressure case

            pright = 101325 if pright is None else pright*6894.76

            ndimPresInit = self.tondimpressure(pinit,pleft,pright)
            ndimPressure = self.solve_ndimpressure_pconst(X,T,ndimPresInit)

            dimPressure = self.todimpressure(ndimPressure,pleft,pright)

            return xaxis/0.3048,dimPressure/6894.76

        ndimPressure = self.solve_ndimpressure_noflux(X,T)
        dimPressure = self.todimpressure(ndimPressure,pleft,pinit)

        return xaxis/0.3048,dimPressure/6894.76

    @staticmethod
    def solve_ndimpressure_pconst(ndimx:numpy.ndarray,ndimtime:float,ndimPressInit:float,terms:int=100):
        """Calculates the non-dimensional pressure for constant pressure
        boundary conditions on both sides.

        ndimx           : non-dimensional distance
        ndimtime        : non-dimensional time
        ndimPressInit   : non-dimensional pressure at initial pressure conditions
        terms           : number of terms after which to terminate the summation

        """

        n,psum1,psum2 = 1,numpy.zeros(ndimx.shape),numpy.zeros(ndimx.shape)

        while n<=terms:

            sinterm1 = numpy.sin(n*numpy.pi*ndimx)
            sinterm2 = numpy.sin(n*numpy.pi*(ndimx-1))

            expterm = numpy.exp(-n**2*numpy.pi**2*ndimtime)
            
            psum1 += (-1)**n/n*sinterm1*expterm
            psum2 += (-1)**n/n*sinterm2*expterm

            n += 1

        return 1-ndimx-ndimPressInit*2/numpy.pi*psum1-(1-ndimPressInit)*2/numpy.pi*psum2

    @staticmethod
    def solve_ndimpressure_noflux(ndimx:numpy.ndarray,ndimtime:float,terms:int=100):
        """Calculates the non-dimensional pressure for constant pressure at upstream
        and no-flux at downstream.

        ndimx      : non-dimensional x
        ndimtime   : non-dimensional time
        terms      : number of terms after which to terminate the summation

        """

        n,psum = 0,numpy.zeros(ndimx.shape)

        while n<terms:
            
            expterm = numpy.exp(-(2*n+1)**2/4*numpy.pi**2*ndimtime)
            sinterm = numpy.sin((2*n+1)*numpy.pi/2*ndimx)
            
            # psum += (-1)**n/(2*n+1)*expterm*costerm
            psum += 1/(2*n+1)*expterm*sinterm
            # psum += 4/((2*n+1)*numpy.pi)*expterm*costerm

            n += 1

        return 1-4/numpy.pi*psum
    
    @staticmethod
    def tondimx(dimx,length):
        """Converts dimensional x to non-dimensional x."""
        return dimx/length

    @staticmethod
    def todimx(ndimx,length):
        """Converts non-dimensional x to dimensional x."""
        return ndimx*length

    @staticmethod
    def tondimtime(dimtime,eta,length):
        """Converts dimensional time to non-dimensional time."""
        return eta*dimtime/length**2

    @staticmethod
    def todimtime(ndimtime,eta,length):
        """Converts non-dimensional time to dimensional time."""
        return ndimtime*length**2/eta

    @staticmethod
    def tondimpressure(dimpressure,pb1,pb2):
        """Converts dimensional pressure to non-dimensional pressure."""
        return (dimpressure-pb2)/(pb1-pb2)

    @staticmethod
    def todimpressure(ndimpressure,pb1,pb2):
        """Converts non-dimensional pressure to dimensional pressure."""
        return ndimpressure*(pb1-pb2)+pb2

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

        data = numpy.loadtxt('Lecture_08_dispersion.txt',skiprows=1)

        time = data[:,0]*3600   # seconds

        cc0_c = concentration(XX,v,L,time)

        cc0 = data[:,1]         # dimensionless

        return numpy.sum((cc0-cc0_c)**2)

    def concentration(XX,v,L,time):

        DL = 10**(-XX)

        term1 = erfc((L-v*time)/(2*numpy.sqrt(DL*time)))
        term2 = erfc((L+v*time)/(2*numpy.sqrt(DL*time)))

        if not numpy.any(term2):
            return 1/2*(term1)
        else:
            term3 = numpy.exp(v*L/DL)
            return 1/2*(term1+term2*term3)