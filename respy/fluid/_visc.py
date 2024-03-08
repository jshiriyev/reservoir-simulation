from inspect import isfunction

import numpy

class CKB():
    """Carr Kobayashi Burrows Method"""

    a0  = -2.46211820
    a1  =  2.970547414
    a2  = -2.86264054e-1
    a3  =  8.05420522e-3
    a4  =  2.80860949
    a5  = -3.49803305
    a6  =  3.60373020e-1
    a7  = -1.044324e-2
    a8  = -7.93385648e-1
    a9  =  1.39643306
    a10 = -1.49144925e-1
    a11 =  4.41015512e-3
    a12 =  8.39387178e-2
    a13 = -1.86408848e-1
    a14 =  2.03367881e-2
    a15 = -6.09579263e-4

    def __init__(self,spgr,temp,N2=0,CO2=0,H2S=0):

        self._spgr = spgr
        self._temp = temp

        mu0 = self.uncorrected(spgr,temp)
        muN = self.N2corr(spgr,N2)
        muC = self.CO2corr(spgr,CO2)
        muS = self.H2Scorr(spgr,H2S)

        mu1 = self.corrected(mu0,muN,muC,muS)

        self.mu1 = mu1

    @property
    def spgr(self):
        return self._spgr

    @property
    def temp(self):
        return self._temp

    def __call__(self,pred,tred):

        return self.mu1*self.ratio(pred,tred)

    def ratio(self,pred,tred):
        """Dempsey (1965) viscosity ratio, mu/mu1
        
        pred    : reduced pressure
        tred    : reduced temperature
        """

        b0 = self.a0,  self.a1,  self.a2,  self.a3
        b1 = self.a4,  self.a5,  self.a6,  self.a7
        b2 = self.a8,  self.a9,  self.a10, self.a11
        b3 = self.a12, self.a13, self.a14, self.a15

        A0 = self._polynomial(*b0,pred)
        A1 = self._polynomial(*b1,pred)
        A2 = self._polynomial(*b2,pred)
        A3 = self._polynomial(*b3,pred)

        poly = self._polynomial(A0,A1,A2,A3,tred)

        return numpy.exp(poly)/tred

    @staticmethod
    def _polynomial(c0,c1,c2,c3,x):
        return c0+c1*x+c2*x**2+c3*x**3

    @staticmethod
    def corrected(mu0,dmuN2,dmuCO2,dmuH2S):
        """Corrected gas viscosity at one atmospheric
        pressure and reservoir temperature, cp"""
        return mu0+dmuN2+dmuCO2+dmuH2S

    @staticmethod
    def uncorrected(spgr,T):
        """Uncorrected gas viscosity, cp

        T   : Temperature in Fahrenheits
        """
        return 8.188e-3-6.15e-3*numpy.log(spgr)\
            +(1.709e-5-2.062e-6*spgr)*T

    @staticmethod
    def N2corr(spgr,mfract):
        """Visocsity correction due to the presence of N2"""
        return mfract*(8.48e-3*numpy.log(spgr)+9.59e-3)

    @staticmethod
    def CO2corr(spgr,mfract):
        """Visocsity correction due to the presence of CO2"""
        return mfract*(9.08e-3*numpy.log(spgr)+6.24e-3)

    @staticmethod
    def H2Scorr(spgr,mfract):
        """Visocsity correction due to the presence of H2S"""
        return mfract*(8.49e-3*numpy.log(spgr)+3.73e-3)

class LGE():
    """Lee Gonzalez Eakin Method"""
    def __init__(self,G,T):

        self._G = G
        self._T = T*(5./9)

        self.A = self.get_A()
        self.B = self.get_B()
        self.C = self.get_C()

    @property
    def G(self):
        return self._G

    @property
    def M(self):
        return 28.964*self.G

    @property
    def _M(self):
        return 28.964*self._G/1000

    @property
    def T(self):
        return self._T*(9./5)

    def get_A(self):
        return (9.379+0.01607*self.M)*self.T**1.5/(209.2+19.26*self.M+self.T)

    def get_B(self):
        return 3.448+986.4/self.T+0.01009*self.M

    def get_C(self):
        return 2.447-0.2224*self.B
    
    def __call__(self,P,Z):
        """Returns Gas Viscosity in cp"""

        rho = (P*self.M)/(Z*10.731577089016*self.T)
        
        return self.A*numpy.exp(self.B*rho**self.C)/10000

if __name__ == "__main__":

    a = LGE(0.8,200)

    print(a._T)

    print(a(500,0.9))

    def L(x):
        return x**2

    c = 2

    print(callable(c))