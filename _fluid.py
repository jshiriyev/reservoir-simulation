from dataclasses import dataclass

import math

from scipy import optimize

@dataclass(frozen=True)
class Fluid:
    
    density		: float
    viscosity   : float

    def __post_init__(self):

        object.__setattr__(self,'_density',self.density*16.0185)
        object.__setattr__(self,'_viscosity',self.viscosity*0.001)

class Gas():

    A1  = 0.3265
    A2  = -1.0700
    A3  = -0.5339
    A4  = 0.01569
    A5  = -0.05165
    A6  = 0.5475
    A7  = -0.7361
    A8  = 0.1844
    A9  = 0.1056
    A10 = 0.6134
    A11 = 0.7210

    def __init__(self,grav:float):
        """
        grav    : gas specific gravity
        """
        self._grav = grav

    @property
    def grav(self):
        return self._grav
    
    @property
    def molmass(self):
        """Molecular Mass, lb/lbmol"""
        return 28.964*self._grav

    def _molmass(self):
        """Molecular Mass, kg/mol"""
        return 28.964*self._grav/1000

    @property
    def pcrit(self):
        """Function to Calculate Gas Critical Pressure in psia
        for high molecular weigth reservoir gases"""
        return 756.8-131*self._grav-3.6*self._grav**2

    @property
    def tcrit(self):
        """Function to Calculate Gas Critical Temperature in °R
        for high molecular weigth reservoir gases"""
        return 169.2+349.5*self._grav-74*self._grav**2

    def update(self,P,T):
        """
        P       Pressure, psia
        T       Temperature, °R
        """
        self.pres = P
        self.temp = T

        self.red_pres = P/self.pcrit
        self.red_temp = T/self.tcrit

        self._zfact = self.get_zfact(self.red_pres,self.red_temp)

    @property
    def zfact(self):
        return self._zfact
    
    @staticmethod
    def get_zfact(Pr,Tr):
        """Function to Calculate Gas Compressibility Factor
        Pr         reduced pressure
        Tr         reduced temperatue
        """
        # rhor = 0.27*(Pr)/(Z*Tr)
        # 
        # 1
        # +(A1+A2/Tr+A3/Tr**3+A4/Tr**4+A5/Tr**5)*rhor
        # +(A6+A7/Tr+A8/Tr**2)*rhor**2
        # -A9*(A7/Tr+A8/Tr**2)*rhor**5
        # +A10*(1+A11*rhor**2)*(rhor**2/Tr**3)*math.exp(-A11*rhor**2)

        a = 1.39*(Tr-0.92)**0.5-0.36*Tr-0.101
        b = (0.62-0.23*Tr)*Pr+(0.066/(Tr-0.86)-0.037)*Pr**2+0.32*Pr**6/(10**(9*(Tr-1)))
        c = (0.132-0.32*math.log(Tr)/math.log(10))
        d = 10**(0.3106-0.49*Tr+0.1824*Tr**2)

        return a+(1-a)*math.exp(-b)+c*Pr**d

    @property
    def rho(self):
        """Returns density value in lb/ft3"""
        return (self.pres*self.molmass)/(self.zfact*10.731577089016*self.temp)

    @property
    def fvf(self):
        """Function to Calculate Gas Formation Volume Factor in ft3/scf"""
        return 0.02827*self.zfact*self.temp/self.pres

    @property
    def comp(self):
        """This method should not be used at Tr<1.4 for 0.4<Pr<3.0"""
        Pc,Tc = self.pcrit,self.tcrit

        Pr,Tr = self.red_pres,self.red_temp

        Z = self.zfact

        rhor = 0.27*(Pr)/(Z*Tr)

        prime = self.A1+self.A2/Tr+self.A3/Tr**3+self.A4/Tr**4
        prime += self.A5/Tr**5+2*rhor*(self.A6+self.A7/Tr+self.A8/Tr**2)
        prime -= 5*rhor**4*self.A9*(self.A7/Tr+self.A8/Tr**2)
        prime += 2*self.A10*rhor/Tr**3*(1+self.A11*rhor**2-self.A11**2*rhor**4)*math.exp(-self.A11*rhor**2)

        return 1/self.pres-0.27/(Z**2*Tr*Pc)*(prime/(1+rhor/Z*prime))

    @property
    def visc(self):
        """Returns Gas Viscosity in cp"""
        A = (9.379+0.01607*self.molmass)*T**1.5/(209.2+19.26*self.molmass+self.temp)
        B = 3.448+986.4/T+0.01009*self.molmass
        C = 2.447-0.2224*B

        rho = (self.press*self.molmass)/(self._zfact*10.731577089016*self.temp)
        
        return A*math.exp(B*rho**C)/10000
    