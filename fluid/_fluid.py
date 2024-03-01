from dataclasses import dataclass

import numpy

from scipy import optimize

from respy.fluid._zfact import unknown
from respy.fluid._zfact import Hall_Yarborough
from respy.fluid._zfact import Dranchuk_Abu_Kassem
from respy.fluid._zfact import Dranchuk_Purvis_Robinson

"""
 1. Watch Python packaging again and setup parent directory accordingly
 2. Complete the exercise with zfactor, get the graph
 3. Complete all zfactor methods to give z and cg
 4. Exercise for single phase compressible flow analytical solution
 5. Finalize linear module one-phase flow
 6. Benchmark numerical solution for compressible flow
 7. Finalize the single phase flow of all kinds, all non-linearities in mind
 8. Two phase flow of water and oil, no capillary pressure
 9. Finalize Buckley Leverett solution
10. Benchmark numerical solution with Buckley Leverett
11. Effective implementation of upwinding
12. Finalize Relative Permeability curves
13. Two phase flow of water and oil, with capillary pressure
14. Finalize capillary pressure models.
15. Finalize two phase flow
16. Finalize three phase flow
"""

@dataclass(frozen=True)
class Fluid:
    
    rho		: float
    visc    : float
    comp    : float = None
    fvf     : float = None

    def __post_init__(self):

        object.__setattr__(self,'_rho',self.rho*16.0185)
        object.__setattr__(self,'_visc',self.visc*0.001)

        if self.comp is not None:
            object.__setattr__(self,'_comp',self.comp/6894.76)

        if self.fvf is not None:
            object.__setattr__(self,'_fvf',self.fvf)

class Gas():

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

    @property
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

    def __call__(self,P,T,method=None):
        """
        P       Pressure, psia
        T       Temperature, °R
        """

        Pr = P/self.pcrit
        Tr = T/self.tcrit

        Z,Zprime = self.get_zfact(Pr,Tr,method=method)

        rho  = self.get_rho(P,T,Z,M)
        visc = self.get_visc(P,T,Z,M)
        comp = self.get_comp(Pr,Tr,Z,Zprime)/self.pcrit
        fvf  = self.get_fvf(P,T,Z)

        return Fluid(rho,visc,comp,fvf)
    
    @staticmethod
    def get_zfact(Pr,Tr,method=None):
        """Function to Calculate Gas Compressibility Factor"""
        return

    @staticmethod
    def get_rho(P,T,Z,M):
        """Returns density value in lb/ft3"""
        return (P*M)/(Z*10.731577089016*T)

    @staticmethod
    def get_visc(P,T,Z,M):
        """Returns Gas Viscosity in cp"""
        A = (9.379+0.01607*M)*T**1.5/(209.2+19.26*M+T)
        B = 3.448+986.4/T+0.01009*M
        C = 2.447-0.2224*B

        rho = (P*M)/(Z*10.731577089016*T)
        
        return A*numpy.exp(B*rho**C)/10000

    @staticmethod
    def get_comp(Pr,Tr,Z,Zprime):
        """Returns cpr=cg*Pc. This method should not be used at Tr<1.4 for 0.4<Pr<3.0"""
        return 1/(1+(0.27*Pr)/(Z**2*Tr)*Zprime)/Pr

    @staticmethod
    def get_fvf(P,T,Z):
        """Function to Calculate Gas Formation Volume Factor in ft3/scf"""
        return 0.02827*Z*T/P

class Oil():

    def __init__(self):

        pass

    def update(self):

        pass

class Water():

    def __init__(self):

        pass

    def update(self):

        pass