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

"""Critical Temperature and Pressure Calculations"""

def Pcrit(grav):
    """Function to Calculate Gas Critical Pressure in psia
    for high molecular weigth reservoir gases

    grav       gas specific gravity
    """
    return 756.8-131*grav-3.6*grav**2

def Tcrit(grav):
    """Function to Calculate Gas Critical Temperature in °R
    for high molecular weigth reservoir gases

    grav       gas specific gravity
    """
    return 169.2+349.5*grav-74*grav**2

"""Reduced Temperature and Pressure Calculations"""

def Pred(grav,P):
    return P/Pcrit(grav)

def Tred(grav,T):
    return T/Tcrit(grav)

"""Compressibility Factor Calculations"""

def zfact(Pr,Tr):
    """Function to Calculate Gas Compressibility Factor
    Pr         reduced pressure
    Tr         reduced temperatue
    """
    a = 1.39*(Tr-0.92)**0.5-0.36*Tr-0.101
    b = (0.62-0.23*Tr)*Pr+(0.066/(Tr-0.86)-0.037)*Pr**2+0.32*Pr**6/(10**(9*(Tr-1)))
    c = (0.132-0.32*math.log(Tr)/math.log(10))
    d = 10**(0.3106-0.49*Tr+0.1824*Tr**2)
    return a+(1-a)*math.exp(-b)+c*Pr**d

"""Gas Density Calculations"""

def grho(grav,P,T,Z):
    """Returns density value in lb/ft3
    
    grav    Specific gravity

    P       Pressure, psia
    T       Temperature, °R
    Z       Compressibility factor
    """
    M = 28.964*grav # lb/lbmol
    R = 10.731577089016 # (psi⋅ft3)/(lbmol⋅°R)

    return (P*M)/(Z*R*T)

"""Gas Formation Volume Factor Calculations"""

def gfvf(P,T,Z):
    """Function to Calculate Gas Formation Volume Factor in ft3/scf

    P       Pressure, psia
    T       Temperature,°R
    Z       Compressibility factor
    """
    return 0.02827*Z*T/P

"""Gas Compressibility Calculations"""

def gcomp(grav,P,T,Z,zprime,pcritical):

    return (pcritical/P-1/Z*zprime)/pcritical

"""Gas Viscosity Calculations"""

def gvisc(grav,P,T,Z): 
    """Function to Calculate Gas Viscosity in cp
    
    grav    Specific gravity

    P       Pressure, psia
    T       Temperature, °R
    Z       Compressibility factor
    """
    M = 28.964*grav
    X = 3.448+986.4/T+0.01009*M
    Y = 2.447-0.2224*X
    rho = (1.4926/1000)*P*M/Z/T
    K = (9.379+0.01607*M)*T**1.5/(209.2+19.26*M+T)
    return K*math.exp(X*rho**Y)/10000