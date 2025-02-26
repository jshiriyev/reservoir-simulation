class Oil():

    """
    12. properties of crude oil systems
    13. crude oil gravity
    14. specific gravity of the solution gas
    15. gas solubility
    16. bubble-point pressure
    17. oil formation volume factor
    18. isothermal compressibility coefficient of crude oil
    19. oil formation volume factor for undersaturated oils
    20. crude oil density
    21. crude oil viscosity
    22. methods of calculating viscosity of the dead oil
    23. methods of calculating the saturated oil viscosity
    24. methods of calculating the viscosity of the undersaturated oil
    25. surface/interfacial tension
    """

    def __init__(self):

        pass

    def get_gravity(self):
        
        pass

    def get_solution_gas_gravity(self):
        
        pass

    def get_gas_solubility(self):
        
        pass

    def get_bpp(self):
        
        pass

    def get_fvf(self):
        
        pass

    def get_isothermal_compressibility(self):
        
        pass

    def get_density(self):
        
        pass

    def get_viscosity(self):
        
        pass

def Rst():
    """Calculates the stock-tank GOR. It should not be used
    if the separator temperature is >140F."""

    A1 =  0.3818
    A2 = -5.506
    A3 =  2.902
    A4 =  1.327
    A5 = -0.7355

    logrst = A1+A2*numpy.log(gamma_oil)+A3*numpy.log(gamma_gas_sp)+\
                A4*numpy.log(Psp)+A5*numpy.log(Tsp)

    return numpy.exp(logrst)

def pbubble(temp):
    """
    Calculates the bubblepoint pressure of the oil at reservoir conditions

    Rs          : solution GOR including stock-tank vent gas
    gamma_gas   : specific gravity of the separator gas

    The equations are valid to 325F.
    """

    Cpb = (Rs/gamma_gas)**0.83*10**(0.00091*temp-0.0125*gamma_API)

    return 18.2*(Cpb-1.4)

def gor(pres,temp):
    """
    Calculates the solution GOR

    pres        : pressure values below bubble point pressure
    temp        : temperature
    gamma_gas   : specific gravity of the separator gas
    """

    Cpb = pres/18.2+1.4

    return (Cpb*10**(-0.00091*temp+0.0125*gamma_API))**(1/0.83)*gamma_gas

def fvf(pres,temp,Rs,gamma_gas,gamma_oil):

    CBob = Rs*(gamma_gas/gamma_oil)**0.5+1.25*temp

    Bob = 0.9759+12e-5*CBob**1.2

    return Bob*numpy.exp(co*(pbubble-pres))

def rho(pres):

    rhor = (rhoSTO+0.01357*Rs*gamma_gas)/fvf

    return rhob*numpy.exp(comp*(pres-pbubble))

def comp():

    A1 =  -1_433.0
    A2 =       5.0
    A3 =      17.2
    A4 =  -1_180.0
    A5 =      12.61
    A6 = 100_000.0

    return (A1+A2*Rs+A3*temp+A4*gamma_gas+A5*gamma_API)/(A6*pres)

def visc():

    pass