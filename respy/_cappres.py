import numpy

class BrooksCorey():

    def __init__(self,irreducible=0,residual=0,lamda=2,entry=4.5):
        """
        BrooksCorey capillary  pressure model

        irreducible     : Irreducible saturation
        residual        : Residual saturation

        lamda           : Empirical model constant
        entry           : Capillary entry pressure

        """

        self.irreducible = irreducible
        self.residual = residual

        self.lamda = lamda
        self.entry = entry

    def drainage(self,saturation):
        """Calculates capillary pressure values for drainage
        scenario, forward calculations."""

        saturation = numpy.asarray(saturation)

        Sstar = (saturation-self.irreducible)/(1-self.irreducible)

        return self.entry*Sstar**(-1/self.lamda)

    def idrainage(self,cappres):
        """Calculates saturation values from capillary pressure for
        drainage scenario, inverse calculations.
        
        cappres     : capillary pressure values

        """

        cappres = numpy.asarray(cappres)

        Sstar = (self.entry/cappres)**(self.lamda)

        return self.irreducible+Sstar*(1-self.irreducible)

    def imbibition(self,saturation):
        """Calculates capillary pressure values for imbibition scenario,
        forward calculations."""

        saturation = numpy.asarray(saturation)

        Se = (saturation-self.irreducible)/(1-self.irreducible-self.residual)

        return self.entry*(Se**(-1/self.lamda)-1)

    def iimbibition(self,cappres):
        """Calculates saturation values for imbibition scenario,
        inverse calculations.

        cappres     : capillary pressure values

        """

        cappres = numpy.asarray(cappres)

        Se = (cappres/self.entry+1)**(-self.lamda)

        return self.irreducible+Se*(1-self.irreducible-self.residual)

class VanGenuchten():

    def __init__(self,irreducible=0,residual=0,n=2,m=2,gamma=0.02):

        self.irreducible = irreducible
        self.residual = residual

        self.n,self.m = n,m

        self.gamma = gamma

    def imbibition(self,sw):

        saturation = (sw-self.irreducible)/(1-self.irreducible-self.residual)

        return 1/self.gamma*(saturation**(-1/self.m)-1)**(1/self.n)

    def iimbibition(self,cappres):
        """Inverse imbibition calculation, saturation from capillary pressure."""

        saturation = ((self.gamma*cappres)**self.n+1)**(-self.m)

        return self.irreducible+saturation*(1-self.irreducible-self.residual)

class JFunction():
    """
    The Leverett J-function is a way to nondimensionalize capillary pressure curves.
    It can applied to determine capillary pressure in a similar rock type
    but different permeability, porosity, or interfacial tension of the fluids.
    """

    @staticmethod
    def direct(pressure,permeability,porosity,interfacial,contact):
        """Calculates J-function from capillary pressure"""

        termPS = numpy.sqrt(permeability/porosity)
        termIT = interfacial*numpy.cos(contact)

        return pressure*termPS/termIT

    @staticmethod
    def inverse(jfunction,permeability,porosity,interfacial,contact):
        """Calculates capillary pressure from J-function"""

        termPS = numpy.sqrt(permeability/porosity)
        termIT = interfacial*numpy.cos(contact)

        return jfunction*termIT/termPS

def ScanCurves(sw,hyst,cappres,epspc=0,residual=0):
    """Initialization of a reservoir is assumed to occur by primary drainage (oleic
    displaces aqueous phase) by migration from a source rock.

    During waterflooding, the capillary pressure transitions from a drainage curve to
    an imbibition curve.

    Capillary scanning curves are utilized in simulators to model this
    transition.

    A capillary scanning curve is a weighted average of the imbibition
    and drainage capillary pressure curves."""

    conds = sw<=hyst

    cappres = numpy.zeros(sw.shape)

    drainage = cappres.drainage(sw)

    imbibition = cappres.imbibition(sw)

    Av = 1-residual-hyst
    Bv = sw-hyst

    fv = (Av+epspc)/(Bv+epspc)*Bv/Av

    cappres[conds] = drainage

    cappres[~conds] = fv*imbibition+(1-fv)*drainage

    return cappres

if __name__ == "__main__":

    cappres = BrooksCorey(0.1,0.4,lamda=2,entry=3.5)

    pcDm = cappres.drainage(0.3)
    pcIm = cappres.imbibition(0.3)

    print(pcDm,pcIm)

    jfunctionD = JFunction.direct(pcDm,100,0.2,20,0)
    jfunctionI = JFunction.direct(pcIm,100,0.2,20,0)

    pcDc = JFunction.inverse(jfunctionD,50,0.15,25,0)
    pcIc = JFunction.inverse(jfunctionI,50,0.15,25,0)

    print(pcDc,pcIc)