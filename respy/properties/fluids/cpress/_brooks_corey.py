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

    def inverse_drainage(self,pressure):
        """Calculates saturation values from capillary pressure for
        drainage scenario, inverse calculations.
        
        pressure     : capillary pressure values

        """

        pressure = numpy.asarray(pressure)

        Sstar = (self.entry/pressure)**(self.lamda)

        return self.irreducible+Sstar*(1-self.irreducible)

    def imbibition(self,saturation):
        """Calculates capillary pressure values for imbibition scenario,
        forward calculations."""

        saturation = numpy.asarray(saturation)

        Se = (saturation-self.irreducible)/(1-self.irreducible-self.residual)

        return self.entry*(Se**(-1/self.lamda)-1)

    def inverse_imbibition(self,pressure):
        """Calculates saturation values for imbibition scenario,
        inverse calculations.

        pressure     : capillary pressure values

        """

        pressure = numpy.asarray(pressure)

        Se = (pressure/self.entry+1)**(-self.lamda)

        return self.irreducible+Se*(1-self.irreducible-self.residual)