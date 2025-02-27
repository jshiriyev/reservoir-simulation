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