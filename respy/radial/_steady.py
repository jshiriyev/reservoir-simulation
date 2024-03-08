class Steady():

    def incompressible(self,rvals,rwell,permeability,height,viscosity):
        """This is a constant flow rate solution of
    incompressible fluid flow through porous media."""

        coeff = (2*numpy.pi*permeability*height)/viscosity

        return self.flowRate/coeff*numpy.log(rvals/rwell)

    def ideal_liquid(self,rvals,rwell,permeability,height,viscosity):
        """This is a constant flow rate solution of
    slightly compressible (ideal liquid) fluid flow through porous media
    at steady state conditions."""

        coeff = (2*numpy.pi*permeability*height)/viscosity

        coeff = self.flowRate/coeff*self.beta*numpy.log(rvals/rwell)

        return numpy.log(coeff+1)/self.beta

    def real_gas(self):

        pass
