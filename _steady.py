import numpy

class Incompressible():
    """This is a constant flow rate solution of
    incompressible fluid flow through porous media."""

    def __init__(self,flowRate):

        self.flowRate = flowRate

    def linear(self,xvals,permeability,width,height,viscosity):

        coeff = (permeability*width*height)/viscosity

        return self.flowRate/coeff*xvals

    def radial(self,rvals,rwell,permeability,height,viscosity):

        coeff = (2*numpy.pi*permeability*height)/viscosity

        return self.flowRate/coeff*numpy.log(rvals/rwell)

class IdealLiquid():
    """This is a constant flow rate solution of
    slightly compressible (ideal liquid) fluid flow through porous media
    at steady state conditions."""

    def __init__(self,flowRate,beta):
        """
        beta  : fluid compressbility constant
        """

        self.flowRate = flowRate

        self.beta = beta

    def linear(self):

        pass

    def radial(self,rvals,rwell,permeability,height,viscosity):

        coeff = (2*numpy.pi*permeability*height)/viscosity

        coeff = self.flowRate/coeff*self.beta*numpy.log(rvals/rwell)

        return numpy.log(coeff+1)/self.beta

class RealGas():

    def __init__(self,flowRate):

        pass