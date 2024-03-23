from dataclasses import dataclass

import numpy

class Steady():

    def __init__(self,rock:Rock,well:Well,fluid:Fluid):
        """Initializes steady state constant rate radial solution of
        diffusivity equation.

        """

        self._well  = well
        self._rock  = rock
        self._fluid = fluid

    def __call__(self,radii):
        """Returns reservoir pressure values for a constant flow rate
        solution of slightly compressible fluid flow through porous media
        at steady state conditions."""

        trans = self.trans(
            self.rock._perm,
            self.rock._height,
            self.fluid._visc,
            self.fluid._fvf,
            radii,
            self.well._radius
            )

        if self.fluid._comp==0:
            return self.well._rate/trans+self.well._press

        return numpy.log(self.well._rate*self.fluid._comp/trans+1)/self.fluid._comp+self.well._press

    @staticmethod
    def trans(k,h,mu,fvf,re,rw):
        return (2*numpy.pi*k*h)/(mu*fvf*numpy.log(re/rw))

    @property
    def well(self):
        return self._well

    @property
    def rock(self):
        return self._rock

    @property
    def fluid(self):
        return self._fluid
    
    
    
    
