import logging

import numpy

from ._reservoir import Reservoir, Boundary, Solution

class Pseudo(Reservoir):

    gamma = 0.5772

    GEOMETRY = {
        "circle": Boundary(31.62,0.1,0.06,0.1),
        "triangle": Boundary(27.6,0.2,0.07,0.09),
        "square": Boundary(30.8828,0.1,0.05,0.09),
        "hexagon": Boundary(31.6,0.1,0.06,0.1),
        }

    def __init__(self,shape,*args,**kwargs):
        # There can be two slightly compressible fluids where the
        # second one is at irreducible saturation, not mobile
        super().__init__(*args,**kwargs)

        self.shape = shape

        self.tmin = None

        self.vpore = None

    @property
    def shape(self):
        return self._shape

    @shape.setter
    def shape(self,value):
        self._shape = value
    
    @property
    def tmin(self):
        return self._tmin/(24*60*60)

    @tmin.setter
    def tmin(self,value):
        self._tmin = self.GEOMETRY[self.shape].time_pss_accurate/(self._hdiff/self._area)

    def tDA(self,values):
        return (self._hdiff/self._area)*values*(24*60*60)

    @property
    def vpore(self):
        return self._vpore/(0.3048**3)
    
    @vpore.setter
    def vpore(self,value):
        self._vpore = self._area*self._height*self.rrock._poro

    def _time_correction(self,values):

        boundary = values>=self.tmin

        if numpy.any(~boundary):
            logging.warning("Not all times satisfy the early time limits!")

        return numpy.where(boundary,values,numpy.nan)

    def solve(self,times,points):

        times = self._time_correction(times)

        result = Solution(times,points)

        CA = self.GEOMETRY[self.shape].factor

        inner = (4*self._area)/(numpy.exp(self.gamma)*CA*self.well._radius**2)

        deltap1 = self._pterm*(1/2*numpy.log(inner)+self.well._skin)

        deltap2 = (self.well._cond*self.fluid._fvf)/(self._vpore*self._tcomp)*result._times

        result._press = self._pinit-deltap1-deltap2

        return result
