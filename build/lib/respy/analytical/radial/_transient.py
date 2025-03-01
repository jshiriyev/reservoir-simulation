import logging

import numpy

from scipy import special

from ._reservoir import Reservoir, Solution

class Transient(Reservoir):
    """Line source solution based on exponential integral."""

    def __init__(self,*args,**kwargs):

        super().__init__(*args,**kwargs)

        self.tmin = None
        self.tmax = None

    @property
    def tmin(self):
        return self._tmin/(24*60*60)
    
    @tmin.setter
    def tmin(self,value):
        """setting mimimum time limit because of the wellbore size"""
        self._tmin = 100*self.well._radius**2/self._hdiff

    @property
    def tmax(self):
        return self._tmax/(24*60*60)

    @tmax.setter
    def tmax(self,tmax=None):
        """setting maximum time limit because of the external flow radius"""
        self._tmax = 0.25*self._radius**2/self._hdiff

    def _time_correction(self,values):

        bound_internal = values>=self.tmin
        bound_external = values<=self.tmax

        if numpy.any(~bound_internal):
            logging.warning("Not all times satisfy the early time limits!")

        if numpy.any(~bound_external):
            logging.warning("Not all times satisfy the late time limits!")

        valids = numpy.logical_and(bound_internal,bound_external)

        return numpy.where(valids,values,numpy.nan)

    def solve(self,times,points):

        times = self._time_correction(times)

        result = Solution(times,points)

        expis = special.expi(-(result._points**2)/(4*self._hdiff*result._times))

        deltap = self._pterm*(-1/2*expis+self.well._skin)

        result._press = self._pinit-deltap

        return result