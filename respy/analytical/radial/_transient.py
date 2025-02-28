import numpy

from scipy import special

from ._reservoir import Reservoir

class Transient(Reservoir):
    """Line source solution based on exponential integral"""

    def __init__(self,rate,rw,*args,**kwargs):

        super().__init__(*args,**kwargs)

        self.rate = rate

        self.rw = rw

        self.tmin = None
        self.tmax = None

    @property
    def rate(self):
        return self._rate
    
    @rate.setter
    def rate(self,value):
        """setting mimimum time limit because of the wellbore size"""
        self._rate = value

    @property
    def rw(self):
        return self._rw
    
    @rw.setter
    def rw(self,value):
        self._rw = value

    @property
    def tmin(self):
        return self._tmin
    
    @tmin.setter
    def tmin(self,value):
        """setting mimimum time limit because of the wellbore size"""
        self._tmin = 100*self._rw**2/self._diffusivity

    @property
    def tmax(self):
        return self._tmax

    @tmax.setter
    def tmax(self,tmax=None):
        """setting maximum time limit because of the external flow radius"""
        self._tmax = 0.25*self._radius**2/self._diffusivity

    @property
    def times(self):
        return self._times

    @times.setter
    def times(self,values):

        bound_int = values>=self.tmin
        bound_ext = values<=self.tmax

        if numpy.any(~bound_int):
            raise Warning("Not all times satisfy the early time limits!")

        if numpy.any(~bound_ext):
            raise Warning("Not all times satisfy the late time limits!")

        validtimes = numpy.logical_and(bound_int,bound_ext)

        values = values[validtimes]

        self._times = values.reshape((1,-1))

    @property
    def observers(self):
        return self._observers

    @observers.setter
    def observers(self,values):
        self._observers = numpy.ravel(values).reshape((-1,1))

    def solve(self,times,observers,pinit=None):

        constant = (self._rate*self.fluid._visc)/(2*numpy.pi*self.rrock._perm*self._height)
        expivals = special.expi(-(self._observers**2)/(4*self._diffusivity*self._times))

        deltap = -1/2*constant*expivals

        if pinit is None:
            return deltap

        return pinit-deltap

