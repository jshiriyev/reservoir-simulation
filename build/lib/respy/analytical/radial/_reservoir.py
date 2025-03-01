from dataclasses import dataclass

import logging

import numpy

class Reservoir():

    def __init__(self,area,height,rrock,fluid,well,pinit=None,tcomp=None):
        """Initialization of base reservoir class for analytical calculations."""
        self.area = area

        self.radius = None
        self.height = height

        self.rrock = rrock
        self.xflow = None
        self.fluid = fluid

        self.well = well

        self.pinit = pinit
        self.tcomp = tcomp

        self.hdiff = None
        self.pterm = None

    @property
    def area(self):
        return self._area/(43560*0.3048**2)

    @area.setter
    def area(self,value):
        self._area = value*43560*0.3048**2

    @property
    def radius(self):
        return self._radius/0.3048

    @radius.setter
    def radius(self,value):
        self._radius = numpy.sqrt(self._area/numpy.pi)

    @property
    def height(self):
        return self._height/0.3048

    @height.setter
    def height(self,value):
        self._height = value*0.3048

    @property
    def xflow(self):
        """Getter for the rock transmissibility in x-direction."""
        return self._xflow/(9.869233e-16*0.3048)

    @xflow.setter
    def xflow(self,value):
        """Setter for the rock transmissibility in x-direction."""
        self._xflow = (self.rrock._xperm*self._height)

    @property
    def pinit(self):
        return None if self._pinit is None else self._pinit/6894.76

    @pinit.setter
    def pinit(self,value):
        self._pinit = None if value is None else value*6894.76
    
    @property
    def tcomp(self):
        """Getter for the total compressibility."""
        return None if self._tcomp is None else self._tcomp*6894.76

    @tcomp.setter
    def tcomp(self,value):
        """Setter for the total compressibility."""
        if value is None:
            try:
                self._tcomp = self.rrock._comp+self.fluid._comp
            except Exception as e:
                logging.warning(f"Missing attribute when calculating total compressibility: {e}")
        else:
            self._tcomp = value/6894.76

    @property
    def hdiff(self):
        return self._hdiff
    
    @hdiff.setter
    def hdiff(self,value):
        self._hdiff = (self.rrock._xperm)/(self.rrock._poro*self.fluid._visc*self._tcomp)

    @property
    def pterm(self):
        return self._pterm/6894.76
    
    @pterm.setter
    def pterm(self,value):
        self._pterm = (self.well._cond)/(2*numpy.pi*self._xflow*self.fluid._mobil)

@dataclass(frozen=True)
class Boundary:

    #shape factor, {C_A} value
    factor: float

    # PSS is exact for higher values
    time_pss_accurate: float = None
    # PSS gives less than 1% error for higher values
    time_pss_error_prone: float = None
    # Use infinite system solution with less than 1 % Error for lesser values
    time_infinite: float = None

class Solution():

    def __init__(self,times,points):

        self.times  = times
        self.points = points

    @property
    def times(self):
        return self._times/(24*60*60)

    @times.setter
    def times(self,values:numpy.ndarray):
        self._times = numpy.ravel(values).reshape((1,-1))*(24*60*60)

    @property
    def points(self):
        return self._points/0.3048

    @points.setter
    def points(self,values:numpy.ndarray):
        self._points = numpy.ravel(values).reshape((-1,1))*0.3048

    @property
    def press(self):
        return self._press/6894.76

    @press.setter
    def press(self,values):
        self._press = values*6894.76