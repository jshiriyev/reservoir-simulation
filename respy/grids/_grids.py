import numpy

from ._base import GridBase

class Grids(GridBase):

    def __init__(self,plat,xdelta,ydelta,zdelta):        
        super().__init__(xdelta,ydelta,zdelta)
        self.plat = plat
        self.volume = None

    @property
    def nums(self):
        """Returns total number of grids."""
        return self.plat.shape[0]

    @property
    def dims(self):
        return int(self.plat.shape[1]/2)

    @property
    def index(self):
        return numpy.arange(self.nums)

    @property
    def volume(self):
        return self._volume*35.3147

    @volume.setter
    def volume(self,value):
        self._volume = numpy.prod((self._xdelta,self._ydelta,self._zdelta),axis=0)