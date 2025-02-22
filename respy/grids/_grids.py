import numpy

from ._base import GridBase

class Grids(GridBase):

    def __init__(self,xdelta,ydelta,zdelta,table):        
        super().__init__(xdelta,ydelta,zdelta)

        self.volume = None

        self.table = table

    @property
    def delta(self):
        """Returns the size of edges in x, y, and z direction, shape = (nums,3)."""
        return numpy.column_stack(
            (self.xdelta,self.ydelta,self.zdelta)
            )

    @property
    def volume(self):
        """Returns the volume of grids in field units."""
        return self._volume*35.3147

    @volume.setter
    def volume(self,value):
        self._volume = numpy.prod((self._xdelta,self._ydelta,self._zdelta),axis=0)

    @property
    def nums(self):
        """Returns total number of grids."""
        return self.table.shape[0]

    @property
    def dims(self):
        """Returns the flow dimensions."""
        return int(self.table.shape[1]/2)

    @property
    def index(self):
        """Returns the indices of grids."""
        return numpy.arange(self.nums)

    @property
    def xmin(self):
        """x-minimum boundary indices"""
        return self.index[self.xminbool]

    @property
    def xpos(self):
        """x-positive neighbors indices"""
        return self.index[self.xposbool]

    @property
    def xneg(self):
        """x-negative neighbors indices"""
        return self.index[self.xnegbool]

    @property
    def xmax(self):
        """x-maximum boundary indices"""
        return self.index[self.xmaxbool]

    @property
    def ymin(self):
        """y-minimum boundary indices"""
        return self.index[self.yminbool]

    @property
    def ypos(self):
        """y-positive neighbors indices"""
        return self.index[self.yposbool]

    @property
    def yneg(self):
        """y-negative neighbors indices"""
        return self.index[self.ynegbool]

    @property
    def ymax(self):
        """y-maximum boundary indices"""
        return self.index[self.ymaxbool]

    @property
    def zmin(self):
        """z-minimum boundary indices"""
        return self.index[self.zminbool]

    @property
    def zpos(self):
        """z-positive neighbors indices"""
        return self.index[self.zposbool]

    @property
    def zneg(self):
        """z-negative neighbors indices"""
        return self.index[self.znegbool]

    @property
    def zmax(self):
        """z-maximum boundary indices"""
        return self.index[self.zmaxbool]
    
    @property
    def xminbool(self):
        """x-minimum boundary boolean"""
        return self.index==self.table[:,0]

    @property
    def xposbool(self):
        """x-positive neighbors boolean"""
        return self.index!=self.table[:,0]

    @property
    def xnegbool(self):
        """x-negative neighbors boolean"""
        return self.index!=self.table[:,1]

    @property
    def xmaxbool(self):
        """x-maximum boundary boolean"""
        return self.index==self.table[:,1]

    @property
    def yminbool(self):
        """y-minimum boundary boolean"""
        return self.index==self.table[:,2] if self.dims>1 else numpy.full(self.nums,True)

    @property
    def yposbool(self):
        """y-positive neighbors boolean"""
        return self.index!=self.table[:,2] if self.dims>1 else numpy.full(self.nums,False)

    @property
    def ynegbool(self):
        """y-negative neighbors boolean"""
        return self.index!=self.table[:,3] if self.dims>1 else numpy.full(self.nums,False)

    @property
    def ymaxbool(self):
        """y-maximum boundary boolean"""
        return self.index==self.table[:,3] if self.dims>1 else numpy.full(self.nums,True)

    @property
    def zminbool(self):
        """z-minimum boundary boolean"""
        return self.index==self.table[:,4] if self.dims>2 else numpy.full(self.nums,True)

    @property
    def zposbool(self):
        """z-positive neighbors boolean"""
        return self.index!=self.table[:,4] if self.dims>2 else numpy.full(self.nums,False)

    @property
    def znegbool(self):
        """z-negative neighbors boolean"""
        return self.index!=self.table[:,5] if self.dims>2 else numpy.full(self.nums,False)

    @property
    def zmaxbool(self):
        """z-maximum boundary boolean"""
        return self.index==self.table[:,5] if self.dims>2 else numpy.full(self.nums,True)