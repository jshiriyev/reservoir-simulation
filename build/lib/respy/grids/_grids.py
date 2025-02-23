import numpy

from ._base import GridBase

class Grids(GridBase):
    """Represents a structured reservoir grid with additional table attribute."""

    def __init__(self,xdelta:numpy.ndarray,ydelta:numpy.ndarray,zdelta:numpy.ndarray,table:numpy.ndarray):
        """
        Initialize the grid with spatial discretization and associated table data.

        Parameters:
        xdelta (np.ndarray): Grid cell sizes in the x-direction (feet).
        ydelta (np.ndarray): Grid cell sizes in the y-direction (feet).
        zdelta (np.ndarray): Grid cell sizes in the z-direction (feet).

        table (np.ndarray): Data table mapping the neighbors of each grid, integers.
        """
        super().__init__(xdelta,ydelta,zdelta)

        self.volume = None # Placeholder for volume calculations

        self.table  = table

    @property
    def delta(self):
        """Returns the cell sizes in x, y, and z direction, shape = (nums,3)."""
        return numpy.column_stack(
            (self.xdelta,self.ydelta,self.zdelta)
            )

    @property
    def volume(self):
        """Returns the volume of grids in field units."""
        return self._volume*35.3147

    @volume.setter
    def volume(self,value):
        """Calculates the volume of each cell."""
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
        """Returns x-minimum boundary indices."""
        return self.index[self.xminbool]

    @property
    def xpos(self):
        """Returns x-positive neighbor indices."""
        return self.index[self.xposbool]

    @property
    def xneg(self):
        """Returns x-negative neighbor indices."""
        return self.index[self.xnegbool]

    @property
    def xmax(self):
        """Returns x-maximum boundary indices."""
        return self.index[self.xmaxbool]

    @property
    def ymin(self):
        """Returns y-minimum boundary indices."""
        return self.index[self.yminbool]

    @property
    def ypos(self):
        """Returns y-positive neighbor indices."""
        return self.index[self.yposbool]

    @property
    def yneg(self):
        """Returns y-negative neighbor indices."""
        return self.index[self.ynegbool]

    @property
    def ymax(self):
        """Returns y-maximum boundary indices."""
        return self.index[self.ymaxbool]

    @property
    def zmin(self):
        """Returns z-minimum boundary indices."""
        return self.index[self.zminbool]

    @property
    def zpos(self):
        """Returns z-positive neighbor indices."""
        return self.index[self.zposbool]

    @property
    def zneg(self):
        """Returns z-negative neighbor indices."""
        return self.index[self.znegbool]

    @property
    def zmax(self):
        """Returns z-maximum boundary indices."""
        return self.index[self.zmaxbool]
    
    @property
    def xminbool(self):
        """Returns x-minimum boundary boolean."""
        return self.index==self.table[:,0]

    @property
    def xposbool(self):
        """Returns x-positive neighbor boolean."""
        return self.index!=self.table[:,0]

    @property
    def xnegbool(self):
        """Returns x-negative neighbor boolean."""
        return self.index!=self.table[:,1]

    @property
    def xmaxbool(self):
        """Returns x-maximum boundary boolean."""
        return self.index==self.table[:,1]

    @property
    def yminbool(self):
        """Returns y-minimum boundary boolean."""
        return self.index==self.table[:,2] if self.dims>1 else numpy.full(self.nums,True)

    @property
    def yposbool(self):
        """Returns y-positive neighbor boolean."""
        return self.index!=self.table[:,2] if self.dims>1 else numpy.full(self.nums,False)

    @property
    def ynegbool(self):
        """Returns y-negative neighbor boolean."""
        return self.index!=self.table[:,3] if self.dims>1 else numpy.full(self.nums,False)

    @property
    def ymaxbool(self):
        """Returns y-maximum boundary boolean."""
        return self.index==self.table[:,3] if self.dims>1 else numpy.full(self.nums,True)

    @property
    def zminbool(self):
        """Returns z-minimum boundary boolean."""
        return self.index==self.table[:,4] if self.dims>2 else numpy.full(self.nums,True)

    @property
    def zposbool(self):
        """Returns z-positive neighbor boolean."""
        return self.index!=self.table[:,4] if self.dims>2 else numpy.full(self.nums,False)

    @property
    def znegbool(self):
        """Returns z-negative neighbor boolean."""
        return self.index!=self.table[:,5] if self.dims>2 else numpy.full(self.nums,False)

    @property
    def zmaxbool(self):
        """Returns z-maximum boundary boolean."""
        return self.index==self.table[:,5] if self.dims>2 else numpy.full(self.nums,True)