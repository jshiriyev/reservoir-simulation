import numpy

from ._base import GridBase
from ._grids import Grids

class GridDelta(GridBase):
	"""Interface to get a three-dimensional rectangular cuboid for cell-based flow simulations."""

	def __init__(self,xdelta:numpy.ndarray,ydelta:numpy.ndarray,zdelta:numpy.ndarray,dims:int=None):
		"""
		Initializes a 3D grid structure with spatial discretization.
		
		Parameters:

		xdelta	: length of grids in x-direction (feet), shape = (xnums,)
		ydelta	: width of grids in y-direction (feet), shape = (ynums,)
		zdelta	: height of grids in z-direction (feet), shape = (znums,)
		
		dims 	: (optional) flow dimension, determines the shape of neighborhood table.

		The object has edge properties to calculate the control
		volume implementation of flow calculations.
		"""
		super().__init__(xdelta,ydelta,zdelta)

		self.dims = dims

		self.length = None # Placeholder for reservoir domain total length calculations.
		self.width = None  # Placeholder for reservoir domain total width calculations.
		self.height = None # Placeholder for reservoir domain total height calculations.

	@property
	def nums(self):
		"""Returns tuple of (xnums,ynums,znums)."""
		return (self.xnums,self.ynums,self.znums)

	@property
	def xnums(self):
		"""Returns the number of grids in x-direction."""
		return self.xdelta.size

	@property
	def ynums(self):
		"""Returns the number of grids in y-direction."""
		return self.ydelta.size

	@property
	def znums(self):
		"""Returns the number of grids in z-direction."""
		return self.zdelta.size

	@property
	def dims(self):
		"""Returns the flow dimensions."""
		return self._dims

	@dims.setter
	def dims(self,value):
		"""Calculates the flow dimensions."""
		if value is not None:
			self._dims = value
		elif self.nums[2]>1:
			self._dims = 3
		elif self.nums[1]>1:
			self._dims = 2
		else:
			self._dims = 1

	@property
	def length(self):
		"""Returns the size of reservoir domain in x-direction."""
		return self._length/0.3048

	@length.setter
	def length(self,value):
		"""Calculates the size of reservoir domain in x-direction."""
		self._length = self._xdelta.sum()

	@property
	def width(self):
		"""Returns the size of reservoir domain in y-direction."""
		return self._width/0.3048

	@width.setter
	def width(self,value):
		"""Calculates the size of reservoir domain in y-direction."""
		self._width = self._ydelta.sum()

	@property
	def height(self):
		"""Returns the size of reservoir domain in z-direction."""
		return self._height/0.3048

	@height.setter
	def height(self,value):
		"""Calculates the size of reservoir domain in z-direction."""
		self._height = self._zdelta.sum()

	def center(self):
		"""Calculates the midpoints of grids."""
		self.xcenter = None
		self.ycenter = None
		self.zcenter = None
		return self

	@property
	def xcenter(self):
		"""Returns the x-center of the grids, shape = (xnums,)."""
		return self._xcenter/0.3048

	@xcenter.setter
	def xcenter(self,value):
		"""Calculates the x-center of the grids, shape = (xnums,)."""
		self._xcenter = self.midpoints(self._xdelta)

	@property
	def ycenter(self):
		"""Returns the y-center of the grids, shape = (ynums,)."""
		return self._ycenter/0.3048

	@ycenter.setter
	def ycenter(self,value):
		"""Calculates the y-center of the grids, shape = (ynums,)."""
		self._ycenter = self.midpoints(self._ydelta)

	@property
	def zcenter(self):
		"""Returns the z-center of the grids, shape = (znums,)."""
		return self._zcenter/0.3048

	@zcenter.setter
	def zcenter(self,value):
		"""Calculates the z-center of the grids, shape = (znums,)."""
		self._zcenter = self.midpoints(self._zdelta)

	@staticmethod
	def midpoints(array:numpy.ndarray):
		"""Calculates the midpoints of the given array."""
		pos = numpy.cumsum(array)
		neg = numpy.roll(pos,1)
		neg[0] = 0
		return (neg+pos)/2

	@property
	def index(self):
		"""Returns the indices of all grids."""
		return numpy.arange(numpy.prod(self.nums),dtype=numpy.int_)

	@property
	def grids(self):
		"""Returns Grids instance necessary for flow calculations."""
		xdelta = numpy.tile(self.xdelta,self.ynums*self.znums)
		ydelta = numpy.tile(numpy.repeat(self.ydelta,self.xnums),self.znums)
		zdelta = numpy.repeat(self.zdelta,self.xnums*self.ynums)

		return Grids(xdelta,ydelta,zdelta,self.table)
	
	@property
	def table(self):
		"""Returns the table of grids that stores neighborhood indices."""
		map_ = numpy.tile(self.index,(self.dims*2,1)).T

		map_[self.index.reshape(-1,self.xnums)[:,1:].ravel(),0] -= 1
		map_[self.index.reshape(-1,self.xnums)[:,:-1].ravel(),1] += 1

		if self.dims>1:
			map_[self.index.reshape(self.znums,-1)[:,self.xnums:],2] -= self.xnums
			map_[self.index.reshape(self.znums,-1)[:,:-self.xnums],3] += self.xnums

		if self.dims>2:
			map_[self.index.reshape(self.znums,-1)[1:,:],4] -= self.xnums*self.ynums
			map_[self.index.reshape(self.znums,-1)[:-1,:],5] += self.xnums*self.ynums

		return map_

if __name__ == "__main__":

	grid = GridDelta((750,1000,1250),(750,1000,1250),(20,),)

	print(grid.xcenter)

	print(grid.xpos)
	print(grid.xposb)

	print(grid.zpos)
	print(grid.zposb)