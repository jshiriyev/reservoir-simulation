import sys

import numpy

class GridDelta():

	def __init__(self,xdelta:tuple,ydelta:tuple,zdelta:tuple,dims=None):
		"""
		Three-dimensional rectangular cuboid initialized with:

		xdelta	: length of grids in ft, shape = (xnums,)
		ydelta	: width of grids in ft, shape = (ynums,)
		zdelta	: height of grids in ft, shape = (znums,)
		
		dims 	: flow dimension, determines the shape of plat.

		The object has edge properties to calculate the control
		volume implementation of flow calculations.
		"""

		self.xdelta = self.get_delta(xdelta)
		self.ydelta = self.get_delta(ydelta)
		self.zdelta = self.get_delta(zdelta)

		self.__dims = dims

	"""Properties essential for the visualization of results:"""

	@property
	def nums(self):
		"""Returns tuple of (xnums,ynums,znums)."""
		return (self.xnums,self.ynums,self.znums)

	@property
	def xnums(self):
		return self.xdelta.size

	@property
	def ynums(self):
		return self.ydelta.size

	@property
	def znums(self):
		return self.zdelta.size

	@property
	def length(self):
		"""Returns the size of reservoir in x-direction."""
		return self.xdelta.sum()

	@property
	def width(self):
		"""Returns the size of reservoir in y-direction."""
		return self.ydelta.sum()

	@property
	def height(self):
		"""Returns the size of reservoir in z-direction."""
		return self.zdelta.sum()

	@property
	def xcenter(self):
		"""Returns the x-center of the grids, shape = (xnums,)"""
		side2 = numpy.cumsum(self.xdelta)
		side1 = numpy.roll(side2,1)

		side1[0] = 0

		return (side1+side2)/2

	@property
	def ycenter(self):
		"""Returns the y-center of the grids, shape = (ynums,)"""
		side2 = numpy.cumsum(self.ydelta)
		side1 = numpy.roll(side2,1)

		side1[0] = 0

		return (side1+side2)/2

	@property
	def zcenter(self):
		"""Returns the z-center of the grids, shape = (znums,)"""
		side2 = numpy.cumsum(self.zdelta)
		side1 = numpy.roll(side2,1)

		side1[0] = 0

		return (side1+side2)/2

	"""Properties essential for the flow simulation:"""

	def __call__(self):

		edge = self.edge*0.3048
		plat = self.plat
		rows = self.rows

		return (edge,plat,rows)

	@property
	def nums(self):
		"""Returns total number of grids."""
		return numpy.prod(self.nums3).item()

	@property
	def dims(self):
		"""Returns the flow dimensions."""
		return self.get_dims(dims=self.__dims,nums3=self.nums3)
	
	@property
	def edge(self):
		"""Returns the size of edges in x, y, and z direction, shape = (nums,3)."""

		x = self.xedge.reshape((-1,1))
		y = self.yedge.reshape((-1,1))
		z = self.zedge.reshape((-1,1))

		return numpy.concatenate((x,y,z),axis=1)

	@property
	def xedge(self):
		"""Returns the size of edges in x-direction, shape = (nums,)."""
		return numpy.tile(self.xdelta,self.ynums*self.znums)
	
	@property
	def yedge(self):
		"""Returns the size of edges in y-direction, shape = (nums,)."""
		return numpy.tile(numpy.repeat(self.ydelta,self.xnums),self.znums)
	
	@property
	def zedge(self):
		"""Returns the size of edges in z-direction, shape = (nums,)."""
		return numpy.repeat(self.zdelta,self.xnums*self.ynums)

	@property
	def rows(self):
		return numpy.arange(self.nums,dtype=numpy.int_)
	
	@property
	def plat(self):
		"""Plat of grids that locates neighbor index information."""

		dims = self.dims

		Nx,Ny,Nz = self.xnums,self.ynums,self.znums

		index = numpy.arange(self.nums)

		plat = numpy.tile(index,(dims*2,1)).T

		plat[index.reshape(-1,Nx)[:,1:].ravel(),0] -= 1
		plat[index.reshape(-1,Nx)[:,:-1].ravel(),1] += 1

		if dims>1:
			plat[index.reshape(Nz,-1)[:,Nx:],2] -= Nx
			plat[index.reshape(Nz,-1)[:,:-Nx],3] += Nx

		if dims>2:
			plat[index.reshape(Nz,-1)[1:,:],4] -= Nx*Ny
			plat[index.reshape(Nz,-1)[:-1,:],5] += Nx*Ny

		return plat

	"""Essential indices for grids:"""

	@property
	def xmin(self):
		"""x-minimum boundary indices"""
		return self.rows[self.xminb]

	@property
	def xpos(self):
		"""x-positive neighbors indices"""
		return self.rows[self.xposb]

	@property
	def xneg(self):
		"""x-negative neighbors indices"""
		return self.rows[self.xnegb]

	@property
	def xmax(self):
		"""x-maximum boundary indices"""
		return self.rows[self.xmaxb]

	@property
	def ymin(self):
		"""y-minimum boundary indices"""
		return self.rows[self.yminb]

	@property
	def ypos(self):
		"""y-positive neighbors indices"""
		return self.rows[self.yposb]

	@property
	def yneg(self):
		"""y-negative neighbors indices"""
		return self.rows[self.ynegb]

	@property
	def ymax(self):
		"""y-maximum boundary indices"""
		return self.rows[self.ymaxb]

	@property
	def zmin(self):
		"""z-minimum boundary indices"""
		return self.rows[self.zminb]

	@property
	def zpos(self):
		"""z-positive neighbors indices"""
		return self.rows[self.zposb]

	@property
	def zneg(self):
		"""z-negative neighbors indices"""
		return self.rows[self.znegb]

	@property
	def zmax(self):
		"""z-maximum boundary indices"""
		return self.rows[self.zmaxb]

	"""Essential booleans for grids:"""
	
	@property
	def xminb(self):
		"""x-minimum boundary boolean"""
		return self.rows==self.plat[:,0]

	@property
	def xposb(self):
		"""x-positive neighbors boolean"""
		return self.rows!=self.plat[:,0]

	@property
	def xnegb(self):
		"""x-negative neighbors boolean"""
		return self.rows!=self.plat[:,1]

	@property
	def xmaxb(self):
		"""x-maximum boundary boolean"""
		return self.rows==self.plat[:,1]

	@property
	def yminb(self):
		"""y-minimum boundary boolean"""
		if self.dims>1:
			return self.rows==self.plat[:,2]
		return numpy.full(self.nums,True)

	@property
	def yposb(self):
		"""y-positive neighbors boolean"""
		if self.dims>1:
			return self.rows!=self.plat[:,2]
		return numpy.full(self.nums,False)

	@property
	def ynegb(self):
		"""y-negative neighbors boolean"""
		if self.dims>1:
			return self.rows!=self.plat[:,3]
		return numpy.full(self.nums,False)

	@property
	def ymaxb(self):
		"""y-maximum boundary boolean"""
		if self.dims>1:
			return self.rows==self.plat[:,3]
		return numpy.full(self.nums,True)

	@property
	def zminb(self):
		"""z-minimum boundary boolean"""
		if self.dims>2:
			return self.rows==self.plat[:,4]
		return numpy.full(self.nums,True)

	@property
	def zposb(self):
		"""z-positive neighbors boolean"""
		if self.dims>2:
			return self.rows!=self.plat[:,4]
		return numpy.full(self.nums,False)

	@property
	def znegb(self):
		"""z-negative neighbors boolean"""
		if self.dims>2:
			return self.rows!=self.plat[:,5]
		return numpy.full(self.nums,False)

	@property
	def zmaxb(self):
		"""z-maximum boundary boolean"""
		if self.dims>2:
			return self.rows==self.plat[:,5]
		return numpy.full(self.nums,True)

	"""Static methods:"""

	@staticmethod
	def get_delta(delta):
		"""Returns flat numpy.ndarray of float dtype."""
		return numpy.asarray(delta).flatten().astype(numpy.float_)

	@staticmethod
	def get_dims(*,dims=None,nums3=None):
		"""Returns flow dimensions that determines shape of plat."""

		if dims is not None:
			return dims

		if nums3[2]>1:
			return 3

		if nums3[1]>1:
			return 2

		return 1

if __name__ == "__main__":

	grid = GridDelta((750,1000,1250),(750,1000,1250),(20,),)

	print(grid.xcenter)

	print(grid.xpos)
	print(grid.xposb)

	print(grid.zpos)
	print(grid.zposb)