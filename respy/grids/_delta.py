import sys

import numpy

if __name__ == "__main__":
	sys.path.append(r'C:\Users\3876yl\Documents\respy')

from respy.grids._cube import RecCube

class GridDelta(RecCube):

	def __init__(self,xdelta:tuple,ydelta:tuple,zdelta:tuple,depth=None,dims=None,**kwargs):
		"""Three-dimensional rectangular cuboid initialized with:

		xdelta	: length of grids in ft, shape = (Nlength,)
		ydelta	: width of grids in ft, shape = (Nwidth,)
		zdelta	: height of grids in ft, shape = (Nheight,)
		
		dims 	: flow dimension, determines the shape of plat.

		Any kargs item should follow the size regulation.

		The object has properties to calculate the control volume implementation.
		"""

		object.__setattr__(self,"xdelta",self.set_delta(xdelta))
		object.__setattr__(self,"ydelta",self.set_delta(ydelta))
		object.__setattr__(self,"zdelta",self.set_delta(zdelta))

		object.__setattr__(self,"dims",self.set_dims(dims))

		rows = numpy.arange(self.numtot,dtype=numpy.int_)

		object.__setattr__(self,"rows",rows)
		object.__setattr__(self,"prop",set())

		for key,value in kwargs.items():
			self.__setattr__(key,value)

	def set_delta(self,delta):

		return numpy.asarray(delta).flatten().astype(numpy.float_)

	def set_dims(self,dims):

		if dims is not None:
			return dims

		if self.znums>1:
			return 3

		if self.ynums>1:
			return 2

		return 1

	@property
	def nums(self):
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
	def numtot(self):
		return numpy.prod(self.nums).item()

	@property
	def length(self):
		return self.xdelta.sum()

	@property
	def width(self):
		return self.ydelta.sum()

	@property
	def height(self):
		return self.zdelta.sum()
	
	@property
	def edge(self):
		return numpy.concatenate((self.xedge,self.yedge,self.zedge),axis=1)

	@property
	def xedge(self):
		return numpy.tile(self.xdelta,self.ynums*self.znums).reshape((-1,1))
	
	@property
	def yedge(self):
		return numpy.tile(
			numpy.repeat(self.ydelta,self.xnums),self.znums).reshape((-1,1))
	
	@property
	def zedge(self):
		return numpy.repeat(self.zdelta,self.xnums*self.ynums).reshape((-1,1))

	@property
	def xcenter(self):

		side2 = numpy.cumsum(self.xdelta)
		side1 = numpy.roll(side2,1)

		side1[0] = 0

		return (side1+side2)/2

	@property
	def ycenter(self):

		side2 = numpy.cumsum(self.ydelta)
		side1 = numpy.roll(side2,1)

		side1[0] = 0

		return (side1+side2)/2

	@property
	def zcenter(self):

		side2 = numpy.cumsum(self.zdelta)
		side1 = numpy.roll(side2,1)

		side1[0] = 0

		return (side1+side2)/2

	@property
	def plat(self):

		dims = self.dims

		Nx,Ny,Nz,index = *self.nums,numpy.arange(self.numtot)

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

if __name__ == "__main__":

	grid = GridDelta((750,1000,1250),(750,1000,1250),(20,),
		perm=numpy.array((1,2,3,4,5,6,7,8,9)).reshape((-1,1)))

	# print(rcube.ymin)

	# print(grid.cube.xarea)

	print(grid.xmin.rows)
	print(grid.xmin.perm)

	print(grid.xcenter)

	# print(grid.xmin)
	# print(grid.xmax)
	# print(grid.ymin)
	# print(grid.ymax)
	# print(grid.zmin)
	# print(grid.zmax)

	# print(grid.xpos)
	# print(grid.xneg)
	# print(grid.ypos)
	# print(grid.yneg)
	# print(grid.zpos)
	# print(grid.zneg)