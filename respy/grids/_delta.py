import sys

import numpy

if __name__ == "__main__":
	sys.path.append(r'C:\Users\3876yl\Documents\respy')

from respy.grids._cube import RecCube
from respy.grids._cvol import ConVolume

class GridDelta():

	def __init__(self,xdelta:tuple,ydelta:tuple,zdelta:tuple,dims=None):
		"""Three-dimensional rectangular cuboid initialized with:

		xdelta	: length of grids in ft, shape = (Nlength,)
		ydelta	: width of grids in ft, shape = (Nwidth,)
		zdelta	: height of grids in ft, shape = (Nheight,)
		
		dims 	: flow dimension, determines the shape of plat.

		The object has properties to calculate the control volume implementation.
		"""

		self.xdelta = self.set_delta(xdelta)
		self.ydelta = self.set_delta(ydelta)
		self.zdelta = self.set_delta(zdelta)

		self.dims = self.set_dims(dims)

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
	def cvol(self):
		return ConVolume(self.edge)
	
	@property
	def cube(self):
		return RecCube(self.edge)

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

	@property
	def xmin(self):
		"""Properties of grids on x-minimum boundary"""
		cvol,xmin_rows = self.cvol,self.plat[:,0]
		return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows==xmin_rows])

	@property
	def xpos(self):
		"""Properties of grids that has x-positive neighbors"""
		cvol,xmin_rows = self.cvol,self.plat[:,0]
		return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows!=xmin_rows])

	@property
	def xneg(self):
		"""Properties of grids that has x-negative neighbors"""
		cvol,xmax_rows = self.cvol,self.plat[:,1]
		return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows!=xmax_rows])

	@property
	def xmax(self):
		"""Properties of grids on x-maximum boundary"""
		cvol,xmax_rows = self.cvol,self.plat[:,1]
		return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows==xmax_rows])

	@property
	def ymin(self):
		"""Properties of grids on y-minimum boundary"""
		if self.dims>1:
			cvol,ymin_rows = self.cvol,self.plat[:,2]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows==ymin_rows])

	@property
	def ypos(self):
		"""Properties of grids that has y-positive neighbors"""
		if self.dims>1:
			cvol,ymin_rows = self.cvol,self.plat[:,2]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows!=ymin_rows])

	@property
	def yneg(self):
		"""Properties of grids that has y-negative neighbors"""
		if self.dims>1:
			cvol,ymax_rows = self.cvol,self.plat[:,3]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows!=ymax_rows])

	@property
	def ymax(self):
		"""Properties of grids on y-maximum boundary"""
		if self.dims>1:
			cvol,ymax_rows = self.cvol,self.plat[:,3]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows==ymax_rows])

	@property
	def zmin(self):
		"""Properties of grids on z-minimum boundary"""
		if self.dims>2:
			cvol,zmin_rows = self.cvol,self.plat[:,4]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows==zmin_rows])

	@property
	def zpos(self):
		"""Properties of grids that has z-positive neighbors"""
		if self.dims>2:
			cvol,zmin_rows = self.cvol,self.plat[:,4]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows!=zmin_rows])

	@property
	def zneg(self):
		"""Properties of grids that has z-negative neighbors"""
		if self.dims>2:
			cvol,zmax_rows = self.cvol,self.plat[:,5]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows!=zmax_rows])

	@property
	def zmax(self):
		"""Properties of grids on z-maximum boundary"""
		if self.dims>2:
			cvol,zmax_rows = self.cvol,self.plat[:,5]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows==zmax_rows])

if __name__ == "__main__":

	grid = GridDelta((750,1000,125),(750,1000,1250),(20,))

	# print(rcube.ymin)

	# print(grid.cube.xarea)

	print(grid.xmin)
	print(grid.xmax)
	print(grid.ymin)
	print(grid.ymax)
	print(grid.zmin)
	print(grid.zmax)

	print(grid.xpos)
	print(grid.xneg)
	print(grid.ypos)
	print(grid.yneg)
	print(grid.zpos)
	print(grid.zneg)