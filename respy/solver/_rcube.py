import sys

if __name__ == "__main__":
    # sys.path.append(r'C:\Users\javid.shiriyev\Documents\respy')
    sys.path.append(r'C:\Users\3876yl\Documents\respy')

import numpy

from numpy import ndarray

class RecCube():
	"""Rectangular Cuboid class"""

	def __init__(self,edge:ndarray,plat:ndarray,rows:ndarray=None,prop:set=None,**kwargs):
		"""
		Rectangular Cuboid class holding cell properties:
		
		edge 	: edge size of each grids, (nums,3)
		
		plat 	: neighbour indices, (nums,2*dims) where
				  dims is the flow dimensions, it can be 1, 2, or 3
		
		rows 	: indices of cells, (nums,)

		prop 	: a set containing names of the cell properties

		**kwargs contains any key-value pair for cell properties.

		"""

		object.__setattr__(self,"edge",edge)
		object.__setattr__(self,"plat",plat)

		if rows is None:
			rows = numpy.arange(
				self.edge[:,0].size,dtype=numpy.int_
				)

		object.__setattr__(self,"rows",rows)

		if prop is None:
			prop = {"edge","plat","rows","prop"}

		object.__setattr__(self,"prop",prop)

		for key,value in kwargs.items():
			self.__setattr__(key,value)

	def __setattr__(self,key,value):

		object.__setattr__(self,key,value)
		
		self.prop.add(key)

	def __getitem__(self,key):

		kwargs = {}

		for item in self.prop:

			prop = getattr(self,item)

			if not hasattr(prop,"shape"):
				kwargs[item] = prop
			elif prop.shape[0]==self.rows.size:
				kwargs[item] = prop[key,]
			else:
				kwargs[item] = prop

		return RecCube(**kwargs)

	@property
	def shape(self):
		return (self.rows.size,)

	@property
	def dims(self):
		return self.plat.shape[1]//2
	
	@property
	def xedge(self):
		return self.edge[:,0]

	@property
	def yedge(self):
		return self.edge[:,1]

	@property
	def zedge(self):
		return self.edge[:,2]

	@property
	def xarea(self):
		return self.yedge*self.zedge

	@property
	def yarea(self):
		return self.zedge*self.xedge

	@property
	def zarea(self):
		return self.xedge*self.yedge

	@property
	def volume(self):
		return self.xedge*self.yedge*self.zedge

	@property
	def xmin(self):
		"""x-minimum boundary boolean"""
		return self[self.rows==self.plat[:,0]]

	@property
	def xpos(self):
		"""x-positive neighbors boolean"""
		return self[self.rows!=self.plat[:,0]]

	@property
	def xneg(self):
		"""x-negative neighbors boolean"""
		return self[self.rows!=self.plat[:,1]]

	@property
	def xmax(self):
		"""x-maximum boundary boolean"""
		return self[self.rows==self.plat[:,1]]

	@property
	def ymin(self):
		"""y-minimum boundary boolean"""
		if self.dims>1:
			return self[self.rows==self.plat[:,2]]
		return self

	@property
	def ypos(self):
		"""y-positive neighbors boolean"""
		if self.dims>1:
			return self[self.rows!=self.plat[:,2]]
		return self[numpy.full(self.rows.shape,False)]

	@property
	def yneg(self):
		"""y-negative neighbors boolean"""
		if self.dims>1:
			return self[self.rows!=self.plat[:,3]]
		return self[numpy.full(self.rows.shape,False)]

	@property
	def ymax(self):
		"""y-maximum boundary boolean"""
		if self.dims>1:
			return self[self.rows==self.plat[:,3]]
		return self

	@property
	def zmin(self):
		"""z-minimum boundary boolean"""
		if self.dims>2:
			return self[self.rows==self.plat[:,4]]
		return self

	@property
	def zpos(self):
		"""z-positive neighbors boolean"""
		if self.dims>2:
			return self[self.rows!=self.plat[:,4]]
		return self[numpy.full(self.rows.shape,False)]

	@property
	def zneg(self):
		"""z-negative neighbors boolean"""
		if self.dims>2:
			return self[self.rows!=self.plat[:,5]]
		return self[numpy.full(self.rows.shape,False)]

	@property
	def zmax(self):
		"""z-maximum boundary boolean"""
		if self.dims>2:
			return self[self.rows==self.plat[:,5]]
		return self

if __name__ == "__main__":

	cube = RecCube(
		edge=numpy.array([[1,3,3],[2,3,3],[5,3,3],[7,3,3]]),
		plat=numpy.array([[0,1],[0,2],[1,3],[2,3]]),
		perm=numpy.array((400,500,600,700)),
		comp=numpy.array((1,)),
		func=lambda x: x**2)

	print(1,cube)

	print(2,cube.xarea)

	print(3,cube.xmin.rows)
	print(4,cube.xmax.rows)
	print(5,cube.ymin.rows)
	print(6,cube.ymax.rows)
	print(7,cube.zmin)
	print(8,cube.zmax)

	print(9,cube.xpos.rows)
	print(10,cube.xneg.rows)
	print(11,cube.ypos.rows)
	print(12,cube.yneg.rows)
	print(13,cube.zpos.rows)
	print(14,cube.zneg.rows)

	print(15,cube.prop)
	print(15,cube.rows)

	sub = cube[:2]

	print(16,sub.xmin.rows)
	print(17,sub.xmax.rows)
	print(18,sub.ymin.rows)
	print(19,sub.ymax.rows)
	print(20,sub.zmin.rows)
	print(21,sub.zmax.rows)

	print(22,sub.xpos.rows)
	print(23,sub.xneg.rows)
	print(24,sub.ypos.rows)
	print(25,sub.yneg.rows)
	print(26,sub.zpos.rows)
	print(27,sub.zneg.perm)

	print(28,sub.xarea)

	print(29,sub.perm)

	print(30,sub.comp)

	print(31,sub.func)