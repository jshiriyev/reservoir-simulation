import numpy

class RecCube():
	"""Reservoir Cuboid class"""

	def __init__(self,edge:numpy.ndarray,plat:numpy.ndarray,rows:numpy.ndarray=None,**kwargs):
		"""
		Reservoir Cuboid class holding cell properties:
		
		edge 	: edge size of each grid, (N,3) numpy.ndarray
		plat 	: neighboring indices, (N,2*dims) numpy.ndarray where dims is 1, 2, or 3
		rows 	: indices of cells, (N,) flat numpy.ndarray

		The value of kwargs should be flat array of numpy.ndarray of size N or 1.

		"""

		object.__setattr__(self,"edge",edge)
		object.__setattr__(self,"plat",plat)

		if rows is None:
			rows = numpy.arange(
				self.edge[:,0].size,dtype=numpy.int_
				)

		object.__setattr__(self,"rows",rows)
		object.__setattr__(self,"prop",set())

		for key,value in kwargs.items():
			self.__setattr__(key,value)

	def __setattr__(self,key,value):
		"""The value should be flat array of numpy.ndarray of size N or 1"""

		object.__setattr__(self,key,value)
		
		self.prop.add(key)

	def __getitem__(self,key):

		kwargs = {}

		for item in self.prop:

			prop = getattr(self,item)

			if callable(prop):
				kwargs[item] = prop
			elif prop.size>1 or self.rows.size<=1:
				kwargs[item] = prop[key]
			else:
				kwargs[item] = prop

		return RecCube(self.edge[key],self.plat[key],self.rows[key],**kwargs)

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

	vol = RecCube(
		numpy.array([[1,3,3],[1,2,3],[1,7,3],[1,5,3]]),
		numpy.array([[0,1],[0,2],[1,3],[2,3]]),
		perm=numpy.array((400,500,600,700)),
		comp=numpy.array((1,)),
		name=lambda x: x**2)

	print(vol)

	print(vol.xarea)

	print(vol.xmin.rows)
	print(vol.xmax.rows)
	print(vol.ymin.rows)
	print(vol.ymax.rows)
	print(vol.zmin.rows)
	print(vol.zmax.rows)

	print(vol.xpos.rows)
	print(vol.xneg.rows)
	print(vol.ypos.rows)
	print(vol.yneg.rows)
	print(vol.zpos.rows)
	print(vol.zneg.rows)

	print(vol[:2])

	vol1 = vol[:2]

	print(vol1.xmin.rows)
	print(vol1.xmax.rows)
	print(vol1.ymin.rows)
	print(vol1.ymax.rows)
	print(vol1.zmin.rows)
	print(vol1.zmax.rows)

	print(1,vol1.xpos.rows)
	print(vol1.xneg.rows)
	print(vol1.ypos.rows)
	print(vol1.yneg.rows)
	print(vol1.zpos.rows)
	print(vol1.zneg.perm)

	print(vol1.xarea)

	print(vol1.perm)

	print(vol1.comp)

	print(vol1.name)