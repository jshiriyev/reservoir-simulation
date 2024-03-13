import numpy

class RecCube():
	"""Collection of the properties to calculate for the reservoir simulation."""

	def __init__(self,edge:numpy.ndarray,rows=None,**kwargs):
		"""Any kwargs values should follow the size regulation."""

		self.edge = edge

		if rows is None:
			rows = numpy.arange(
				self.edge[:,0].size,dtype=numpy.int_
				)

		self.rows = rows

		self.prop = set()

		self.set_prop(**kwargs)

	def __getitem__(self,key):

		kwargs = {}

		for item in self.prop:
			kwargs[item] = getattr(self,item)[key]

		return ConVolume(self.edge[key],self.rows[key],**kwargs)

	def set_prop(self,**kwargs):
		"""Any kwargs values should follow the size regulation."""

		for key,value in kwargs.items():
			setattr(self,key,value)
			self.prop.add(key)

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
		"""Properties of grids on x-minimum boundary"""
		cvol,xmin_rows = self.cube.cvol,self.plat[:,0]
		return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows==xmin_rows])

	@property
	def xpos(self):
		"""Properties of grids that has x-positive neighbors"""
		cvol,xmin_rows = self.cube.cvol,self.plat[:,0]
		return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows!=xmin_rows])

	@property
	def xneg(self):
		"""Properties of grids that has x-negative neighbors"""
		cvol,xmax_rows = self.cube.cvol,self.plat[:,1]
		return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows!=xmax_rows])

	@property
	def xmax(self):
		"""Properties of grids on x-maximum boundary"""
		cvol,xmax_rows = self.cube.cvol,self.plat[:,1]
		return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows==xmax_rows])

	@property
	def ymin(self):
		"""Properties of grids on y-minimum boundary"""
		if self.dims>1:
			cvol,ymin_rows = self.cube.cvol,self.plat[:,2]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows==ymin_rows])

	@property
	def ypos(self):
		"""Properties of grids that has y-positive neighbors"""
		if self.dims>1:
			cvol,ymin_rows = self.cube.cvol,self.plat[:,2]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows!=ymin_rows])

	@property
	def yneg(self):
		"""Properties of grids that has y-negative neighbors"""
		if self.dims>1:
			cvol,ymax_rows = self.cube.cvol,self.plat[:,3]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows!=ymax_rows])

	@property
	def ymax(self):
		"""Properties of grids on y-maximum boundary"""
		if self.dims>1:
			cvol,ymax_rows = self.cube.cvol,self.plat[:,3]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows==ymax_rows])

	@property
	def zmin(self):
		"""Properties of grids on z-minimum boundary"""
		if self.dims>2:
			cvol,zmin_rows = self.cube.cvol,self.plat[:,4]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows==zmin_rows])

	@property
	def zpos(self):
		"""Properties of grids that has z-positive neighbors"""
		if self.dims>2:
			cvol,zmin_rows = self.cube.cvol,self.plat[:,4]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows!=zmin_rows])

	@property
	def zneg(self):
		"""Properties of grids that has z-negative neighbors"""
		if self.dims>2:
			cvol,zmax_rows = self.cube.cvol,self.plat[:,5]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows!=zmax_rows])

	@property
	def zmax(self):
		"""Properties of grids on z-maximum boundary"""
		if self.dims>2:
			cvol,zmax_rows = self.cube.cvol,self.plat[:,5]
			return RecCube(cvol=cvol,rows=cvol.rows[cvol.rows==zmax_rows])

if __name__ == "__main__":

	vol = ConVolume(
		numpy.array([[1,2,3],[1,2,3],[1,2,3],[1,2,3]]),
		perm=numpy.array((400,500,600,700)))

	print(type(vol[:2]))

	vol1 = vol[:2]

	print(vol1.xarea)

	print(vol1.perm)