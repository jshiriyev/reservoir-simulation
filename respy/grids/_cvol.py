import numpy

class ConVolume():
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

if __name__ == "__main__":

	vol = ConVolume(
		numpy.array([[1,2,3],[1,2,3],[1,2,3],[1,2,3]]),
		perm=numpy.array((400,500,600,700)))

	print(type(vol[:2]))

	vol1 = vol[:2]

	print(vol1.xarea)

	print(vol1.perm)