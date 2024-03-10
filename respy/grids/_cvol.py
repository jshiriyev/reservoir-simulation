import numpy

class ConVolume():
	"""Collection of the properties to calculate for the reservoir simulation."""

	def __init__(self,edge:numpy.ndarray,rows=None):

		self.edge = edge

		if rows is None:
			rows = numpy.arange(
				self.edge[:,0].size,dtype=numpy.int_
				)

		self.rows = rows

	def __getitem__(self,key):

		return ConVolume(self.edge[key],self.rows[key])

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

	vol = ConVolume(numpy.array([[1,2,3],[1,2,3],[1,2,3],[1,2,3]]))

	print(type(vol[:2]))

	vol1 = vol[:2]

	print(vol1.xarea)