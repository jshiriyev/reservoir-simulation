import numpy

class GridBase():

	def __init__(self,xdelta,ydelta,zdelta):
		"""Initialize grid deltas in feet. Internally stored in meters."""
		self.xdelta = xdelta # ft
		self.ydelta = ydelta # ft
		self.zdelta = zdelta # ft
		
	@property
	def xdelta(self):
		return self._xdelta/0.3048

	@xdelta.setter
	def xdelta(self,value):
		self._xdelta = numpy.asarray(value).flatten().astype(numpy.float64)*0.3048

	@property
	def ydelta(self):
		return self._ydelta/0.3048

	@ydelta.setter
	def ydelta(self,value):
		self._ydelta = numpy.asarray(value).flatten().astype(numpy.float64)*0.3048

	@property
	def zdelta(self):
		return self._zdelta/0.3048

	@zdelta.setter
	def zdelta(self,value):
		self._zdelta = numpy.asarray(value).flatten().astype(numpy.float64)*0.3048