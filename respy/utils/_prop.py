import numpy

class Prop():
	"""
	This class helps to treat constant values as
	function of *args and **kwargs while returning
	their constant value only.
	"""

	def __init__(self,data=None):

		self.data = data

	@property
	def isconst(self):
		"""Checks whether data is callable or not."""

		if self.data is None:
			return
		
		if callable(self.data):
			return False

		return True

	def __call__(self,*args,**kwargs):

		if callable(self.data):
			return self.data(*args,**kwargs)
		
		return self.reshape(*args,**kwargs)

	def reshape(self,*args,**kwargs):

		data = numpy.asarray(self.data)

		if data.shape==args[0].shape:
			return data

		data = data.flatten()

		if data.size==1:
			return data*numpy.ones_like(args[0])

		return data.reshape((-1,1))

if __name__ == "__main__":

	t = lambda v: v**2

	x = Prop(1)
	y = Prop(t)
	z = Prop(None)

	print(x(numpy.array([4])))
	print(y(numpy.array([4])))
	print(z(numpy.array([5])))

	print(x.isconst)
	print(y.isconst)
	print(z.isconst)