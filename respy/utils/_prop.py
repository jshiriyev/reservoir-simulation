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
		
		return self.data

if __name__ == "__main__":

	t = lambda v: v**2

	x = Prop(1)
	y = Prop(t)
	z = Prop(None)

	print(x(4))
	print(y(4))
	print(z(5))

	print(x.isconst)
	print(y.isconst)
	print(z.isconst)