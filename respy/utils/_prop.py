class Prop():

	def __init__(self,data=None,conv=1):

		self.data = data
		self.conv = conv

	@property
	def const(self):

		if self.data is None:
			return
		
		if callable(self.data):
			return False

		return True

	def __call__(self,*args,**kwargs):

		try:
			multp = kwargs.pop("multp")
		except KeyError:
			multp = None

		if callable(self.data):
			data = self.data(*args,**kwargs)
		else:
			data = self.data

		if multp is None:
			return self.data
		
		if multp is False:
			return self.data/self.conv

		return self.data*self.conv

if __name__ == "__main__":

	x = Prop(1)

	# print(x.__value)

	# print(x.get())
	print(x(23))

	# for d in dir(Prop(5,2)):
	# 	print(d)