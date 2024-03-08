from inspect import isfunction

class Prop():

	def __init__(self,data=None,conv=1):

		self.conv = conv
		self.data = data

	def set(self,*args,**kwargs):

		if isfunction(self.data):
			return self.data(*args,**kwargs)*self.conv

		return self.data*self.conv

	def get(self,*args,**kwargs):

		if isfunction(self.data):
			return self.data(*args,**kwargs)/self.conv

		return self.data/self.conv