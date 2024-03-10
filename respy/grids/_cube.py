import sys

import numpy

if __name__ == "__main__":
	sys.path.append(r'C:\Users\3876yl\Documents\respy')

from respy.grids._cvol import ConVolume

class RecCube(numpy.ndarray):
	"""It is a numpy.ndarray interface for control volume."""

	def __new__(cls,*args,cvol:ConVolume=None,rows=None):
		"""
		edge 	: edge of the rectnagular cuboids, Nx3 numpy ndarray

		If edge is None, ConVolume instance should be provided:

		cvol    : ConVolume instance

		If rows is None, ConVolume's rows property is used:

		rows 	: Row indices to access from ConVolume
		"""

		if len(args)>1:
			raise Warning("RecCube does not accept more than one positional argument")
		
		if len(args)==1:
			cvol = ConVolume(args[0])
		
		if rows is None:
			rows = cvol.rows

		obj = numpy.asarray(rows,dtype=numpy.int_).view(cls)

		obj.cvol = cvol

		return obj

	def __array_finalize__(self,obj):

		if obj is None: return

		self.cvol = getattr(obj,'cvol',None)

	def __getattr__(self,key):
		
		return getattr(self.cvol,key)[self]

if __name__ == "__main__":

	import numpy as np

	edge = np.array([[1,2,4],[5,6,7],[7,8,9],[10,11,12]])

	cube1 = RecCube(edge)

	print(cube1.xarea)

	cvol = ConVolume(edge)

	cube2 = RecCube(cvol=cvol,rows=(1,2))

	print(cube2.xarea)

	# grid = RectCubGrid((750,1000,1250),(750,1000,1250),(20,))

	# print(grid.plat)

	# print(grid.size)

	# print(grid.size_test)

	# print(grid.size_test.ypos)

	# print(grid.xneg1._area)

	# print(area.xmin)