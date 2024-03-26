import sys

if __name__ == "__main__":
	# sys.path.append(r'C:\Users\javid.shiriyev\Documents\respy')
	sys.path.append(r'C:\Users\3876yl\Documents\respy')

import numpy

from respy.grids._delta import GridDelta

class GridRegular(GridDelta):

	def __init__(self,size:tuple,nums3:tuple,dims:int=None):
		"""
		Three-dimensional regular rectangular cuboid initialized with:

		size 	: reservoir dimensions (length, width, height)

		length 	: length of reservoir in ft, (x direction)
		width 	: width of reservoir in ft, (y direction)
		height 	: height of reservoir in ft, (z direction)

		nums3 	: number of grids, (xnums, ynums, znums)

		dims 	: dims dimension it can be 1, 2, or 3

		The object has edge properties to calculate the control
		volume implementation of flow calculations.
		"""

		xdelta = numpy.repeat(size[0]/nums3[0],nums3[0])
		ydelta = numpy.repeat(size[1]/nums3[1],nums3[1])
		zdelta = numpy.repeat(size[2]/nums3[2],nums3[2])

		super().__init__(xdelta,ydelta,zdelta,dims)

if __name__ == "__main__":

	grid = GridRegular((1000,2000,30),(3,4,1))

	print(grid.dims)

	print(grid.xmin)

