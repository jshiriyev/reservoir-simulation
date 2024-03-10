import sys

import numpy

if __name__ == "__main__":
	sys.path.append(r'C:\Users\3876yl\Documents\respy')

from respy.grids._delta import GridDelta

class GridRegular(GridDelta):

	def __init__(self,size:tuple,nums:tuple,dims:int=None,**kwargs):
		"""
		size 	: reservoir dimensions (length x width x height)

		length 	: length of reservoir in ft, (x direction)
		width 	: width of reservoir in ft, (y direction)
		height 	: height of reservoir in ft, (z direction)

		num 	: number of grids, (Nlength, Nwidth, Nheight)

		dims 	: dims dimension it can be 1, 2, or 3

		Any kargs item should follow the size regulation.
		"""

		xdelta = numpy.repeat(size[0]/nums[0],nums[0])
		ydelta = numpy.repeat(size[1]/nums[1],nums[1])
		zdelta = numpy.repeat(size[2]/nums[2],nums[2])

		super().__init__(xdelta,ydelta,zdelta,dims,**kwargs)

if __name__ == "__main__":

	grid = GridRegular((1000,2000,30),(3,4,1))

	print(grid.cube)

