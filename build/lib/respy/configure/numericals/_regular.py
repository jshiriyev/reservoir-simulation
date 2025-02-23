import numpy

from ._time import Time

def regular(step:float,nums:int) -> Time:
	"""
	Creates a Time instance with regularly spaced time steps.

    Parameters:
        step (float) : The size of each time step.
        nums (int) 	 : The number of time steps.

    Returns:
        Time : A Time object with regularly spaced time steps.
	"""
	return Time(numpy.full(nums,step))