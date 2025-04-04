import numpy as np

class Time():
	"""Time class for the reservoir simulation"""

	def __init__(self,steps:np.ndarray):
		"""
		Initializes the time settings for the simulator.
		
		Parameters
		----------
		steps   : time steps in days.
		"""
		self.steps = steps
		self.total = None
		self.times = None
		
	@property
	def steps(self):
		"""Getter for time steps"""
		return self._steps/(24*60*60)

	@steps.setter
	def steps(self,value):
		"""Setter for times steps"""
		self._steps = np.asarray(value).flatten().astype(np.float64)*(24*60*60)

	@property
	def nums(self):
		"""Getter for the number of time steps"""
		return self._steps.size

	@property
	def total(self):
		"""Getter for the total simulation time"""
		return self._total/(24*60*60)

	@total.setter
	def total(self,value):
		"""Setter for the total simulation time"""
		self._total = np.sum(self._steps)

	def __iter__(self):
		"""
		It iterates over the time steps, and yields:

		curr_time (float) 	: The current time of the iteration
		step_size (float) 	: The step size of the iteration
		"""
		for curr_time,step_size in zip(self._times[:-1],self._steps):

			yield curr_time,step_size

	def __repr__(self):
		"""Representation of Time class instances"""
		return (f"Time(steps={self.steps}, total={self.total})")

	@property
	def times(self):
		"""Getter for times array"""
		return self._times/(24*60*60)

	@times.setter
	def times(self,value):
		"""Setter for times array"""
		self._times = np.insert(np.cumsum(self._steps),0,0)

def bysteps(steps:np.ndarray) -> Time:
	"""
	Creates a Time instance with irregularly spaced time steps.

	Parameters:
        steps (array) : The size of time step for the simulation.

    Returns:
        Time : A Time object with regularly spaced time steps.

	"""
	return Time(np.asarray(steps))

def regular(step:float,nums:int) -> Time:
	"""
	Creates a Time instance with regularly spaced time steps.

    Parameters:
        step (float) : The size of each time step.
        nums (int) 	 : The number of time steps.

    Returns:
        Time : A Time object with regularly spaced time steps.

	"""
	return Time(np.full(nums,step))

if __name__ == "__main__":

	t = Time((1,1,1,1,1,1,1,1,1,1))

	print(t.times)
	print(t.steps)
	print(t.total)

	print(t.nums)

	for a in t:
		print(a[0],a[1]/86400,a[2]/86400,t.times[a[0]+1])

	print(t.times.size)