import numpy

class Time():
	"""Time class for the reservoir simulation"""

	def __init__(self,steps:numpy.ndarray):
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
		self._steps = numpy.asarray(value).flatten().astype(numpy.float64)*(24*60*60)

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
		self._total = numpy.sum(self._steps)

	def __iter__(self):
		"""
		It iterates over the time steps, and yields:

		index (int) 		: The time iteration index
		curr_time (float) 	: The current time of the iteration
		step_size (float) 	: The step size of the iteration
		"""
		for index,(curr_time,step_size) in enumerate(zip(self._times[:-1],self._steps)):

			yield index,curr_time,step_size

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
		self._times = numpy.insert(numpy.cumsum(self._steps),0,0)
    

if __name__ == "__main__":

	t = Time((1,1,1,1,1,1,1,1,1,1))

	print(t.times)
	print(t.steps)
	print(t.total)

	print(t.nums)

	for a in t:
		print(a[0],a[1]/86400,a[2]/86400,t.times[a[0]+1])

	print(t.times.size)