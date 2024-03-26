import numpy

class Time():

	def __init__(self,start:float,total:float,**kwargs):
		"""Defining the time settings for the simulator

		start   : first time step defined in days
		total   : total simulation time defined in days
		"""

		self._start = start*24*60*60
		self._total = total*24*60*60

		self.set_times(**kwargs)

	def set_times(self,method="linspace",**kwargs):
		"""Sets time steps for reservoir simulation."""
		self._times = getattr(self,method)(**kwargs)

		self._steps = self.get_steps(self._times)

	def linspace(self,nums=None):
		"""Returns linearly spaced time data."""
		if nums is None:
			return numpy.arange(
				self._start,self._total+self._start/2,self._start)

		return numpy.linspace(
			self._start,self._total,nums)

	def logspace(self,nums=50):
		"""Returns logarithmically spaced time data."""
		return numpy.logspace(
			numpy.log10(self._start),numpy.log10(self._total),nums)

	@property
	def start(self):
		return self._start/(24*60*60)

	@property
	def total(self):
		return self._total/(24*60*60)

	@property
	def times(self):
		return self._times/(24*60*60)

	@property
	def steps(self):
		return self._steps/(24*60*60)

	@property
	def nums(self):
		return self._times.size

	"""Static methods:"""
	
	@staticmethod
	def get_steps(times):
		"""Returns the time steps in the given time array."""
		prev_times = numpy.roll(times,1)
		prev_times[0] = 0
		return times-prev_times

if __name__ == "__main__":

	t = Time(0.1,10,method='logspace')

	print(t.steps)

	print(t.steps.sum())