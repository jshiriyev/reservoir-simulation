import numpy

class Time():
	"""Time class for the reservoir simulation"""

	def __init__(self,delta:float,total:float):
		"""Defining the time settings for the simulator

		delta   : first time step defined in days
		total   : total simulation time defined in days
		"""

		self._delta = delta*24*60*60
		self._total = total*24*60*60

		self._times = self.get_times()
		self._steps = self.get_steps()

	def get_times(self):
		"""Returns linearly spaced time data."""
		return numpy.arange(
			self._total+self._delta/2,step=self._delta)

	def get_steps(self):
		"""Returns the time steps in the given time array."""
		return self._times[1:]-self._times[:-1]

	def __iter__(self):
		"""It starts from 0 time and iterates till the last time step."""

		zipped = zip(self._times,self._steps)

		for index,(time,step) in enumerate(zipped):

			yield index,time,step

	@property
	def delta(self):
		return self._delta/(24*60*60)

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
		return self._steps.size

if __name__ == "__main__":

	t = Time(1,10)

	print(t.times)
	print(t.steps)

	print(t.nums)

	for a in t:
		print(a[0],[A/(24*60*60) for A in a[1:]])