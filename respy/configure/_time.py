import numpy

class Time():
	"""Time class for the reservoir simulation"""
	def __init__(self,delta:float,total:float):
		"""Defining the time settings for the simulator

		delta   : first time step defined in days
		total   : total simulation time defined in days
		"""

		self.delta = delta
		self.total = total

		self.times = None
		self.steps = None

	@property
	def delta(self):
		return self._delta/(24*60*60)

	@delta.setter
	def delta(self,value):
		self._delta = value*24*60*60

	@property
	def total(self):
		return self._total/(24*60*60)

	@total.setter
	def total(self,value):
		self._total = value*24*60*60

	@property
	def times(self):
		return self._times/(24*60*60)

	@times.setter
	def times(self,value):
		"""Sets linearly spaced time data."""
		self._times = np.arange(
			self._total+self._delta/2,step=self._delta)

	@property
	def steps(self):
		return self._steps/(24*60*60)

	@steps.setter
	def steps(self,value):
		"""Sets the time steps in the given time array."""
		self._steps = self._times[1:]-self._times[:-1]

	def __iter__(self):
		"""It starts from 0 time and iterates till the last time step."""

		zipped = zip(self._times,self._steps)

		for index,(time,step) in enumerate(zipped):

			yield index,time,step

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