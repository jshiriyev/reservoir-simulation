import numpy

class Time():

	def __init__(self,step:float,total:float=None,nstep:int=1):
		"""Defining the numerical time settings

		step    : time step defined in days
		total   : total simulation time defined in days
		nstep   : number of total time steps
		"""

		self._step = step*24*60*60

		if total is None:
			self._total = self._tstep*nstep
			self._nstep = nstep
		else:
			self._total = total*24*60*60
			self._nstep = int(self._total/self._step)

		self._time = self.uniform()

	def uniform(self):
		return numpy.arange(
			self._step,self._total+self._step/2,self._step)

	@property
	def step(self):
		return self._step/(24*60*60)

	@property
	def total(self):
		return self._ttime/(24*60*60)

	@property
	def nstep(self):
		return self._nstep
	
	@property
	def time(self):
		return self._time/(24*60*60)

if __name__ == "__main__":

	t = Time(0.1,10)

	print(t.step)