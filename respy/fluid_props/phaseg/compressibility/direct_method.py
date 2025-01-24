import numpy

from scipy import optimize

class DirectMethod():
	"""Direct Method: Provides direct function to calculate
	gas compressibility fractor (z-factor);
	"""

	def __init__(self,crit,temp):
		"""Z factor class that can be used for ResPy, returning compressibility
		factors for pressures when called. Initialization parameters are:

		crit    : tuple of (pcrit in psi, tcrit in Rankine)

		temp    : temperature value (Rankine) which will be used to calculate
		          fluid properties when the class is called.

		- Initialize the class with critical pressure and temperature,
		- Provide temperature at which properties will be calculated during the iterations,
		assuming isothermal conditions
		- Call the class to calculate z-factor and its prime at different pressure values.

		"""

		pcrit,tcrit = crit

		self.pcrit = pcrit
		self.tcrit = tcrit

		self.temp = temp

		self.a = self.get_a(self.treduced)
		self.c = self.get_c(self.treduced)
		self.d = self.get_d(self.treduced)

	@property
	def pcrit(self):
		"""Critical Pressure in psi, underscore parameter is in SI units."""
		return self._pcrit/6894.76

	@pcrit.setter
	def pcrit(self,value:float):
		"""Critical Pressure in SI units, Pascal."""
		self._pcrit = value*6894.76

	@property
	def tcrit(self):
		"""Critical Temperature in Rankine, underscore parameter is in SI units."""
		return self._tcrit*(9./5)

	@tcrit.setter
	def tcrit(self,value:float):
		"""Critical Temperature in SI units, Kelvin."""
		self._tcrit = value*(5./9)

	@property
	def temp(self):
		"""Temperature in Rankine, underscore parameter is in SI units."""
		return self._temp*(9./5)

	@temp.setter
	def temp(self,value:float):
		"""Temperature in SI units, Kelvin."""
		self._temp = value*(5./9)

	@property
	def treduced(self):
		"""Returns reduced temperature (class property)."""
		return self.temp/self.tcrit

	def preduced(self,press:numpy.ndarray):
		"""Returns reduced pressure values for input pressure values in psi."""
		return numpy.asarray(press)/self.pcrit

	@staticmethod
	def get_a(treduced):
		return 1.39*(treduced-0.92)**0.5-0.36*treduced-0.101

	@staticmethod
	def get_b(preduced,treduced):
		return (0.62-0.23*treduced)*preduced+(0.066/(treduced-0.86)-0.037)*preduced**2+0.32*preduced**6/(10**(9*(treduced-1)))

	@staticmethod
	def get_c(treduced):
		return (0.132-0.32*numpy.log10(treduced))

	@staticmethod
	def get_d(treduced):
		return 10**(0.3106-0.49*treduced+0.1824*treduced**2)

	@staticmethod
	def get_e(preduced,treduced):
		# e is defined as the derivative of b w.r.t Pr
		return (0.62-0.23*treduced)+(0.132/(treduced-0.86)-0.074)*preduced+1.92*preduced**5/(10**(9*(treduced-1)))

	def __call__(self,press:numpy.ndarray,derivative:bool=False):

		preduced = self.preduced(press)

		b = self.get_b(preduced,self.treduced)
		e = self.get_e(preduced,self.treduced)

		if derivative:
			return self.zvalue(preduced,b),self.zprime(preduced,b,e)

		return self.zvalue(preduced,b)

	def zvalue(self,preduced,b):
		"""Internal function to calculate z factor when the class is called."""
		return self.a+(1-self.a)*numpy.exp(-b)+self.c*preduced**self.d

	def zprime(self,preduced,b,e):
		"""Internal function to calculate z prime when the class is called."""
		return (self.a-1)*numpy.exp(-b)*e+self.c*self.d*preduced**(self.d-1)