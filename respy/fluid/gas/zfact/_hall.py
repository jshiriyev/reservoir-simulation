import numpy

from scipy import optimize

class HY():
	"""
	Hall Yarborough method for z factor calculation; not recommended
	for if pseudo-reduced temperature is less than one.
	"""

	def __init__(self,crit,temp):
		"""Z factor class that can be used for ResPy, returning compressibility
			factors for pressures when called. Initialization parameters are:

		crit    : tuple of (pcrit in psi, tcrit in Rankine)

		temp    : temperature value (Rankine) which will be used to calculate
				  fluid properties when the class is called.
		"""

		pcrit,tcrit = crit

		self._pcrit = pcrit*6894.76
		self._tcrit = tcrit*(5./9)

		self._temp = temp*(5./9)

		if self.tred<1:
			raise Warning("Hall Yarborough method is not recommended for Tr less than one.")

		self.X2 = self.get_X2(self.tred)
		self.X3 = self.get_X3(self.tred)
		self.X4 = self.get_X4(self.tred)

	@property
	def pcrit(self):
		"""Critical Pressure in psi"""
		return self._pcrit/6894.76

	@property
	def tcrit(self):
		"""Critical Temperature in Rankine"""
		return self._tcrit*(9./5)

	@property
	def temp(self):
		return self._temp*(9./5)

	@property
	def tred(self):
		"""Returns reduced temperature (class property)."""
		return self.temp/self.tcrit

	@staticmethod
	def get_X2(tred):
		return 14.76/tred-9.76/tred**2+4.58/tred**3

	@staticmethod
	def get_X3(tred):
		return 90.7/tred-242.2/tred**2+42.4/tred**3

	@staticmethod
	def get_X4(tred):
		return 2.18+2.82/tred

	def __call__(self,press:numpy.ndarray):

		pred = self.pred(press)

		X1 = 0.06125*pred/self.tred*numpy.exp(-1.2*(1-1/self.tred)**2)

		residual = lambda Y: self.zfunc(Y)*Y-X1
		resprime = lambda Y: self.zfunc(Y)+self.zprime(Y)*Y

		Y0 = 0.0125*pred/self.tred*numpy.exp(-1.2*(1-1/self.tred)**2) # initial guess for Y

		Y = optimize.newton(residual,Y0,fprime=resprime)

		return X1/Y,self.zprime(Y)

	def pred(self,press):
		"""Returns reduced pressure values for input pressure values in psi."""
		return numpy.asarray(press/self.pcrit)

	def zfunc(self,Y):
		return (1+Y+Y**2+Y**3)/(1-Y)**3-self.X2*Y+self.X3*Y**(self.X4-1)

	def zprime(self,Y):
		return 4*(1+Y+Y**2)/(1-Y)**4-self.X2+self.X3*(self.X4-1)*Y**(self.X4-2)

if __name__=="__main__":

	import matplotlib.pyplot as plt

	import numpy as np

	Pr,Tr = 2.99,1.52

	# print(Orkahub_Energy(Pr,Tr,derivative=True))
	# print(Hall_Yarborough(Pr,Tr,derivative=True))
	# print(Dranchuk_Abu_Kassem(Pr,Tr,derivative=True))
	# print(Dranchuk_Purvis_Robinson(Pr,Tr,derivative=True))

	Pr = np.linspace(0.2,3,200)

	z1,d1 = Orkahub_Energy(Pr,Tr,derivative=True)
	z2,d2 = Hall_Yarborough(Pr,Tr,derivative=True)
	z3,d3 = Dranchuk_Abu_Kassem(Pr,Tr,derivative=True)
	z4,d4 = Dranchuk_Purvis_Robinson(Pr,Tr,derivative=True)

	# plt.plot(Pr,z1,label = "Orkahub_Energy")
	# plt.plot(Pr,z2,label = "Hall_Yarborough")
	# plt.plot(Pr,z3,label = "Dranchuk_Abu_Kassem")
	# plt.plot(Pr,z4,label = "Dranchuk_Purvis_Robinson")

	cr1 = 1/(1+(0.27*Pr)/(z1**2*Tr)*d1)/Pr
	cr2 = 1/(1+(0.27*Pr)/(z2**2*Tr)*d2)/Pr
	cr3 = 1/(1+(0.27*Pr)/(z3**2*Tr)*d3)/Pr
	cr4 = 1/(1+(0.27*Pr)/(z4**2*Tr)*d4)/Pr

	plt.plot(Pr,cr1,label = "Orkahub_Energy")
	plt.plot(Pr,cr2,label = "Hall_Yarborough")
	plt.plot(Pr,cr3,label = "Dranchuk_Abu_Kassem")
	plt.plot(Pr,cr4,label = "Dranchuk_Purvis_Robinson")

	plt.legend()

	plt.show()



	