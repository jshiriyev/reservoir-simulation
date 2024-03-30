import numpy

from scipy import optimize

class DPR():
	"""
	Dranchuk-Purvis-Robinson Method class for calculation of z-factor;
	Recommended for applications where:
	0.2 < Pr < 3.0 and 1.05 < Tr < 3.0.
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

		self.T1 = self.get_T1(self.tred)
		self.T2 = self.get_T2(self.tred)
		self.T3 = self.get_T3(self.tred)
		self.T4 = self.get_T4(self.tred)

		self.A8 = 0.68446549

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
	def get_T1(tred):
		A1 =  0.31506237
		A2 = -1.0467099
		A3 = -0.57832720
		return A1+A2/tred+A3/tred**3

	@staticmethod
	def get_T2(tred):
		A4 =  0.53530771
		A5 = -0.61232032
		return A4+A5/tred

	@staticmethod
	def get_T3(tred):
		A5 = -0.61232032
		A6 = -0.10488813
		return A5*A6/tred

	@staticmethod
	def get_T4(tred):
		A7 = 0.68157001
		return A7/tred**3

	def __call__(self,press:numpy.ndarray):

		pred = self.pred(press)

		if self.tred<1.05 or self.tred>3 or numpy.any(pred<0.2) or numpy.any(pred>3.0):
			raise Warning("Dranchuk-Purvis-Robinson Method is not recommended")

		T5 = (0.27*pred)/self.tred

		residual = lambda rhor: self.zfunc(rhor)-T5/rhor
		resprime = lambda rhor: self.zprime(rhor)+T5/rhor**2

		rhor = optimize.newton(residual,T5,fprime=resprime)

		return T5/rhor,self.zprime(rhor)

	def pred(self,press):
		"""Returns reduced pressure values for input pressure values in psi."""
		return numpy.asarray(press/self.pcrit)

	def zfunc(self,rhor):
		return 1+self.T1*rhor+self.T2*rhor**2+self.T3*rhor**5\
				+self.T4*rhor**2*(1+self.A8*rhor**2)*numpy.exp(-self.A8*rhor**2)

	def zprime(self,rhor):
		return self.T1+2*self.T2*rhor+5*self.T3*rhor**4\
			+2*self.T4*rhor*(1+self.A8*rhor**2-self.A8**2*rhor**4)*numpy.exp(-self.A8*rhor**2)

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



	