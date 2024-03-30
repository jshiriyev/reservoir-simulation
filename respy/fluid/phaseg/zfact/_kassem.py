import numpy

from scipy import optimize

class DAK():
	"""
	Dranchuk-Abu-Kassem method represents the Standing-Katz correlation:
	0.2 < Pr < 15 and 0.7 < Tr < 3.0 within 1 % error
	 15 < Pr < 30 within 3% error.
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

		self.R1 = self.get_R1(self.tred)
		self.R2 = self.get_R2(self.tred)
		self.R3 = self.get_R3(self.tred)
		self.R4 = self.get_R4(self.tred)

		self.A11 = 0.7210

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
	def get_R1(tred):
		A1,A2,A3,A4,A5 = 0.3265,-1.0700,-0.5339,0.01569,-0.05165
		return A1+A2/tred+A3/tred**3+A4/tred**4+A5/tred**5

	@staticmethod
	def get_R2(tred):
		A6,A7,A8 = 0.5475,-0.7361,0.1844
		return A6+A7/tred+A8/tred**2

	@staticmethod
	def get_R3(tred):
		A7,A8,A9 = -0.7361,0.1844,0.1056
		return A9*(A7/tred+A8/tred**2)

	@staticmethod
	def get_R4(tred):
		A10 = 0.6134
		return A10/tred**3

	def __call__(self,press:numpy.ndarray):

		pred = self.pred(press)

		if numpy.any(pred>30):
			raise Warning("Dranchuk-Abu-Kassem method is not recommended.")
		elif numpy.any(pred<0.2):
			raise Warning("Dranchuk-Abu-Kassem method is not recommended.")
		elif numpy.any(pred<15) and (self.tred<0.7 or self.tred>3.0):
			raise Warning("Dranchuk-Abu-Kassem method is not recommended.")
			
		R5 = (0.27*pred)/self.tred

		residual = lambda rhor: self.zfunc(rhor)-R5/rhor
		resprime = lambda rhor: self.zprime(rhor)+R5/rhor**2

		rhor = optimize.newton(residual,R5,fprime=resprime)

		return R5/rhor,self.zprime(rhor)

	def pred(self,press):
		"""Returns reduced pressure values for input pressure values in psi."""
		return numpy.asarray(press/self.pcrit)

	def zfunc(self,rhor):
		return 1+self.R1*rhor+self.R2*rhor**2-self.R3*rhor**5\
				+self.R4*(1+self.A11*rhor**2)*rhor**2*numpy.exp(-self.A11*rhor**2)

	def zprime(self,rhor):
		return self.R1+2*self.R2*rhor-5*self.R3*rhor**4\
			+2*self.R4*rhor*(1+self.A11*rhor**2*(1-self.A11*rhor**2))*numpy.exp(-self.A11*rhor**2)

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



	