import numpy

from scipy import optimize

class OE():
	"""Orkahub Energy: Provides Direct Function to Calculate
	Gas Compressibility Factor"""

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

		self.a = self.get_a(self.tred)
		self.c = self.get_c(self.tred)
		self.d = self.get_d(self.tred)

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
	def get_a(tred):
		return 1.39*(tred-0.92)**0.5-0.36*tred-0.101

	@staticmethod
	def get_b(pred,tred):
		return (0.62-0.23*tred)*pred+(0.066/(tred-0.86)-0.037)*pred**2+0.32*pred**6/(10**(9*(tred-1)))

	@staticmethod
	def get_c(tred):
		return (0.132-0.32*numpy.log10(tred))

	@staticmethod
	def get_d(tred):
		return 10**(0.3106-0.49*tred+0.1824*tred**2)

	@staticmethod
	def get_e(pred,tred):
		# e is defined as the derivative of b w.r.t Pr
		return (0.62-0.23*tred)+(0.132/(tred-0.86)-0.074)*pred+1.92*pred**5/(10**(9*(tred-1)))

	def __call__(self,press:numpy.ndarray):

		pred = self.pred(press)

		b = self.get_b(pred,self.tred)
		e = self.get_e(pred,self.tred)

		return self.zfunc(pred,b),self.zprime(pred,b,e)

	def pred(self,press):
		"""Returns reduced pressure values for input pressure values in psi."""
		return numpy.asarray(press/self.pcrit)

	def zfunc(self,pred,b):
		return self.a+(1-self.a)*numpy.exp(-b)+self.c*pred**self.d

	def zprime(self,pred,b,e):
		return (self.a-1)*numpy.exp(-b)*e+self.c*self.d*pred**(self.d-1)

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



	