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



	