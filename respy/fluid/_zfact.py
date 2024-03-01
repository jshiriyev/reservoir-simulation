import numpy

from scipy import optimize

def unknown(Pr:numpy.ndarray,Tr:float):
	"""Function to Calculate Gas Compressibility Factor

	Pr	: reduced pressure
	Tr	: reduced temperatue
	"""
	Pr = numpy.asarray(Pr)

	a = 1.39*(Tr-0.92)**0.5-0.36*Tr-0.101
	b = (0.62-0.23*Tr)*Pr+(0.066/(Tr-0.86)-0.037)*Pr**2+0.32*Pr**6/(10**(9*(Tr-1)))
	c = (0.132-0.32*numpy.log(Tr)/numpy.log(10))
	d = 10**(0.3106-0.49*Tr+0.1824*Tr**2)

	return a+(1-a)*numpy.exp(-b)+c*Pr**d

HY = """
Hall Yarborough method is not recommended for application\n\
if pseudo-reduced temperature is less than one.
"""

def Hall_Yarborough(Pr:numpy.ndarray,Tr:float):
	"""
	Hall-Yarborough Method

	Pr	: reduced pressure
	Tr	: reduced temperatue
	"""

	Pr = numpy.asarray(Pr)

	if Tr<1:
		raise Warning(HY)

	X1 = -0.06125*Pr/Tr*numpy.exp(-1.2*(1-1/Tr)**2)
	X2 = 14.76/Tr-9.76/Tr**2+4.58/Tr**3
	X3 = 90.7/Tr-242.2/Tr**2+42.4/Tr**3
	X4 = 2.18+2.82/Tr

	residual = lambda Y: X1+(Y+Y**2+Y**3+Y**4)/(1-Y)**3-X2*Y**2+X3*Y**X4
	resprime = lambda Y: (1+4*Y+4*Y**2-4*Y**3+Y**4)/(1-Y)**4-2*X2*Y+X3*X4*Y**(X4-1)

	Y0 = 0.0125*Pr/Tr*numpy.exp(-1.2*(1-1/Tr)**2) # initial guess for Y

	Y = optimize.newton(residual,Y0,fprime=resprime)

	return -X1/Y

DAK = """
Dranchuk-Abu-Kassem method represents the Standing-Katz correlation:\n\
0.2 < Pr < 15 and 0.7 < Tr < 3.0 within 1 % error \n\
 15 < Pr < 30 within 3% error.
"""

def Dranchuk_Abu_Kassem(Pr:numpy.ndarray,Tr:float,*,derivative=False):
	"""
	Dranchuk-Abu-Kassem Method for calculation of z-factor.

	Pr	: reduced pressure
	Tr	: reduced temperatue
	
	If derivative is True, the derivative of z w.r.t rhor is also returned.
	"""

	Pr = numpy.asarray(Pr)

	if numpy.any(Pr>30):
		raise Warning(DAK)
	elif numpy.any(Pr<0.2):
		raise Warning(DAK)
	elif numpy.any(Pr<15) and (Tr<0.7 or Tr>3.0):
		raise Warning(DAK)
		
	A1  =  0.3265
	A2  = -1.0700
	A3  = -0.5339
	A4  =  0.01569
	A5  = -0.05165
	A6  =  0.5475
	A7  = -0.7361
	A8  =  0.1844
	A9  =  0.1056
	A10 =  0.6134
	A11 =  0.7210

	R1 = A1+A2/Tr+A3/Tr**3+A4/Tr**4+A5/Tr**5
	R2 = A6+A7/Tr+A8/Tr**2
	R3 = A9*(A7/Tr+A8/Tr**2)
	R4 = A10/Tr**3
	R5 = (0.27*Pr)/Tr

	zfunc = lambda rhor: 1+R1*rhor+R2*rhor**2-R3*rhor**5\
			+R4*(1+A11*rhor**2)*rhor**2*numpy.exp(-A11*rhor**2)

	prime = lambda rhor: R1+2*R2*rhor-5*R3*rhor**4\
			+2*R4*rhor*(1+A11*rhor**2*(1-A11*rhor**2))*numpy.exp(-A11*rhor**2)

	residual = lambda rhor: zfunc(rhor)-R5/rhor
	resprime = lambda rhor: prime(rhor)+R5/rhor**2

	rhor = optimize.newton(residual,R5,fprime=resprime)

	if derivative:
		return R5/rhor,prime(rhor)

	return R5/rhor

DPR = """
Dranchuk-Purvis-Robinson method is recommended for applications where\n\
0.2 < Pr < 3.0 and 1.05 < Tr < 3.0.
"""

def Dranchuk_Purvis_Robinson(Pr,Tr):
	"""
	Dranchuk-Purvis-Robinson Method

	Pr	: reduced pressure
	Tr	: reduced temperatue
	"""

	if Tr<1.05 or Tr>3 or Pr<0.2 or Pr>3.0:
		raise Warning(DPR)

	A1 =  0.31506237
	A2 = -1.0467099
	A3 = -0.57832720
	A4 =  0.53530771
	A5 = -0.61232032
	A6 = -0.10488813
	A7 =  0.68157001
	A8 =  0.68446549

	T1 = A1+A2/Tr+A3/Tr**3
	T2 = A4+A5/Tr
	T3 = A5*A6/Tr
	T4 = A7/Tr**3
	T5 = (0.27*Pr)/Tr

	zfunc = lambda rhor: 1+T1*rhor+T2*rhor**2+T3*rhor**5\
		+(T4*rhor**2*(1+A8*rhor**2)*math.exp(-A8*rhor**2))

	residual = lambda rhor: zfunc(rhor)-T5/rhor

	rhor = optimize.newton(residual,T5)

	return T5/rhor

if __name__=="__main__":

	import matplotlib.pyplot as plt

	import numpy as np

	Pr,Tr = 2.99,1.52

	# print(unknown(Pr,Tr))
	# print(Hall_Yarborough(Pr,Tr))
	# print(Dranchuk_Abu_Kassem(Pr,Tr))
	# print(Dranchuk_Purvis_Robinson(Pr,Tr))

	Pr = np.linspace(0.2,8,200)

	# z1 = unknown(Pr,Tr)
	# z2 = np.zeros(Pr.shape)
	z3,d3 = Dranchuk_Abu_Kassem(Pr,Tr,derivative=True)
	# z4 = np.zeros(Pr.shape)

	# plt.plot(Pr,z1,label = "unknown")
	# plt.plot(Pr,z2,label = "Hall_Yarborough")
	plt.plot(Pr,z3,label = "Dranchuk_Abu_Kassem")
	# plt.plot(Pr,z4,label = "Dranchuk_Purvis_Robinson")

	cr3 = 1/(1+(0.27*Pr)/(z3**2*Tr)*d3)/Pr

	plt.plot(Pr,cr3,label = "Reduced Compressibility")

	plt.legend()

	plt.show()



	