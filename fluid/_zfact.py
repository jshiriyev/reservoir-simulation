"""
Hall-Yarborough
Dranchuk-Abu-Kassem
Dranchuk-Purvis-Robinson
"""

from scipy import optimize

def unknown(Pr,Tr):
	"""Function to Calculate Gas Compressibility Factor
	Pr         reduced pressure
	Tr         reduced temperatue
	"""

	a = 1.39*(Tr-0.92)**0.5-0.36*Tr-0.101
	b = (0.62-0.23*Tr)*Pr+(0.066/(Tr-0.86)-0.037)*Pr**2+0.32*Pr**6/(10**(9*(Tr-1)))
	c = (0.132-0.32*math.log(Tr)/math.log(10))
	d = 10**(0.3106-0.49*Tr+0.1824*Tr**2)

	return a+(1-a)*math.exp(-b)+c*Pr**d

def HY():

	pass

def DAK():

	rhor = 0.27*(Pr)/(Z*Tr)

	1
	+(A1+A2/Tr+A3/Tr**3+A4/Tr**4+A5/Tr**5)*rhor
	+(A6+A7/Tr+A8/Tr**2)*rhor**2
	-A9*(A7/Tr+A8/Tr**2)*rhor**5
	+A10*(1+A11*rhor**2)*(rhor**2/Tr**3)*math.exp(-A11*rhor**2)

def DPR():

	return