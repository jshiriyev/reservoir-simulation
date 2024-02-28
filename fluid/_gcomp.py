import math

A1  = 0.3265
A2  = -1.0700
A3  = -0.5339
A4  = 0.01569
A5  = -0.05165
A6  = 0.5475
A7  = -0.7361
A8  = 0.1844
A9  = 0.1056
A10 = 0.6134
A11 = 0.7210

def get_gcomp(P,T,Pc,Tc,Z):
    """This method should not be used at Tr<1.4 for 0.4<Pr<3.0"""

    Pr,Tr = P/Pc,T/Tc
    
    rhor = 0.27*(Pr)/(Z*Tr)

    prime = A1+A2/Tr+A3/Tr**3+A4/Tr**4
    prime += A5/Tr**5+2*rhor*(A6+A7/Tr+A8/Tr**2)
    prime -= 5*rhor**4*A9*(A7/Tr+A8/Tr**2)
    prime += 2*A10*rhor/Tr**3*(1+A11*rhor**2-A11**2*rhor**4)*math.exp(-A11*rhor**2)

    return 1/pres-0.27/(Z**2*Tr*Pc)*(prime/(1+rhor/Z*prime))