"""
THIS IS THE ONE TO RUN

MAKE SURE TO HAVE
single_phase_flow_analytical.py
single_phase_flow_numerical.py

in the same directory.
"""


import matplotlib.pyplot as plt

from A132_single_phase_flow_analytical import CoreFlow as analytical
from A132_single_phase_flow_numerical import CoreFlow as numerical

k = 10      #mD
phi = 0.2   #ndim
mu = 0.1    #cp
ct = 1e-5   #1/psi

coeff = 6894.76*9.869233e-16/1e-3 # conversion to SI unit

eta = k/(phi*mu*ct)*coeff # m2/second

L = 2       #m

Pi = 689476 #Pascal
Pr = 101325 #Pascal

t_1 = 0.001 #second
t_2 = 0.010 #second
t_3 = 0.050 #second
t_4 = 0.100 #second
t_5 = 0.300 #second
t_6 = 0.500 #second
t_7 = 1.000 #second

xa,P1a = analytical.solve(L,t_1,Pi,Pr,eta,200)
# xa,P2a = analytical.solve(L,t_2,Pi,Pr,eta,200)
xa,P3a = analytical.solve(L,t_3,Pi,Pr,eta,200)
# xa,P4a = analytical.solve(L,t_4,Pi,Pr,eta,200)
# xa,P5a = analytical.solve(L,t_5,Pi,Pr,eta,200)
# xa,P6a = analytical.solve(L,t_6,Pi,Pr,eta,200)
xa,P7a = analytical.solve(L,t_7,Pi,Pr,eta,200)

xn,P1n = numerical.solve(L,t_1,Pi,Pr,eta,20,10)
# xn,P2n = numerical.solve(L,t_2,Pi,Pr,eta,20,10)
xn,P3n = numerical.solve(L,t_3,Pi,Pr,eta,20,10)
# xn,P4n = numerical.solve(L,t_4,Pi,Pr,eta,20,10)
# xn,P5n = numerical.solve(L,t_5,Pi,Pr,eta,20,10)
# xn,P6n = numerical.solve(L,t_6,Pi,Pr,eta,20,10)
xn,P7n = numerical.solve(L,t_7,Pi,Pr,eta,20,10)

plt.plot(xa,P1a/6894.76,color='black',linewidth=0.9,linestyle='solid')
# plt.plot(xa,P2a/6894.76,color='black',linewidth=0.8,linestyle=(0,(5,1)))
plt.plot(xa,P3a/6894.76,color='black',linewidth=0.7,linestyle=(0,(5,1)))
# plt.plot(xa,P4a/6894.76,color='black',linewidth=0.6,linestyle='dashed')
# plt.plot(xa,P5a/6894.76,color='black',linewidth=0.5,linestyle='dashed')
# plt.plot(xa,P6a/6894.76,color='black',linewidth=0.4,linestyle='dashed')
plt.plot(xa,P7a/6894.76,color='black',linewidth=0.3,linestyle='dashed')

plt.scatter(xn,P1n/6894.76,color='red')
# plt.scatter(xn,P2n/6894.76,color='red')
plt.scatter(xn,P3n/6894.76,color='red')
# plt.scatter(xn,P4n/6894.76,color='red')
# plt.scatter(xn,P5n/6894.76,color='red')
# plt.scatter(xn,P6n/6894.76,color='red')
plt.scatter(xn,P7n/6894.76,color='red')

plt.ylabel('Pressure [psi]')
plt.xlabel('x-axis [m]')

plt.show()
