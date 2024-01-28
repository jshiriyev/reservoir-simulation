import numpy as np
import matplotlib.pyplot as plt

from scipy.special import expi

N = 10000

r = np.linspace(0.25,1000,N)

k = 80

co, cw, cf = 5e-6, 3e-6, 3e-6
So, Sw = 0.75, 0.25

ct = co*So+cw*Sw+cf

h = 50

B = 1.2

rw = 0.25

mu = 3

phi = 0.18

Pi = 5000

q = -300*B

t1 = 1
t2 = 10
t3 = 100

p1 = Pi-70.6*(q*mu)/(k*h)*expi(-39.5*(phi*mu*ct*r**2)/(k*t1))
p2 = Pi-70.6*(q*mu)/(k*h)*expi(-39.5*(phi*mu*ct*r**2)/(k*t2))  
p3 = Pi-70.6*(q*mu)/(k*h)*expi(-39.5*(phi*mu*ct*r**2)/(k*t3))

plt.plot(r,p1,'b',label='P(1 day)')
plt.plot(r,p2,'g',label='P(10 days)')
plt.plot(r,p3,'r',label='P(100 days)')

plt.legend()

plt.xlabel('distance [ft]')
plt.ylabel('pressure [psi]')

plt.grid()

plt.show()
