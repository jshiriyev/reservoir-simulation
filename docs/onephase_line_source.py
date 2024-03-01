import matplotlib.pyplot as plt

import _setup

mu = 0.72   # cp

kx = 0.1    # md
ky = 0.1    # md
kz = 0.1    # md

ct = 1.5e-5 # 1/psi

pi = 3000   # psi

re = 3000   # ft
rw = 0.5    # ft

h = 150     # ft

q = 20      # STB/day

Bo = 1.475  # bbl/STB

phi = 0.23

time = 3    # hours

# UNIT CONVERSION TO SI UNITS

mu = mu*0.001

kx = kx*0.986923E-15
ky = ky*0.986923E-15
kz = kz*0.986923E-15

ct = ct/6894.76

pi = pi*6894.76

re = re*0.3048
rw = rw*0.3048

h = h*0.3048
q = q*Bo*1.85185E-6

time = time*3600

rpgreen = rectangular_parallelepiped(pi,2*re,2*re,h,phi,ct,mu,kx,ky,kz,q,time)

rpgreen.set_observers()

pressure = rpgreen.line(re,re)

rpgreen.x = rpgreen.x/0.3048
pressure = pressure/6894.76

plt.plot(rpgreen.x,pressure)

plt.show()