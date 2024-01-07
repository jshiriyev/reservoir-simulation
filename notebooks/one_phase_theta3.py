import matplotlib.pyplot as plt

import _setup

X = np.linspace(0,10,10000,dtype=np.float64)
T = np.linspace(0,0.2,20,dtype=np.float64)[1:]

theta = elliptictheta3(X,T,Nmax=1000)

iterator1 = theta.function()
iterator2 = theta.integral()

for t in T:

    print("Current time step is {}".format(t))

    func = next(iterator1)
    intg = next(iterator2)

plt.plot(X,func,label="function")
plt.plot(X,intg,label="integral")

plt.legend()

plt.show()