import matplotlib.pyplot as plt

import numpy as np

from borepy.pormed import OnePhase

from borepy.gmodel.mesh import OneDimCuboid

Pu = 200
Pd = 14.7

Nt = 3

dt = 1e-5

grids = OneDimCuboid(length=1,numtot=4,csa=0.02)

grids.porosity = 0.2
grids.viscosity = 0.01

grids.pressure_initial = np.ones((grids.numtot,1))*14.7

grids.set_permeability(permeability=50,homogeneous=True,isotropic=True)

Tm = OnePhase.transmissibility(grids)*(1.06235016e-14/1.4503774389728e-7*60*60*24)

T = OnePhase.transmatrix(grids,Tm)
J = np.zeros((grids.numtot,grids.numtot))
Q = np.zeros((grids.numtot,1))

for i in range(grids.numtot):

	if i==0:
		J[i,i] = 2*Tm[i,0]
		Q[i,0] = Pu*2*Tm[i,0]

	if i==grids.numtot-1:
		J[i,i] = 2*Tm[i,1]
		Q[i,0] = Pd*2*Tm[i,1]

P = OnePhase.picardsolver(grids,dt,Nt,T,J,Q)

plt.scatter(grids.xaxis,grids.pressure_initial)
plt.scatter(grids.xaxis,P[:,0])
plt.scatter(grids.xaxis,P[:,1])
plt.scatter(grids.xaxis,P[:,2])

plt.show()