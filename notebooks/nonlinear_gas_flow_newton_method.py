import matplotlib.pyplot as plt

import numpy

import numpy as np

class OneDimCuboid():

    def __init__(self,length,numtot,csa):

        self.length = length

        self.numtot = numtot

        self.csa = csa

        self.num = (numtot,1,1)

        self._index()

        self._size()

        self._area()

    def _index(self):

        idx = numpy.arange(self.numtot)

        self.index = numpy.tile(idx,(7,1)).T

        self.index[idx.reshape(-1,self.num[0])[:,1:].ravel(),1] -= 1
        self.index[idx.reshape(-1,self.num[0])[:,:-1].ravel(),2] += 1
        self.index[idx.reshape(self.num[2],-1)[:,self.num[0]:],3] -= self.num[0]
        self.index[idx.reshape(self.num[2],-1)[:,:-self.num[0]],4] += self.num[0]
        self.index[idx.reshape(self.num[2],-1)[1:,:],5] -= self.num[0]*self.num[1]
        self.index[idx.reshape(self.num[2],-1)[:-1,:],6] += self.num[0]*self.num[1]

    def _size(self):

        self.size = numpy.zeros((self.numtot,3))

        self.size[:,0] = self.length/self.num[0]
        self.size[:,1] = self.csa
        self.size[:,2] = 1

        self.xaxis = numpy.linspace(0,self.length,self.numtot,endpoint=False)+self.size[:,0]/2

    def _area(self):

        self.area = numpy.zeros((self.numtot,6))

        self.area[:,0] = self.size[:,1]*self.size[:,2]
        self.area[:,1] = self.size[:,1]*self.size[:,2]
        self.area[:,2] = self.size[:,2]*self.size[:,0]
        self.area[:,3] = self.size[:,2]*self.size[:,0]
        self.area[:,4] = self.size[:,0]*self.size[:,1]
        self.area[:,5] = self.size[:,0]*self.size[:,1]

    def set_permeability(self,permeability,homogeneous=True,isotropic=True):

        self.permeability = numpy.zeros((self.numtot,3))

        if homogeneous and isotropic:
            
            self.permeability[:] = permeability

        elif homogeneous and not isotropic:

            self.permeability[:] = permeability

        elif not homogeneous and isotropic:

            self.permeability[:,0] = permeability
            self.permeability[:,1] = permeability
            self.permeability[:,2] = permeability

        else:

            self.permeability[:] = permeability

def transmissibility(grids):
    
    dx_m = (grids.size[:,0]+grids.size[grids.index[:,1],0])/2
    dx_p = (grids.size[:,0]+grids.size[grids.index[:,2],0])/2
    dy_m = (grids.size[:,1]+grids.size[grids.index[:,3],1])/2
    dy_p = (grids.size[:,1]+grids.size[grids.index[:,4],1])/2
    dz_m = (grids.size[:,2]+grids.size[grids.index[:,5],2])/2
    dz_p = (grids.size[:,2]+grids.size[grids.index[:,6],2])/2

    kx_m = (2*dx_m)/(grids.size[:,0]/grids.permeability[:,0]+
                     grids.size[grids.index[:,1],0]/grids.permeability[grids.index[:,1],0])
    kx_p = (2*dx_p)/(grids.size[:,0]/grids.permeability[:,0]+
                     grids.size[grids.index[:,2],0]/grids.permeability[grids.index[:,2],0])
    ky_m = (2*dy_m)/(grids.size[:,1]/grids.permeability[:,1]+
                     grids.size[grids.index[:,3],1]/grids.permeability[grids.index[:,3],1])
    ky_p = (2*dy_p)/(grids.size[:,1]/grids.permeability[:,1]+
                     grids.size[grids.index[:,4],1]/grids.permeability[grids.index[:,4],1])
    kz_m = (2*dz_m)/(grids.size[:,2]/grids.permeability[:,2]+
                     grids.size[grids.index[:,5],2]/grids.permeability[grids.index[:,5],2])
    kz_p = (2*dz_p)/(grids.size[:,2]/grids.permeability[:,2]+
                     grids.size[grids.index[:,6],2]/grids.permeability[grids.index[:,6],2])

    transmissibility = np.zeros((grids.numtot,6))

    transmissibility[:,0] = (kx_m*grids.area[:,0])/(grids.viscosity*dx_m)
    transmissibility[:,1] = (kx_p*grids.area[:,1])/(grids.viscosity*dx_p)
    transmissibility[:,2] = (ky_m*grids.area[:,2])/(grids.viscosity*dy_m)
    transmissibility[:,3] = (ky_p*grids.area[:,3])/(grids.viscosity*dy_p)
    transmissibility[:,4] = (kz_m*grids.area[:,4])/(grids.viscosity*dz_m)
    transmissibility[:,5] = (kz_p*grids.area[:,5])/(grids.viscosity*dz_p)

    return transmissibility

def transmatrix(grids,transmissibility):

    T = np.zeros((grids.numtot,grids.numtot))

    for i in range(grids.numtot):

        if i>0:
            T[i,i-1] = -transmissibility[i,0]
            T[i,i] += transmissibility[i,0]

        if i<grids.numtot-1:
            T[i,i+1] = -transmissibility[i,1]
            T[i,i] += transmissibility[i,1]

    return T

def newtonsolver(grids,timestep,timesteps,T,J,Q):

    array = np.zeros((grids.numtot,timesteps))

    P = grids.pressure_initial

    for j in range(timesteps):

        Pk = P.copy()

        error = 1

        k = 1

        while error>1e-6:

            A = np.diag(100/Pk.flatten())

            F = -np.matmul(T+J+A,Pk)+np.matmul(A,P)+Q

            error = np.linalg.norm(F.flatten(),2)

            JACOB = -(T+J)+np.matmul(-A,np.diag(P.flatten()/Pk.flatten()))

            Pk += np.linalg.solve(JACOB,-F)

            print(f"iteration #{k}: {error = }")
            print(f"{Pk}\n")#{F}\n

            k += 1

        P = Pk.copy()

        array[:,j] = P.flatten()

    return array

Pu = 200
Pd = 14.7

Nt = 3

dt = 1e-5

grids = OneDimCuboid(length=1,numtot=4,csa=0.02)

grids.porosity = 0.2
grids.viscosity = 0.01

grids.pressure_initial = np.ones((grids.numtot,1))*14.7

grids.set_permeability(permeability=50,homogeneous=True,isotropic=True)

Tm = transmissibility(grids)*(1.06235016e-14/1.4503774389728e-7*60*60*24)

T = transmatrix(grids,Tm)
J = np.zeros((grids.numtot,grids.numtot))
Q = np.zeros((grids.numtot,1))

for i in range(grids.numtot):

	if i==0:
		J[i,i] = 2*Tm[i,0]
		Q[i,0] = Pu*2*Tm[i,0]

	if i==grids.numtot-1:
		J[i,i] = 2*Tm[i,1]
		Q[i,0] = Pd*2*Tm[i,1]

P = newtonsolver(grids,dt,Nt,T,J,Q)

plt.scatter(grids.xaxis,grids.pressure_initial)
plt.scatter(grids.xaxis,P[:,0])
plt.scatter(grids.xaxis,P[:,1])
plt.scatter(grids.xaxis,P[:,2])

plt.show()