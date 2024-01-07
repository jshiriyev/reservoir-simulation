import numpy as np

from scipy.sparse import csr_matrix as csr
from scipy.sparse import diags

from scipy.sparse.linalg import spsolve as sps

class CoreFlow():

    @staticmethod
    def solve(length,time,pinit,pright,eta,ngrids,nsteps):
        """
        lenght  : length of core sample
        time    : time at which to calculate pressure
        pinit   : initial pressure
        pright  : constant pressure boundary implemented at right hand side
        eta     : hydraulic diffusivity
        ngrids  : number of pressure calculation points
        nsteps  : number of time steps to be used in finite difference solution
        """

        x = np.arange(length/ngrids/2,length,length/ngrids)

        T = CoreFlow.tondimtime(time,eta,length)

        P = CoreFlow.ndimpressure(ngrids,T/nsteps,nsteps)
        
        P = CoreFlow.todimpressure(P,pinit,pright)

        return x,P
    
    @staticmethod
    def tondimx(x,length):
        """Converts dimensional x to non-dimensional x."""

        return x/length

    @staticmethod
    def todimx(x,length):
        """Converts non-dimensional x to dimensional x."""

        return x*length

    @staticmethod
    def tondimtime(time,eta,length):
        """Converts dimensional time to non-dimensional time."""

        return eta*time/length**2

    @staticmethod
    def todimtime(time,eta,length):
        """Converts non-dimensional time to dimensional time."""
        
        return time*length**2/eta

    @staticmethod
    def tondimpressure(pressure,pinit,pright):
        """Converts dimensional pressure to non-dimensional pressure."""

        return (pressure-pright)/(pinit-pright)

    @staticmethod
    def todimpressure(pressure,pinit,pright):
        """Converts non-dimensional pressure to dimensional pressure."""

        return pressure*(pinit-pright)+pright

    @staticmethod
    def ndimpressure(ngrids:int,deltat:int,nsteps:int):
        """Calculates the non-dimensional pressure with finite difference method.

        """

        shape = (ngrids,ngrids)

        deltax = 1/ngrids

        Gmatrix = csr(shape)

        indices = np.array([i for i in range(ngrids-1)],dtype="int16")

        coeffis = [deltat/deltax**2 for _ in range(ngrids-1)]

        Gmatrix -= csr((coeffis,(indices,indices+1)),shape=shape)
        Gmatrix += csr((coeffis,(indices,indices)),shape=shape)
        
        Gmatrix -= csr((coeffis,(indices+1,indices)),shape=shape)
        Gmatrix += csr((coeffis,(indices+1,indices+1)),shape=shape)

        Gmatrix += csr(([2*deltat/deltax**2],([ngrids-1],[ngrids-1])),shape=shape)

        Gmatrix += diags([1 for _ in range(ngrids)],shape=shape)

        pinit = np.ones((ngrids,))

        pn = pinit

        for n in range(nsteps):

            pn = sps(Gmatrix,pn)

        return pn
