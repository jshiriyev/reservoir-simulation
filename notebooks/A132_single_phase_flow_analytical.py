import numpy as np

class CoreFlow():

    @staticmethod
    def solve(length,time,pinit,pright,eta,ngrids):
        """
        lenght  : length of core sample
        time    : time at which to calculate pressure
        pinit   : initial pressure
        pright  : constant pressure boundary implemented at right hand side
        eta     : hydraulic diffusivity
        ngrids  : number of pressure calculation points
        """

        x = np.arange(length/ngrids/2,length,length/ngrids)
        
        X = CoreFlow.tondimx(x,length)

        T = CoreFlow.tondimtime(time,eta,length)
        
        P = CoreFlow.ndimpressure(X,T)
        
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
    def ndimpressure(ndimx:np.ndarray,ndimtime:float,terms:int=100):
        """Calculates the non-dimensional pressure given with equation 5 on the word document.

        ndimx      : non-dimensional x
        ndimtime   : non-dimensional time
        terms      : number of terms after which to terminate the summation

        """

        psum = np.zeros(ndimx.shape)

        for n in range(terms):

            fstterm = (-1)**n/(2*n+1)
            
            expterm = np.exp(-(2*n+1)**2*np.pi**2*ndimtime/4)
            costerm = np.cos((2*n+1)*np.pi/2*ndimx)
            
            psum += fstterm*expterm*costerm

        return psum*4/np.pi
