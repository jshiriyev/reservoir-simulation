import numpy as np

class BrooksCorey():

    def __init__(self,swr:float=0.,sor:float=0.,lamda:float=2.,entry:float=4.5):
        """
        BrooksCorey capillary  pressure model.

        Parameters:
        ----------
        swr     : irreducible water saturation
        sor     : residual oil saturation

        lamda   : empirical model constant
        entry   : capillary entry pressure

        """
        self.swr = swr
        self.sor = sor

        self.lamda = lamda
        self.entry = entry

    @property
    def drain(self):
        """Returns an instance of the Drainage class."""
        return Drainage(self.swr,self.sor,self.lamda,self.entry)

    @property
    def imbibe(self):
        """Returns an instance of the Imbibition class."""
        return Imbibition(self.swr,self.sor,self.lamda,self.entry)

class Drainage(BrooksCorey):
    """
    Drainage process in the Brooks-Corey model.
    """

    def __init__(self,*args,**kwargs):

        super().__init__(*args,**kwargs)

    def pc(self,sw:np.ndarray):
        """Calculates capillary pressure from saturation.

        sw     : water saturation values

        """
        ss = (np.ravel(sw)-self.swr)/(1-self.swr)
        ss = np.clip(ss,0,1)

        return self.entry*ss**(-1/self.lamda)

    def sw(self,pc:np.ndarray):
        """Calculates saturation from capillary pressure, inverse calculations.
        
        pc     : capillary pressure values

        """
        ss = (self.entry/np.ravel(pc))**(self.lamda)

        return self.swr+ss*(1-self.swr)

class Imbibition(BrooksCorey):
    """
    Imbibition process in the Brooks-Corey model.
    """

    def __init__(self,*args,**kwargs):

        super().__init__(*args,**kwargs)

    def pc(self,sw:np.ndarray):
        """Calculates capillary pressure from saturation.

        sw     : water saturation values

        """
        se = (np.ravel(sw)-self.swr)/(1-self.swr-self.sor)
        se = np.clip(se,0,1)

        return self.entry*(se**(-1/self.lamda)-1)

    def sw(self,pc:np.ndarray):
        """Calculates saturation from capillary pressure, inverse calculations.

        pc     : capillary pressure values

        """
        se = (np.ravel(pc)/self.entry+1)**(-self.lamda)

        return self.swr+se*(1-self.swr-self.sor)

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    ow = BrooksCorey(0.1,0.2,2,3.5)

    sw = np.linspace(0,1,1000)

    plt.semilogy(sw,ow.drain.pc(sw),label="Drainage Curve")
    plt.semilogy(sw,ow.imbibe.pc(sw),label="Imbibition Curve")

    plt.vlines(0.1,1e-3,1e3,'k',linestyle='--',linewidth=1.)
    plt.vlines(0.8,1e-3,1e3,'k',linestyle='--',linewidth=1.)

    plt.xlim((0,1.))
    plt.ylim((1e-3,1e3))

    plt.legend()

    plt.show()