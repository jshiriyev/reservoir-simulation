import numpy as np

class JFunction:
    """
    The Leverett J-function is a way to nondimensionalize capillary pressure curves.
    It can applied to determine capillary pressure in a similar rock type
    but different permeability, porosity, or interfacial tension of the fluids.

    """

    def __init__(self,cpress,perm,poro,inten,angle=0.):
        """Initializes the J-function class.

        Parameters:
        ----------
        cpress : capillary pressure values given for the following properties

        perm   : permeability
        poro   : porosity
        inten  : interfacial tension
        angle  : contact angle

        """
        self.cpress = cpress

        self.perm = perm
        self.poro = poro

        self.inten = inten
        self.angle = angle

        self.jvals = None

    @property
    def jvals(self):
        """Getter for the j-values."""
        return self._jvals

    @jvals.setter
    def jvals(self,value):
        """Setter for the j-values."""
        self._jvals = self.get_jvals(self.cpress,self.perm,self.poro,self.inten,self.angle)
    
    def get(self,*args,**kwargs):
        """Get new capillary pressure values for the new rock and fluid properties."""
        return self.get_cpress(self.jvals,*args,**kwargs)

    @staticmethod
    def get_jvals(cpress,perm,poro,inten,angle=0.):
        """Calculates J-function from capillary pressure."""
        return cpress*np.sqrt(perm/poro)/(inten*np.cos(angle))

    @staticmethod
    def get_cpress(jvals,perm,poro,inten,angle=0.):
        """Calculates capillary pressure from J-function."""
        return jvals*(inten*np.cos(angle))/np.sqrt(perm/poro)

if __name__ == "__main__":

    from _brooks_corey import BrooksCorey

    ow = BrooksCorey(0.1,0.4,2,3.5)

    pc1 = ow.drain.pc(0.3)
    pc2 = ow.imbibe.pc(0.3)

    jf = JFunction(pc1,100,0.2,20)
    print(pc1,jf.get(50,0.15,25))

    jf = JFunction(pc2,100,0.2,20)
    print(pc2,jf.get(50,0.15,25))