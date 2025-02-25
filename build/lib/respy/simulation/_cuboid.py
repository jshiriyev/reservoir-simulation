import logging

import numpy

class Cuboid():
    """Rectangular Cuboid class"""

    def __init__(self,grids,rrock,fluid,tcomp=None):
        """Initialization of block (cell) calculation class."""
        self.grids = grids
        self.rrock = rrock
        self.fluid = fluid

        self._tcomp,self.tcomp = None,tcomp

    def __getattr__(self,key):
        """Delegates attribute access to the underlying grid object."""
        return getattr(self.grids,key)

    def __call__(self,rrock=None,fluid=None,tcomp=None):
        """Returns the instance with updated rock and fluid properties."""
        self.rrock = rrock # reservoir rock properties
        self.fluid = fluid # reservoir fluid properties
        self.tcomp = tcomp # total compressibility

    @property
    def rrock(self):
        return self._rrock
    
    @rrock.setter
    def rrock(self,value):
        """Sets the reservoir rock permeability and transmissibility into the cell."""
        if value is None:
            return

        self._rrock = value

        # Reset flow properties to indicate they need recalculation
        self.xflow = None # Block transmissibility in x-direction
        self.yflow = None # Block transmissibility in y-direction
        self.zflow = None # Block transmissibility in z-direction

    @property
    def xflow(self):
        """Getter for the rock transmissibility in x-direction."""
        return self._xflow/(9.869233e-16*0.3048)
    
    @xflow.setter
    def xflow(self,value):
        """Setter for the rock transmissibility in x-direction."""
        self._xflow = (self.rrock._xperm*self._xarea)/(self._xdelta)

    @property
    def yflow(self):
        """Getter for the rock transmissibility in y-direction."""
        return self._yflow/(9.869233e-16*0.3048)

    @yflow.setter
    def yflow(self,value):
        """Setter for the rock transmissibility in y-direction."""
        self._yflow = (self.rrock._yperm*self._yarea)/(self._ydelta)

    @property
    def zflow(self):
        """Getter for the rock transmissibility in z-direction."""
        return self._zflow/(9.869233e-16*0.3048)

    @zflow.setter
    def zflow(self,value):
        """Setter for the rock transmissibility in z-direction."""
        self._zflow = (self.rrock._zperm*self._zarea)/(self._zdelta)

    @property
    def fluid(self):
        return self._fluid
    
    @fluid.setter
    def fluid(self,value):
        """Sets fluid potential and mobility into the cell."""
        if value is None:
            return

        self._fluid = value
        
        # Reset fluid properties to indicate they need recalculation
        self.hhead = None # Hydrostatic head placeholder
        self.power = None # Fluid Potential (hydrostatic head + fluid pressure)

    @property
    def hhead(self):
        """Getter for the hydrostatic head of grids in psi."""
        return self._hhead/6894.76

    @hhead.setter
    def hhead(self,value):
        """Setter for the fluid's hydrostatic head in Pa."""
        self._hhead = self.fluid._grad*self._depths

    @property
    def power(self):
        """Getter for the fluid's phase potential at grids in psi."""
        return self._power/6894.76
    
    @power.setter
    def power(self,value):
        """Setter for the fluid's phase potential at grids in Pa."""
        self._power = self.fluid._grad*self._depths if self.fluid._press is None else self.fluid._press+self.fluid._grad*self._depths

    @property
    def tcomp(self):
        """Getter for the total compressibility."""
        return None if self._tcomp is None else self._tcomp*6894.76

    @tcomp.setter
    def tcomp(self,value):
        """Setter for the total compressibility."""
        if value is None:
            try:
                self._tcomp = self.rrock._comp+self.fluid._comp
            except Exception as e:
                logging.warning(f"Missing attribute when calculating total compressibility: {e}")
        else:
            self._tcomp = value/6894.76

if __name__ == "__main__":

    pass