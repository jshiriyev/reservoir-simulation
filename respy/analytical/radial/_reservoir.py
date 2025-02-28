import numpy

class Reservoir():

    def __init__(self,radius,height,rrock,fluid,tcomp=None):
        """Initialization of matrix building class."""
        self.radius = radius
        self.height = height

        self.rrock = rrock
        self.fluid = fluid

		self.tcomp = tcomp

		self.diffusivity = None

    @property
    def radius(self):
    	return self._radius/0.3048

   	@radius.setter
   	def radius(self,value):
   		self._radius = value*0.3048

   	@property
    def height(self):
    	return self._height/0.3048

   	@height.setter
   	def height(self,value):
   		self._radius = value*0.3048

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

   	@property
    def diffusivity(self):
        return self._diffusivity
    
    @diffusivity.setter
    def diffusivity(self,value):
        self._diffusivity = self.rrock._perm/(self.rrock._poro*self.fluid._visc*self._tcomp)
    