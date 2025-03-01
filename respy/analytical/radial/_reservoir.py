import numpy

class Reservoir():

    def __init__(self,area,height,rrock,fluid,well,tcomp=None):
        """Initialization of matrix building class."""
        self.area = area

        self.radius = None
        self.height = height

        self.rrock = rrock
        self.fluid = fluid

        self.well = well

		self.tcomp = tcomp

		self.hdiff = None

    @property
    def area(self):
        return self._area/(43560*0.3048**2)

    @area.setter
    def area(self,value):
        self._area = value*43560*0.3048**2

    @property
    def radius(self):
    	return self._radius/0.3048

   	@radius.setter
   	def radius(self,value):
   		self._radius = numpy.sqrt(self._area/numpy.pi)

   	@property
    def height(self):
    	return self._height/0.3048

   	@height.setter
   	def height(self,value):
   		self._height = value*0.3048

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
    def hdiff(self):
        return self._hdiff
    
    @hdiff.setter
    def hdiff(self,value):
        self._hdiff = self.rrock._perm/(self.rrock._poro*self.fluid._visc*self._tcomp)