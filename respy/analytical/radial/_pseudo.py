from dataclasses import dataclass

@dataclass(frozen=True)
class Boundary:

    #shape factor, {C_A} value
    factor: float

    # PSS is exact for higher values
    time_accurate: float
    # PSS gives less than 1% error for higher values
    time_error_prone: float = None
    # Use infinite system solution with less than 1 % Error for lesser values
    time_infinite: float = None

class Pseudo():

    gamma = 1.781

    GEOMETRY = {
        "circle": Boundary(31.62,0.1),
        "triangle": Boundary(27.6,0.2),
        "square": Boundary(30.8828,0.1),
        "hexagon": Boundary(31.6,0.1),
        }

    def __init__(self,*args,**kwargs):
        # There can be two slightly compressible fluids where the
        # second one is at irreducible saturation, not mobile

        super().__init__(*args,**kwargs)

    def initialize(self,pressure0,Swirr=0,ctotal=None):

        steady.initialize(self,pressure0,Swirr,ctotal)

    @property
    def tmin(self):
        return self._tmin

    @tmin.setter
    def tmin(self,value):
        """setting mimimum time limit because of the wellbore size"""
        self._tmin = self.exactdimtime*self.PorRock.area/self.diffusivity

    @property
    def times(self):
        return self._times

    @times.setter
    def times(self,values):

        bound = values>=self.tmin

        if numpy.any(~bound):
            raise Warning("Not all times satisfy the early time limits!")

        values = values[bound]

        self._times = values.reshape((1,-1))

    @property
    def vpore(self):
        return self._vpore
    
    @vpore.setter
    def vpore(self,value):
        self._vpore = self.res._area*self.res._height*self.rock._poro

    def solve(self):

        if not hasattr(self,"times"):
            self.set_times()

        rateNumer = self.Well.limits*self.Fluids.viscosity[0]
        rateDenom = self.PorRock.permeability*self.PorRock.thickness

        constRate = (rateNumer)/(2*np.pi*rateDenom)

        inner1 = (4*self.PorRock.area)/(self.gamma*self.shapefactor*self.Well.radii[0]**2)
        inner2 = 2*np.pi*self.diffusivity*self.times/self.PorRock.area

        self.deltap = constRate*(1/2*np.log(inner1)+inner2)

        self.pressure = self.pressure0-self.deltap
