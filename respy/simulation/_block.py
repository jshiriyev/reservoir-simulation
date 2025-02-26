import numpy

from ._cuboid import Cuboid
from ._vector import Vector

class Block(Cuboid):

    def __init__(self,grids,rrock,fluid,tcomp=None):
        """Initialization of block (cell) calculation class."""
        super().__init__(grids,rrock,fluid,tcomp)

    def vector(self,tstep=1.,wells=None,edges=None):
        """Returns vector to be used in the building."""
        avect = self.storage(tstep) # accumulation vector
        tvect = self.fluxes() # inter-block transmissibility vector
        wvect = self.boreholes(wells) # wells' block productivity vector list
        evect = self.edges(edges) # edges' block productivity vector list

        return Vector(avect,*tvect,wvect,evect)

    def storage(self,_tstep:float):
        """Returns accumulation multiplied compressibility (A.ct).
        Input time-step, _tstep, must be in SI units."""
        acomp = (self._volume*self.rrock._poro*self._tcomp)/(_tstep)
        return acomp

    def fluxes(self):
        """Returns xflux, yflux, and zflux in tuple."""
        fluid = (self._mobil,self._power)

        xvect = Mean.diffuse(self._xflow,*fluid,self._xneg,self._xpos)
        yvect = Mean.diffuse(self._yflow,*fluid,self._yneg,self._ypos)
        zvect = Mean.diffuse(self._zflow,*fluid,self._zneg,self._zpos)

        return xvect,yvect,zvect

    def boreholes(self,wells):
        """Returns productivity for all active wells."""
        prods,wells = [],() if wells is None else wells

        for well in wells:

            xsize = self._xdelta[list(well.block)]
            ysize = self._ydelta[list(well.block)]
            zsize = self._zdelta[list(well.block)]

            xperm = self.rrock._xperm[list(well.block)]
            yperm = self.rrock._yperm[list(well.block)]
            zperm = self.rrock._zperm[list(well.block)]

            if well.axis=="x":
                k1,k2,w1,w2,w3 = yperm,zperm,ysize,zsize,xsize
            elif well.axis=='y':
                k1,k2,w1,w2,w3 = zperm,xperm,zsize,xsize,ysize
            elif well.axis=='z':
                k1,k2,w1,w2,w3 = xperm,yperm,xsize,ysize,zsize

            wprop = (well._radius,well.skin)

            prods.append(Mean.potency(k1,k2,w1,w2,w3,*wprop)*self._mobil[list(well.block)])
        
        return prods

    def edges(self,edges):
        """Returns transmissibility for all active boundaries."""
        prods,edges = [],() if edges is None else edges

        for edge in edges:

            bools = getattr(self,f"_{edge.face}")
            flows = getattr(self,f"_{edge.axis}flow")[bools]

            prods.append(flows*self._mobil[bools])

        return prods

class Mean:

    @staticmethod
    def weighted(term1,term2):
        """Returns weighted mean of two terms, term1 and term 2."""
        return (term1+term2)/2

    @staticmethod
    def harmonic(term1,term2):
        """Returns harmonic mean of two terms, term1 and term 2."""
        return (2*term1*term2)/(term1+term2)

    @staticmethod
    def geometric(term1,term2):
        """Returns geometric mean of two terms, term1 and term 2."""
        return numpy.sqrt(term1*term2)

    @staticmethod
    def upwinded(term1,term2,ppot1,ppot2):
        """Computes the upwinded mean of two terms based on phase potential (ppot) values."""
        return numpy.where(ppot1<ppot2,term2,term1)

    @staticmethod
    def diffuse(flow_,mobil,power,_neg,_pos):
        """Helper function to compute mean harmonic and upwinded mobility product."""
        mean_condunce = Mean.harmonic(flow_[_neg],flow_[_pos])
        mean_mobility = Mean.upwinded(mobil[_neg],mobil[_pos],power[_neg],power[_pos])
        return mean_condunce*mean_mobility

    @staticmethod
    def radius(perm1,perm2,delta1,delta2):
        """Returns equivalent radius of grids that contains well."""
        sqrt21 = numpy.power(perm2/perm1,1/2)*numpy.power(delta1,2)
        sqrt12 = numpy.power(perm1/perm2,1/2)*numpy.power(delta2,2)
        quar21 = numpy.power(perm2/perm1,1/4)
        quar12 = numpy.power(perm1/perm2,1/4)
        return 0.28*numpy.sqrt(sqrt21+sqrt12)/(quar21+quar12)

    @staticmethod
    def potency(perm1,perm2,delta1,delta2,delta3,well_radius,well_skin):
        """Returns well's productivity value in a given grid."""
        meanper = Mean.geometric(perm1,perm2)
        meanrad = Mean.radius(perm1,perm2,delta1,delta2)
        logterm = numpy.log(meanrad/well_radius)
        return (2*numpy.pi*delta3*meanper)/(logterm+well_skin)

if __name__ == "__main__":

    pass