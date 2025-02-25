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
        evect = self.boundaries(edges) # edges' block productivity vector list

        return Vector(avect,*tvect,wvect,evect)

    def storage(self,tstep:float):
        """Returns accumulation multiplied compressibility (A.ct)."""
        return (self.volume*self.rrock.poro*self.tcomp)/tstep

    def fluxes(self):
        """Returns Tx, Ty, and Tz in the form of flat arrays."""
        tx = self.flux(self.xflow,self._xneg,self._xpos)
        ty = self.flux(self.yflow,self._yneg,self._ypos)
        tz = self.flux(self.zflow,self._zneg,self._zpos)

        return tx,ty,tz

    def flux(self,flow,bneg,bpos):
        """Helper function to compute mean harmonic and upwinded mobility product."""
        rmean = self.mean_harmonic(flow[bneg],flow[bpos])
        fmean = self.upwinding(
            self.mobil[bneg],self.mobil[bpos],self.power[bneg],self.power[bpos])

        return rmean*fmean

    @staticmethod
    def mean_weighted(term1,term2):
        """Returns weighted mean of two terms, term1 and term 2."""
        return (term1+term2)/2

    @staticmethod
    def mean_harmonic(term1,term2):
        """Returns harmonic mean of two terms, term1 and term 2."""
        return (2*term1*term2)/(term1+term2)

    @staticmethod
    def mean_geometric(term1,term2):
        """Returns geometric mean of two terms, term1 and term 2."""
        return numpy.sqrt(term1*term2)

    @staticmethod
    def upwinding(term1,term2,ppot1,ppot2):
        """Computes the upwinded mean of two terms based on phase potential (ppot) values."""
        return numpy.where(ppot1<ppot2,term2,term1)

    def boreholes(self,wells):
        """Returns productivity for all active wells."""
        wells = () if wells is None else wells
        return [self.borehole(w) for w in wells]

    def borehole(self,well):
        """Returns well transmissibility values for the given well condition."""

        xperm = self.rrock.xperm[list(well.block)]
        yperm = self.rrock.yperm[list(well.block)]
        zperm = self.rrock.zperm[list(well.block)]

        xsize = self.grids.xdelta[list(well.block)]
        ysize = self.grids.ydelta[list(well.block)]
        zsize = self.grids.zdelta[list(well.block)]

        if well.axis=="x":
            k1,k2,k3 = yperm,zperm,xperm
            w1,w2,w3 = ysize,zsize,xsize
        elif well.axis=='y':
            k1,k2,k3 = zperm,xperm,yperm
            w1,w2,w3 = zsize,xsize,ysize
        elif well.axis=='z':
            k1,k2,k3 = xperm,yperm,zperm
            w1,w2,w3 = xsize,ysize,zsize

        well.prod = (2*numpy.pi*w3*numpy.sqrt(k1*k2))/(
            numpy.log(self.requivalent(k1,k2,w1,w2)/well.radius)+well.skin
            )*self.mobil[list(well.block)]
        
        return well

    @staticmethod
    def requivalent(perm1,perm2,delta1,delta2):
        """Returns equivalent radius of grids that contains well."""
        sqrt21 = numpy.power(perm2/perm1,1/2)*numpy.power(delta1,2)
        sqrt12 = numpy.power(perm1/perm2,1/2)*numpy.power(delta2,2)

        quar21 = numpy.power(perm2/perm1,1/4)
        quar12 = numpy.power(perm1/perm2,1/4)

        return 0.28*numpy.sqrt(sqrt21+sqrt12)/(quar21+quar12)

    def boundaries(self,edges):
        """Returns transmissibility for all active boundaries."""
        edges = () if edges is None else edges
        return [self.edge_productivity(edge) for edge in edges]

    def edge_productivity(self,edge):
        """Returns exterior block transmissibility values for the
        given boundary condition."""

        face = getattr(self,edge.face)
        flow = getattr(face,f"{edge.axis}flow")

        edge._prod = flow*face.mobil

        return edge

if __name__ == "__main__":

    pass