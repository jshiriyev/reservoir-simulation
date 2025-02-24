import numpy

from ._cuboid import Cuboid
from ._vector import Vector

class Block(Cuboid):

    def __init__(self,grids,rrock=None,fluid=None,tcomp=None):
        """Initialization of block (cell) calculation class."""
        super().__init__(grids,rrock,fluid,tcomp)

    def vector(self,tstep=1.,wells=None,edges=None):
        """Returns vector to be used in the building."""
        avect = self.accumulation(tstep) # accumulation vector
        tvect = self.transmissibility() # inter-block transmissibility vector
        wvect = self.boreholes(wells) # wells' block productivity vector list
        evect = self.boundaries(edges) # edges' block productivity vector list

        return Vector(avect,*tvect,wvect,evect)

    def accumulation(self,tstep:float):
        """Returns accumulation multiplied compressibility (A.ct)."""
        return (self.volume*self.rrock.poro*self.tcomp)/tstep

    def transmissibility(self):
        """Returns Tx, Ty, and Tz in the form of flat arrays."""
        tx = self.flux(self.xflow,self.xnegbool,self.xposbool)
        ty = self.flux(self.yflow,self.ynegbool,self.yposbool)
        tz = self.flux(self.zflow,self.znegbool,self.zposbool)

        return tx,ty,tz

    def flux(self,flow,negbool,posbool):
        """Helper function to compute mean harmonic and upwinded mobility product."""
        rmean = self.mean_harmonic(flow[negbool],flow[posbool])
        fmean = self.upwinding(
            self.mobil[negbool],self.mobil[posbool],self.power[negbool],self.power[posbool])

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
        return [self.well_productivity(well) for well in wells]

    def well_productivity(self,well):
        """Returns well transmissibility values for the given well condition."""

        if well.axis=="x":
            k1 = self.rrock.yperm[list(well.block)]
            k2 = self.rrock.zperm[list(well.block)]
            w1 = self.grids.ydelta[list(well.block)]
            w2 = self.grids.zdelta[list(well.block)]
            w3 = self.grids.xdelta[list(well.block)]
        elif well.axis=='y':
            k1 = self.rrock.xperm[list(well.block)]
            k2 = self.rrock.zperm[list(well.block)]
            w1 = self.grids.xdelta[list(well.block)]
            w2 = self.grids.zdelta[list(well.block)]
            w3 = self.grids.ydelta[list(well.block)]
        elif well.axis=='z':
            k1 = self.rrock.xperm[list(well.block)]
            k2 = self.rrock.yperm[list(well.block)]
            w1 = self.grids.xdelta[list(well.block)]
            w2 = self.grids.ydelta[list(well.block)]
            w3 = self.grids.zdelta[list(well.block)]

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