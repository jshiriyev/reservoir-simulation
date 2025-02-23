import numpy

from ._vector import Vector

class Block():

    def __init__(self,grids):
        """Initialization of block (cell) calculation class."""
        self.grids = grids

    def __getattr__(self,key):
        """Delegates attribute access to the underlying grid object."""
        return getattr(self.grids,key)

    def __call__(self,rrock,fluid,wconds=None,bconds=None,tstep=1.):
        """Returns pressure updated Vector instance."""

        self.rrock = rrock # permeability and transmissibility
        self.fluid = fluid # hydrostatic head, potential and mobility

        # accumulation & inter-block transmissibility
        tvect = self.cell_vector(tstep)

        # well block productivity
        wvect = self.well_vector(wconds)

        # boundary block transmissibility
        bvect = self.edge_vector(bconds)

        return Vector(*tvect,wvect,bvect)

    @rrock.setter
    def rrock(self,rrock):
        """Sets reservoir rock permeability and transmissibility into the cell."""
        self.xflow = rrock
        self.yflow = rrock
        self.zflow = rrock

    @fluid.setter
    def fluid(self,fluid):
        """Sets fluid potential and mobility into the cell."""
        self.hhead = fluid
        self.power = fluid
        self.mobil = fluid

    @xflow.setter
    def xflow(self,rrock):
        """Setter for the rock transmissibility in x-direction."""
        self._xflow = (rrock._xperm*self._xarea)/(self._xdelta)

    @yflow.setter
    def yflow(self,rrock):
        """Setter for the rock transmissibility in y-direction."""
        self._yflow = (rrock._yperm*self._yarea)/(self._ydelta)

    @zflow.setter
    def zflow(self,rrock):
        """Setter for the rock transmissibility in z-direction."""
        self._zflow = (rrock._zperm*self._zarea)/(self._zdelta)

    @hhead.setter
    def hhead(self,fluid):
        """Setter for the fluid's hydrostatic head."""
        self._hhead = fluid._grad*self._depth

    @power.setter
    def power(self,fluid):
        """Setter for the fluid's phase potential."""
        self._power = fluid._press+self._hhead

    @mobil.setter
    def mobil(self,fluid):
        """Setter for the fluid's phase mobility."""
        self._mobil = (fluid._rperm)/(fluid._visc*fluid._fvf)

    def cell_vector(self,tstep):
        """Returns A.ct, Tx, Ty, and Tz in the form of flat arrays."""

        acvect = self.accumulation(tstep)

        xneg,xpos = self.xneg,self.xpos
        yneg,ypos = self.yneg,self.ypos
        zneg,zpos = self.zneg,self.zpos

        xrmean = self.mean_harmonic(xneg.xflow,xpos.xflow)
        yrmean = self.mean_harmonic(yneg.yflow,ypos.yflow)
        zrmean = self.mean_harmonic(zneg.zflow,zpos.zflow)

        xfmean = self.upwinding(
            xneg.mobil,xpos.mobil,xneg.power,xpos.power)

        yfmean = self.upwinding(
            yneg.mobil,ypos.mobil,yneg.power,ypos.power)

        zfmean = self.upwinding(
            zneg.mobil,zpos.mobil,zneg.power,zpos.power)

        return (acvect,xrmean*xfmean,yrmean*yfmean,zrmean*zfmean,self.hhead)

    def accumulation(self,tstep:float):
        """Returns accumulation multiplied compressibility (A.ct)."""
        return (self.volume*self.rrock.poro)/(tstep*86400)

    @property
    def compressibility(self):
        """Returns total compressibility."""
        if "tcomp" in self.prop:
            return self.tcomp
        return self.rcomp+self.fcomp

    def well_vector(self,wconds):
        """Returns productivity for all active wells."""
        wconds = () if wconds is None else wconds
        return [self.edge_productivity(wcond) for wcond in wconds]

    def edge_vector(self,bconds):
        """Returns transmissibility for all active boundaries."""
        bconds = () if bconds is None else bconds
        return [self.edge_productivity(bcond) for bcond in bconds]

    def edge_productivity(self,well):
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

    def edge_productivity(self,edge):
        """Returns exterior block transmissibility values for the
        given boundary condition."""

        face = getattr(self,edge.face)
        flow = getattr(face,f"{edge.axis}flow")

        edge._prod = flow*face.mobil

        return edge

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

    @staticmethod
    def requivalent(perm1,perm2,delta1,delta2):
        """Returns equivalent radius of grids that contains well."""
        sqrt21 = numpy.power(perm2/perm1,1/2)*numpy.power(delta1,2)
        sqrt12 = numpy.power(perm1/perm2,1/2)*numpy.power(delta2,2)

        quar21 = numpy.power(perm2/perm1,1/4)
        quar12 = numpy.power(perm1/perm2,1/4)

        return 0.28*numpy.sqrt(sqrt21+sqrt12)/(quar21+quar12)

if __name__ == "__main__":

    pass