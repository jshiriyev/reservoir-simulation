import numpy

from respy.solver._rcube import RecCube

from respy.solver._vector import Vector

class Block(RecCube):

    def __init__(self,grid,**kwargs):
        """Initialization of block (cell) calculation class."""

        super().__init__(*grid(),**kwargs)

    def __call__(self,rrock,fluid,wconds=None,bconds=None,tdelta=1.):
        """Returns pressure updated Vector instance."""

        self.set_rrock(rrock) # permeability & transmissibility
        self.set_fluid(fluid) # potential (power) & mobility

        # accumulation & inter-block transmissibility
        tvect = self.get_tvect(tdelta)

        # well block productivity
        wvect = self.get_wvect(wconds)

        # boundary block transmissibility
        bvect = self.get_bvect(bconds)

        return Vector(*tvect,wvect,bvect)

    def set_rrock(self,rrock):
        """Sets reservoir rock permeability and transmissibility into the cell."""

        self.xperm = rrock._xperm
        self.yperm = rrock._yperm
        self.zperm = rrock._zperm

        self.xflow = self.get_rock_xflow()
        self.yflow = self.get_rock_yflow()
        self.zflow = self.get_rock_zflow()

        self.eporo = rrock._poro
        self.rcomp = rrock._comp

    def set_fluid(self,fluid):
        """Sets fluid potential and mobility into the cell."""

        self.fcomp = fluid._comp

        self.fhead = self.get_fluid_hhead(fluid)
        self.power = self.get_fluid_power(fluid)
        self.mobil = self.get_fluid_mobil(fluid)

    def get_tvect(self,tdelta):
        """Returns A.ct, Tx, Ty, and Tz in the form of flat arrays."""

        acvect = self.get_block_accum(tdelta)

        xneg,xpos = self.xneg,self.xpos
        yneg,ypos = self.yneg,self.ypos
        zneg,zpos = self.zneg,self.zpos

        xrmean = self.get_harmonic_mean(xneg.xflow,xpos.xflow)
        yrmean = self.get_harmonic_mean(yneg.yflow,ypos.yflow)
        zrmean = self.get_harmonic_mean(zneg.zflow,zpos.zflow)

        xfmean = self.get_upwinding_mean(
            xneg.mobil,xpos.mobil,xneg.power,xpos.power)

        yfmean = self.get_upwinding_mean(
            yneg.mobil,ypos.mobil,yneg.power,ypos.power)

        zfmean = self.get_upwinding_mean(
            zneg.mobil,zpos.mobil,zneg.power,zpos.power)

        return (acvect,xrmean*xfmean,yrmean*yfmean,zrmean*zfmean,self.fhead)

    def get_wvect(self,wconds):
        """Returns productivity for all active wells."""
        wconds = () if wconds is None else wconds
        return [self.get_wprod(wcond) for wcond in wconds]

    def get_bvect(self,bconds):
        """Returns transmissibility for all active boundaries."""
        bconds = () if bconds is None else bconds
        return [self.get_bprod(bcond) for bcond in bconds]

    def get_rock_xflow(self):
        """Returns rock transmissibility in x-direction."""
        return (self.xperm*self.xarea)/(self.xedge)

    def get_rock_yflow(self):
        """Returns rock transmissibility in y-direction."""
        return (self.yperm*self.yarea)/(self.yedge)

    def get_rock_zflow(self):
        """Returns rock transmissibility in z-direction."""
        return (self.zperm*self.zarea)/(self.zedge)

    def get_fluid_hhead(self,fluid):
        """Returns fluid hydrostatic head."""
        if "depth" not in self.prop:
            return 0

        return fluid._grad*self.depth

    def get_fluid_power(self,fluid):
        """Returns fluid phase potential."""
        return self.fhead+fluid._press

    def get_fluid_mobil(self,fluid):
        """Returns fluid phase mobility."""
        return (fluid._rperm)/(fluid._visc*fluid._fvf)

    def get_block_accum(self,tdelta):
        """Returns accumulation multiplied compressibility (A.ct)."""
        return self.volume*self.eporo*self.ccomp/tdelta

    def get_wprod(self,wcond):
        """Returns well transmissibility values for the given well condition."""

        well = self[wcond.block]

        if wcond.axis == "x":
            dhk = well.xedge*numpy.sqrt(well.yperm*well.zperm)
            req = self.get_equiv_radius(well.yperm,well.zperm,well.yedge,well.zedge)
        elif wcond.axis == "y":
            dhk = well.yedge*numpy.sqrt(well.zperm*well.xperm)
            req = self.get_equiv_radius(well.zperm,well.xperm,well.zedge,well.xedge)
        elif wcond.axis == "z":
            dhk = well.zedge*numpy.sqrt(well.xperm*well.yperm)
            req = self.get_equiv_radius(well.xperm,well.yperm,well.xedge,well.yedge)

        wcond._prod = (2*numpy.pi*dhk)/(numpy.log(req/wcond.radius)+wcond.skin)
        wcond._prod *= well.mobil
        
        return wcond

    def get_bprod(self,bcond):
        """Returns exterior block transmissibility values for the
        given boundary condition."""

        face = getattr(self,bcond.face)
        flow = getattr(face,f"{bcond.axis}flow")

        bcond._prod = flow*face.mobil

        return bcond

    @property
    def ccomp(self):
        """Returns total compressibility."""
        if "tcomp" in self.prop:
            return self.tcomp
        return self.rcomp+self.fcomp

    @staticmethod
    def get_weighted_mean(term1,term2):
        """Returns weighted mean of two terms, term1 and term 2."""
        return (term1+term2)/2

    @staticmethod
    def get_harmonic_mean(term1,term2):
        """Returns harmonic mean of two terms, term1 and term 2."""
        return (2*term1*term2)/(term1+term2)

    @staticmethod
    def get_geometric_mean(term1,term2):
        """Returns geometric mean of two terms, term1 and term 2."""
        return numpy.sqrt(term1*term2)

    @staticmethod
    def get_upwinding_mean(term1,term2,ppot1,ppot2):
        """Returns upwinding mean of two terms, term1 and term 2
        for the given phase potential values, ppot1 and ppot2."""
        term = term1.copy()

        if term.shape[0]==ppot1.size:
            term[ppot1<ppot2] = term2

        return term

    @staticmethod
    def get_equiv_radius(perm1,perm2,edge1,edge2):
        """Returns equivalent radius of grids that contains well."""
        sqrt21 = numpy.power(perm2/perm1,1/2)*numpy.power(edge1,2)
        sqrt12 = numpy.power(perm1/perm2,1/2)*numpy.power(edge2,2)

        quar21 = numpy.power(perm2/perm1,1/4)
        quar12 = numpy.power(perm1/perm2,1/4)

        return 0.28*numpy.sqrt(sqrt21+sqrt12)/(quar21+quar12)

if __name__ == "__main__":

    pass