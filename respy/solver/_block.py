import sys

if __name__ == "__main__":
    # sys.path.append(r'C:\Users\javid.shiriyev\Documents\respy')
    sys.path.append(r'C:\Users\3876yl\Documents\respy')

import numpy

from respy.solver._rcube import RecCube

from respy.solver._vector import Vector

class Block(RecCube):

    def __init__(self,grid,wconds,bconds,depth=None):
        """Initialization of block calculation class."""

        # grid interface
        object.__setattr__(self,"grid",grid)

        # initializing the parent class
        super().__init__(*grid(),depth=depth)

        # well and boundary condition interfaces
        object.__setattr__(self,"wconds",wconds)
        object.__setattr__(self,"bconds",bconds)

    def __call__(self,rrock,fluid):
        """Returns pressure updated Vector instance."""

        # 1. Rock calculations
        rcell = self.get_rcell(rrock) # rock storage and flow capacity

        # 2. Fluid calculations
        fcell = self.get_fcell(fluid) # mobility and potential

        # 3. Inter-block transmissibility calculations
        rmean = self.get_rmean(rcell) # inter-block mean rock transmissibility
        fmean = self.get_fmean(fcell) # inter-block mean fluid mobility
        tmean = self.get_tmean(rmean,fmean) # inter-block mean transmissibility

        # 4. Well block calculations
        wvect = [self.get_wvect(rrock,fcell,wcond) for wcond in self.wconds]

        # 5. Boundary block calculations
        bvect = [self.get_bvect(rcell,fcell,bcond) for bcond in self.bconds]

        return Vector(rcell.pvol,*tmean,wvect,bvect)

    def get_rcell(self,rrock):
        """Returns rock cell transmissibility."""

        pvol = self.get_rpvol(rrock)
        xdir = self.get_rxdir(rrock)
        ydir = self.get_rydir(rrock)
        zdir = self.get_rzdir(rrock)

        return block(pvol=pvol,xdir=xdir,ydir=ydir,zdir=zdir)

    def get_rpvol(self,rrock):
        """Returns rock pore volume."""
        return self.volume*rrock._poro

    def get_rxdir(self,rrock):
        """Returns rock cell transmissibility in x-direction."""
        return (rrock._xperm*self.xarea)/(self.xedge)

    def get_rydir(self,rrock):
        """Returns rock cell transmissibility in y-direction."""
        return (rrock._yperm*self.yarea)/(self.yedge)

    def get_rzdir(self,rrock):
        """Returns rock cell transmissibility in z-direction."""
        return (rrock._zperm*self.zarea)/(self.zedge)

    def get_rmean(self,rcell):
        """Returns rock inter-block transmissibility term"""

        xface = self.get_harmonic_mean(
            rcell.xdir[self.grid.xneg],
            rcell.xdir[self.grid.xpos])

        yface = self.get_harmonic_mean(
            rcell.ydir[self.grid.yneg],
            rcell.ydir[self.grid.ypos])

        zface = self.get_harmonic_mean(
            rcell.zdir[self.grid.zneg],
            rcell.zdir[self.grid.zpos])

        return block(xface=xface,yface=yface,zface=zface)

    def get_fcell(self,fluid,rperm=1):
        """Returns fluid cell potential and mobility terms."""

        # hydrostatic head calculations
        hhead = 0 if self.depth is None else fluid._rho*9.807*self.depth

        # fluid phase potential calculations
        power = fluid._press+hhead

        # fluid phase mobility initialization
        mobil = numpy.empty(self.shape)

        # fluid phase mobility calculations
        mobil[:] = (rperm)/(fluid._visc*fluid._fvf)

        return block(power=power,mobil=mobil)

    def get_fmean(self,fcell):
        """Returns fluid mobility term after upwinding"""

        xface = self.get_upwinding_mean(
            fcell.mobil[self.grid.xneg],
            fcell.power[self.grid.xneg],
            fcell.mobil[self.grid.xpos],
            fcell.power[self.grid.xpos])

        yface = self.get_upwinding_mean(
            fcell.mobil[self.grid.yneg],
            fcell.power[self.grid.yneg],
            fcell.mobil[self.grid.ypos],
            fcell.power[self.grid.ypos])

        zface = self.get_upwinding_mean(
            fcell.mobil[self.grid.zneg],
            fcell.power[self.grid.zneg],
            fcell.mobil[self.grid.zpos],
            fcell.power[self.grid.zpos])

        return block(xface=xface,yface=yface,zface=zface)

    def get_tmean(self,rmean,fmean):

        xvect = rmean.xface*fmean.xface
        yvect = rmean.yface*fmean.yface
        zvect = rmean.zface*fmean.zface

        return (xvect,yvect,zvect)

    def get_wvect(self,rrock,fcell,wcond):
        """Returns well transmissibility values for the given well condition."""

        dx = self.xedge[wcond.block]
        dy = self.yedge[wcond.block]
        dz = self.zedge[wcond.block]

        kx = rrock._xperm[wcond.block]
        ky = rrock._yperm[wcond.block]
        kz = rrock._zperm[wcond.block]

        #fluid mobility values for the well blocks
        fm = fcell.mobil[wcond.block]

        if wcond.axis == "x":
            dhk = dx*numpy.sqrt(ky*kz)
            req = self.get_equiv_radius(dy,ky,dz,kz)
        elif wcond.axis == "y":
            dhk = dy*numpy.sqrt(kz*kx)
            req = self.get_equiv_radius(dz,kz,dx,kx)
        elif wcond.axis == "z":
            dhk = dz*numpy.sqrt(kx*ky)
            req = self.get_equiv_radius(dx,kx,dy,ky)

        return (2*numpy.pi*dhk)*fm/(numpy.log(req/wcond.radius)+wcond.skin)

    def get_bvect(self,rcell,fcell,bcond):
        """Returns exterior block transmissibility values for the
        given boundary condition."""

        rock = getattr(rcell,f"{bcond.face[0]}dir")
        rows = getattr(self.grid,bcond.face)

        return rock[rows]*fcell.mobil[rows]

    """Static methods:"""

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
    def get_upwinding_mean(term1,ppot1,term2,ppot2):
        """Returns upwinding mean of two terms, term1 and term 2
        for the given phase potential values, ppot1 and ppot2."""
        term = term1.copy()

        term[ppot1<ppot2] = term2

        return term

    @staticmethod
    def get_equiv_radius(edge1,perm1,edge2,perm2):
        """Returns equivalent radius of grids that contains well."""
        sqrt21 = numpy.power(perm2/perm1,1/2)*numpy.power(edge1,2)
        sqrt12 = numpy.power(perm1/perm2,1/2)*numpy.power(edge2,2)

        quar21 = numpy.power(perm2/perm1,1/4)
        quar12 = numpy.power(perm1/perm2,1/4)

        return 0.28*numpy.sqrt(sqrt21+sqrt12)/(quar21+quar12)

class block():

    def __init__(**kwargs):

        for key,value in kwargs.items():
            setattr(self,key,value)