import sys

if __name__ == "__main__":
    # sys.path.append(r'C:\Users\javid.shiriyev\Documents\respy')
    sys.path.append(r'C:\Users\3876yl\Documents\respy')

import numpy

from respy.solver._rcube import RecCube

from respy.solver._vector import Vector

class Block():

    def __init__(self,grid,wconds,bconds):
        """Initialization of block calculation class."""

        # cell geometry
        self.grid   = grid

        # well and boundary conditions
        self.wconds = wconds
        self.bconds = bconds
        
        # rectangular cuboid for block calculations
        self.cube  = self.get_rcube(grid)

    def __call__(self,press,rrock,fluid):
        """Returns pressure updated Vector instance."""

        # 1. Rock calculations: a) storage, b) flow capacity, and c) face hydraulic mean
        rpvol = self.get_rpvol(rrock) # rock pore volume
        rcell = self.get_rcell(rrock) # rock cell transmissibility
        rmean = self.get_rmean(rcell) # rock mean transmissibility

        # 2. Fluid calculations: a) mobility, b) potential, and c) face upwinding
        fcell = self.get_fcell(fluid) # fluid cell mobility
        fppot = self.get_fppot(fluid,press) # fluid cell potential
        fmean = self.get_fmean(fcell,fppot) # fluid mean mobility

        # 3. Face transmissibility calculations
        xvect = rmean._xface*fmean._xface
        yvect = rmean._yface*fmean._yface
        zvect = rmean._zface*fmean._zface

        # 4. Well block calculations
        wvect = [self.get_wpart(rrock,fcell,wcond) for wcond in self.wconds]

        # 5. Boundary block calculations
        bvect = [self.get_bpart(rcell,fcell,bcond) for bcond in self.bconds]

        return Vector(rpvol,xvect,yvect,zvect,wvect,bvect)

    def get_rpvol(self,rrock):
        """returns rock pore volume"""
        return self.cube._volume*rrock._poro

    def get_rcell(self,rrock):
        """returns rock transmissibility term"""

        xdir = self.get_rxdir(rrock)
        ydir = self.get_rydir(rrock)
        zdir = self.get_rzdir(rrock)

        return block(xdir=xdir,ydir=ydir,zdir=zdir)

    def get_rxdir(self,rrock):
        """returns rock x-transmissibility term"""
        return (rrock._xperm*self.cube._xarea)/(self.cube._xedge)

    def get_rydir(self,rrock):
        """returns rock y-transmissibility term"""
        return (rrock._yperm*self.cube._yarea)/(self.cube._yedge)

    def get_rzdir(self,rrock):
        """returns rock z-transmissibility term"""
        return (rrock._zperm*self.cube._zarea)/(self.cube._zedge)

    def get_rmean(self,rbloc):
        """returns rock face (inter-block) transmissibility term"""

        xface = self.get_harmonic_mean(
            rbloc._xdir[self.grid.xneg],
            rbloc._xdir[self.grid.xpos],
            )

        yface = self.get_harmonic_mean(
            rbloc._ydir[self.grid.yneg],
            rbloc._ydir[self.grid.ypos],
            )

        zface = self.get_harmonic_mean(
            rbloc._zdir[self.grid.zneg],
            rbloc._zdir[self.grid.zpos],
            )

        return block(xface=xface,yface=yface,zface=zface)

    def get_fcell(self,fluid,rperm=1):
        """returns fluid mobility term"""

        mobil = numpy.empty(self.cube._volume)

        mobil[:] = (rperm)/(fluid._visc*fluid._fvf)

        return block(mobil=mobil)

    def get_fppot(self,fluid,press):
        """returns fluid phase potential in the cell"""
        return press+fluid._rho*9.807*self.rrock._depth

    def get_fmean(self,fpart,fppot):
        """returns fluid mobility term after upwinding"""
        mobil = fpart._mobil

        fpart._xface = self.get_upwinding_mean(
            mobil[self.grid.xneg],
            mobil[self.grid.xpos],
            fppot[self.grid.xneg],
            fppot[self.grid.xpos],
            )

        fpart._yface = self.get_upwinding_mean(
            mobil[self.grid.yneg],
            mobil[self.grid.ypos],
            fppot[self.grid.yneg],
            fppot[self.grid.ypos],
            )

        fpart._zface = self.get_upwinding_mean(
            mobil[self.grid.zneg],
            mobil[self.grid.zpos],
            fppot[self.grid.zneg],
            fppot[self.grid.zpos],
            )

        return fpart

    def get_wpart(self,rrock,fpart,wcond):
        """Returns well transmissibility values for the given well condition."""

        rcube = self.cube

        dx = rcube.xedge[wcond.block]
        dy = rcube.yedge[wcond.block]
        dz = rcube.zedge[wcond.block]

        kx = rrock._xperm[wcond.block]
        ky = rrock._yperm[wcond.block]
        kz = rrock._zperm[wcond.block]

        #fluid mobility calculations
        fm = fpart._mobil[wcond.block]

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

    def get_bpart(self,rpart,fpart,bcond):
        """Returns exterior block transmissibility values for the
        given boundary condition."""

        if bcond.face[0] == 'x':
            rrock = rpart._xdir
        elif bcond.face[0] == 'y':
            rrock = rpart._ydir
        elif bcond.face[0] == 'z':
            rrock = rpart._zdir

        rows = getattr(self.grid,bcond.face)

        return rrock[rows]*fpart._mobil[rows]

    """Static methods:"""

    @staticmethod
    def get_rcube(grid):
        """Returns an instance of RecCube."""
        return RecCube(*grid())

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
    def get_upwinding_mean(term1,term2,phi1,phi2):
        """Returns upwinding mean of two terms, term1 and term 2
        for the given potential values, phi1 and phi2."""
        term = term1.copy()

        term[phi1<phi2] = term2

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