import sys

if __name__ == "__main__":
    # sys.path.append(r'C:\Users\javid.shiriyev\Documents\respy')
    sys.path.append(r'C:\Users\3876yl\Documents\respy')

import numpy

from respy.solver._vector import Vector

class Block():

    def __init__(self,cube,rrock,fluid,wconds,bconds):
        """Initialization of block calculation class."""

        # cell geometry
        self.cube   = cube

        # rock and fluid properties
        self.rrock  = rrock
        self.fluid  = fluid

        # well and boundary conditions
        self.wconds = wconds
        self.bconds = bconds
    
    def __call__(self,press=None):
        """Returns pressure updated Vector instance."""

        # updating rock and fluid properties with input pressure values
        rrock = self.rrock(press) if callable(self.rrock) else self.rrock
        fluid = self.fluid(press) if callable(self.fluid) else self.fluid

        # 1. Rock calculations: a) storage, b) flow capacity, and c) face hydraulic mean
        svect = self.get_rpvol(rrock) 
        rpart = self.get_rpart(rrock)
        rpart = self.get_rmean(rpart)

        # 2. Fluid calculations: a) mobility, b) potential, and c) face upwinding
        fpart = self.get_fpart(fluid)
        fppot = self.get_fppot(fluid,press)
        fpart = self.get_fmean(fpart,fppot)

        # 3. Face transmissibility calculations
        xvect = rpart._xface*fpart._xface
        yvect = rpart._yface*fpart._yface
        zvect = rpart._zface*fpart._zface

        # 4. Well block calculations
        wvect = [self.get_wpart(rrock,fpart,wcond) for wcond in self.wconds]

        # 5. Boundary block calculations
        bvect = [self.get_bpart(rpart,fpart,bcond) for bcond in self.bconds]

        return Vector(rrock,fluid,svect,xvect,yvect,zvect,wvect,bvect)

    def get_rpvol(self,rrock):
        """returns rock pore volume"""
        return self.cube._volume*rrock._poro

    def get_rpart(self,rrock):
        """returns rock transmissibility term"""
        class rpart: pass

        rpart._xdir = self.get_rxdir(rrock)
        rpart._ydir = self.get_rydir(rrock)
        rpart._zdir = self.get_rzdir(rrock)

        return rpart

    def get_rxdir(self,rrock):
        """returns rock x-transmissibility term"""
        return (rrock._xperm*self.cube._xarea)/(self.cube._xedge)

    def get_rydir(self,rrock):
        """returns rock y-transmissibility term"""
        return (rrock._yperm*self.cube._yarea)/(self.cube._yedge)

    def get_rzdir(self,rrock):
        """returns rock z-transmissibility term"""
        return (rrock._zperm*self.cube._zarea)/(self.cube._zedge)

    def get_rmean(self,rpart):
        """returns rock face (inter-block) transmissibility term"""

        rpart._xface = self.get_harmonic_mean(
            rpart._xdir[self.cube.xneg.rows],
            rpart._xdir[self.cube.xpos.rows],
            )

        rpart._yface = self.get_harmonic_mean(
            rpart._ydir[self.cube.yneg.rows],
            rpart._ydir[self.cube.ypos.rows],
            )

        rpart._zface = self.get_harmonic_mean(
            rpart._zdir[self.cube.zneg.rows],
            rpart._zdir[self.cube.zpos.rows],
            )

        return rpart

    def get_fpart(self,fluid,rperm=1):
        """returns fluid mobility term"""
        class fpart: pass

        fpart._mobil = numpy.empty(self.cube._volume)

        fpart._mobil[:] = (rperm)/(fluid._visc*fluid._fvf)

        return fpart

    def get_fppot(self,fluid,press):
        """returns fluid phase potential in the cell"""
        return press+fluid._rho*9.807*self.rrock._depth

    def get_fmean(self,fpart,fppot):
        """returns fluid mobility term after upwinding"""
        mobil = fpart._mobil

        fpart._xface = self.get_upwinding_mean(
            mobil[self.cube.xneg.rows],
            mobil[self.cube.xpos.rows],
            fppot[self.cube.xneg.rows],
            fppot[self.cube.xpos.rows],
            )

        fpart._yface = self.get_upwinding_mean(
            mobil[self.cube.yneg.rows],
            mobil[self.cube.ypos.rows],
            fppot[self.cube.yneg.rows],
            fppot[self.cube.ypos.rows],
            )

        fpart._zface = self.get_upwinding_mean(
            mobil[self.cube.zneg.rows],
            mobil[self.cube.zpos.rows],
            fppot[self.cube.zneg.rows],
            fppot[self.cube.zpos.rows],
            )

        return fpart

    def get_wpart(self,rrock,fpart,wcond):
        """Returns well transmissibility values for the given well condition."""

        dx = self.cube._xedge[wcond.block]
        dy = self.cube._yedge[wcond.block]
        dz = self.cube._zedge[wcond.block]

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
            rock = rpart._xdir
        elif bcond.face[0] == 'y':
            rock = rpart._ydir
        elif bcond.face[0] == 'z':
            rock = rpart._zdir

        rows = getattr(self.cube,bcond.face).rows

        return rock[rows]*fpart._mobil[rows]

    @staticmethod
    def get_weighted_mean(perm1,perm2):
        return (perm1+perm2)/2

    @staticmethod
    def get_harmonic_mean(perm1,perm2):
        return (2*perm1*perm2)/(perm1+perm2)

    @staticmethod
    def get_geometric_mean(perm1,perm2):
        return numpy.sqrt(perm1*perm2)

    @staticmethod
    def get_upwinding_mean(term1,term2,phi1,phi2):
        
        term = term1.copy()

        term[phi1<phi2] = term2

        return term

    @staticmethod
    def get_equiv_radius(edge1,perm1,edge2,perm2):

        sqrt21 = numpy.power(perm2/perm1,1/2)*numpy.power(edge1,2)
        sqrt12 = numpy.power(perm1/perm2,1/2)*numpy.power(edge2,2)

        quar21 = numpy.power(perm2/perm1,1/4)
        quar12 = numpy.power(perm1/perm2,1/4)

        return 0.28*numpy.sqrt(sqrt21+sqrt12)/(quar21+quar12)

    @property
    def static(self):
        return (callable(self.rrock),callable(self.fluid))
    