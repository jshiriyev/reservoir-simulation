class Block():

    def __init__(self,cube,wconds,bconds):

        self.cube = cube

        self.wconds = wconds
        self.bconds = bconds

    @property
    def xrock(self):
        return self.get_block_xrock(self.cube)

    @property
    def yrock(self):
        return self.get_block_yrock(self.cube)

    @property
    def zrock(self):
        return self.get_block_zrock(self.cube)

    @property
    def fluid(self):
        return get_block_fluid(self.cube)
    
    def update(self):

        self.cube.xtrans = self.fluid*self.xrock
        self.cube.ytrans = self.fluid*self.yrock
        self.cube.ztrans = self.fluid*self.zrock

    def set_trans(self):

        xtrans = Tranny.get_harmonic_mean(self.cube.xneg.xtrans,self.cube.xpos.xtrans)
        ytrans = Tranny.get_harmonic_mean(self.cube.yneg.ytrans,self.cube.ypos.ytrans)
        ztrans = Tranny.get_harmonic_mean(self.cube.zneg.ztrans,self.cube.zpos.ztrans)

        wtrans = [wrrock*wfluid for wrrock,wfluid in zip(self._wrrock,self._wfluid)]
        btrans = [brrock*bfluid for brrock,bfluid in zip(self._brrock,self._bfluid)]

        return Vectors(xtrans,ytrans,ztrans,wtrans,btrans)

    @staticmethod
    def get_block_xrock(cube):
        return (cube._xperm*cube._xarea)/(cube._xedge)

    @staticmethod
    def get_block_yrock(cube):
        return (cube._yperm*cube._yarea)/(cube._yedge)

    @staticmethod
    def get_block_zrock(cube):
        return (cube._zperm*cube._zarea)/(cube._zedge)

    @staticmethod
    def get_block_fluid(cube):
        return (cube._rel0)/(cube._visc*cube._fvf)

    @staticmethod
    def get_block_wrock(cond):
        """Returns static transmissibility values for the given
        vertical well."""

        dx = cond.cube._xedge
        dy = cond.cube._yedge
        dz = cond.cube._zedge

        kx = cond.cube._xperm
        ky = cond.cube._yperm
        kz = cond.cube._zperm

        if cond.axis == "x":
            dhk = dx*numpy.sqrt(ky*kz)
            req = Tranny.get_equiv_radius(dy,ky,dz,kz)
        elif cond.axis == "y":
            dhk = dy*numpy.sqrt(kz*kx)
            req = Tranny.get_equiv_radius(dz,kz,dx,kx)
        elif cond.axis == "z":
            dhk = dz*numpy.sqrt(kx*ky)
            req = Tranny.get_equiv_radius(dx,kx,dy,ky)

        return (2*numpy.pi*dhk)/(numpy.log(req/cond.radius)+cond.skin)

    @staticmethod
    def get_equiv_radius(edge1,perm1,edge2,perm2):

        sqrt21 = numpy.power(perm2/perm1,1/2)*numpy.power(edge1,2)
        sqrt12 = numpy.power(perm1/perm2,1/2)*numpy.power(edge2,2)

        quar21 = numpy.power(perm2/perm1,1/4)
        quar12 = numpy.power(perm1/perm2,1/4)

        return 0.28*numpy.sqrt(sqrt21+sqrt12)/(quar21+quar12)

    @staticmethod
    def get_block_brock(cond):
        """Returns static transmissibility values for the given
        surface on the exterior boundary."""
        return getattr(self,f"get_block_{cond.face[0]}rock")(cond.block)

    @staticmethod
    def get_weighted_mean(perm1,perm2):
        return (perm1+perm2)/2

    @staticmethod
    def get_harmonic_mean(perm1,perm2):
        return (2*perm1*perm2)/(perm1+perm2)

    @staticmethod
    def get_geometric_mean(perm1,perm2):
        return numpy.sqrt(perm1*perm2)