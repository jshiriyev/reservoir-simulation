class Tranny():

    def __init__(self,cube,wconds,bconds):

        self.cube = cube

        self.wconds = wconds
        self.bconds = bconds

        for wcond in self.wconds:
            wcond.cube = self.cube[wcond.block]

        for bcond in self.bconds:
            bcond.cube = getattr(getattr(self.cube,bcond.face),"rows")

    def set_trans(self):

        self._xtrans = self._xrrock*self._xfluid

        if self.cube.dims>1:
            self._ytrans = self._yrrock*self._yfluid

        if self.cube.dims>2:
            self._ztrans = self._zrrock*self._zfluid

        self._wtrans = [wrrock*wfluid for wrrock,wfluid in zip(self._wrrock,self._wfluid)]
        self._btrans = [brrock*bfluid for brrock,bfluid in zip(self._brrock,self._bfluid)]

    def set_rrock(self):
        """Self assigns static transmissibility values."""

        self._xrrock = Tranny.get_rrock_axis('x')

        if self.cube.dims>1:
            self._yrrock = Tranny.get_rrock_axis('y')

        if self.cube.dims>2:
            self._zrrock = Tranny.get_rrock_axis('z')

        self._wrrock = [Tranny.get_rrock_well(wcond) for wcond in self.wconds]
        self._brrock = [Tranny.get_rrock_face(bcond) for bcond in self.bconds]

    def set_fluid(self):

        self._xfluid = Tranny.get_fluid_axis('x')

        if self.cube.dims>1:
            self._yfluid = Tranny.get_fluid_axis('y')

        if self.cube.dims>2:
            self._zfluid = Tranny.get_fluid_axis('z')

        self._wfluid = [Tranny.get_fluid_well(wcond) for wcond in self.wconds]
        self._bfluid = [Tranny.get_fluid_face(bcond) for bcond in self.bconds]

    @staticmethod
    def get_rrock_axis(cube,axis):
        """Returns static transmissibility values for the given
        direction and inner faces (interfaces)."""

        edge_neg = getattr(getattr(cube,f"{axis}neg"),f"_{axis}edge")
        edge_pos = getattr(getattr(cube,f"{axis}pos"),f"_{axis}edge")

        area_neg = getattr(getattr(cube,f"{axis}neg"),f"_{axis}area")
        area_pos = getattr(getattr(cube,f"{axis}pos"),f"_{axis}area")
        
        perm_neg = getattr(getattr(cube,f"{axis}neg"),f"_{axis}perm")
        perm_pos = getattr(getattr(cube,f"{axis}pos"),f"_{axis}perm")

        block_neg = Tranny.get_block_tranny(edge_neg,area_neg,perm_neg)
        block_pos = Tranny.get_block_tranny(edge_pos,area_pos,perm_pos)

        return Tranny.get_harmonic_mean(block_neg,block_pos)

    @staticmethod
    def get_rrock_well(cond):
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
    def get_rrock_face(cond):
        """Returns static transmissibility values for the given
        surface on the exterior boundary."""
        
        edge = getattr(cond.block,f"_{cond.face[0]}edge")
        area = getattr(cond.block,f"_{cond.face[0]}area")
        perm = getattr(cond.block,f"_{cond.face[0]}perm")

        return Tranny.get_block_tranny(edge,area,perm)

    @staticmethod
    def get_block_tranny(edge,area,perm):
        return (perm*area)/edge

    @staticmethod
    def get_equiv_radius(edge1,perm1,edge2,perm2):

        sqrt21 = numpy.power(perm2/perm1,1/2)*numpy.power(edge1,2)
        sqrt12 = numpy.power(perm1/perm2,1/2)*numpy.power(edge2,2)

        quar21 = numpy.power(perm2/perm1,1/4)
        quar12 = numpy.power(perm1/perm2,1/4)

        return 0.28*numpy.sqrt(sqrt21+sqrt12)/(quar21+quar12)

    @staticmethod
    def get_harmonic_mean(perm1,perm2):
        return (2*perm1*perm2)/(perm1+perm2)