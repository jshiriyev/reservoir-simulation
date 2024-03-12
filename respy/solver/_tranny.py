class Tranny:

    @staticmethod
    def set_static(self):
        """Self assigns static transmissibility values."""

        self._staticx = Static.get_axis('x')

        if self.grid.dims>1:
            self._staticy = Static.get_axis('y')

        if self.grid.dims>2:
            self._staticz = Static.get_axis('z')

        for bcond in self.bconds:
            bcond.block = getattr(self.grid,bcond.face)

        self._staticw = [Static.get_well(wcond) for wcond in self.wconds]
        self._staticb = [Static.get_face(bcond) for bcond in self.bconds]

    @staticmethod
    def get_axis(cube,axis):
        """Returns static transmissibility values for the given
        direction and inner faces (interfaces)."""

        edge_neg = getattr(getattr(cube,f"{axis}neg"),f"_{axis}edge")
        edge_pos = getattr(getattr(cube,f"{axis}pos"),f"_{axis}edge")

        area_neg = getattr(getattr(cube,f"{axis}neg"),f"_{axis}area")
        area_pos = getattr(getattr(cube,f"{axis}pos"),f"_{axis}area")
        
        perm_neg = getattr(getattr(cube,f"{axis}neg"),f"_{axis}perm")
        perm_pos = getattr(getattr(cube,f"{axis}pos"),f"_{axis}perm")

        block_neg = Static.get_block_tranny(edge_neg,area_neg,perm_neg)
        block_pos = Static.get_block_tranny(edge_pos,area_pos,perm_pos)

        return Static.get_harmonic_mean(block_neg,block_pos)

    @staticmethod
    def get_well(cube,cond):
        """Returns static transmissibility values for the given
        vertical well."""

        dx = cube._xedge[cond.block]
        dy = cube._yedge[cond.block]
        dz = cube._zedge[cond.block]

        kx = cube._xperm[cond.block]
        ky = cube._yperm[cond.block]
        kz = cube._zperm[cond.block]

        if cond.axis == "x":
            dkh = dx*numpy.sqrt(ky*kz)
            req = Static.get_equiv_radius(dy,ky,dz,kz)
        elif cond.axis == "y":
            dkh = dy*numpy.sqrt(kz*kx)
            req = Static.get_equiv_radius(dz,kz,dx,kx)
        elif cond.axis == "z":
            dkh = dz*numpy.sqrt(kx*ky)
            req = Static.get_equiv_radius(dx,kx,dy,ky)

        return (2*numpy.pi*dkh)/(numpy.log(req/cond.radius)+cond.skin)

    @staticmethod
    def get_face(cond):
        """Returns static transmissibility values for the given
        surface on the exterior boundary."""
        
        edge = getattr(cond.block,f"_{cond.face[0]}edge")
        area = getattr(cond.block,f"_{cond.face[0]}area")
        perm = getattr(cond.block,f"_{cond.face[0]}perm")

        return Static.get_block_tranny(edge,area,perm)

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