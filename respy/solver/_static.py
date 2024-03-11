class Static:

    @staticmethod
    def set_static(self):
        """Self assigns static transmissibility values."""

        self._staticx = self.get_static_axis('x')

        if self.grid.flodim>1:
            self._staticy = self.get_static_axis('y')

        if self.grid.flodim>2:
            self._staticz = self.get_static_axis('z')

        for bcond in self.bconds:
            bcond.block = getattr(self.grid,bcond.face)

        self._staticw = [self.get_static_well(wcond) for wcond in self.wconds]
        self._staticb = [self.get_static_face(bcond) for bcond in self.bconds]

    @staticmethod
    def get_axis(self,axis):
        """Returns static transmissibility values for the given
        direction and inner faces (interfaces)."""

        dims_neg = getattr(getattr(self.grid,f"{axis}neg"),f"_{axis}dims")
        dims_pos = getattr(getattr(self.grid,f"{axis}pos"),f"_{axis}dims")

        area_neg = getattr(getattr(self.grid,f"{axis}neg"),f"_{axis}area")
        area_pos = getattr(getattr(self.grid,f"{axis}pos"),f"_{axis}area")
        
        perm_neg = getattr(getattr(self.grid,f"{axis}neg"),f"_{axis}perm")
        perm_pos = getattr(getattr(self.grid,f"{axis}pos"),f"_{axis}perm")

        block_neg = self.get_block_tranny(dims_neg,area_neg,perm_neg)
        block_pos = self.get_block_tranny(dims_pos,area_pos,perm_pos)

        return self.get_harmonic_mean(block_neg,block_pos)

    @staticmethod
    def get_well(self,cond):
        """Returns static transmissibility values for the given
        vertical well."""

        dx = self.grid._xdims[cond.block]
        dy = self.grid._ydims[cond.block]
        dz = self.grid._zdims[cond.block]

        kx = self.grid._xperm[cond.block]
        ky = self.grid._yperm[cond.block]
        kz = self.grid._zperm[cond.block]

        if cond.axis == "x":
            dkh = dx*numpy.sqrt(ky*kz)
            req = self.get_equiv_radius(dy,ky,dz,kz)
        elif cond.axis == "y":
            dkh = dy*numpy.sqrt(kz*kx)
            req = self.get_equiv_radius(dz,kz,dx,kx)
        elif cond.axis == "z":
            dkh = dz*numpy.sqrt(kx*ky)
            req = self.get_equiv_radius(dx,kx,dy,ky)

        return (2*numpy.pi*dkh)/(numpy.log(req/cond.radius)+cond.skin)

    @staticmethod
    def get_face(self,cond):
        """Returns static transmissibility values for the given
        surface on the exterior boundary."""
        
        dims = getattr(cond.block,f"_{cond.face[0]}dims")
        area = getattr(cond.block,f"_{cond.face[0]}area")
        perm = getattr(cond.block,f"_{cond.face[0]}perm")

        return self.get_block_tranny(dims,area,perm)

    @staticmethod
    def get_block_tranny(dims,area,perm):
        return (perm*area)/dims

    @staticmethod
    def get_equiv_radius(dims1,perm1,dims2,perm2):

        sqrt21 = numpy.power(perm2/perm1,1/2)*numpy.power(dims1,2)
        sqrt12 = numpy.power(perm1/perm2,1/2)*numpy.power(dims2,2)

        quar21 = numpy.power(perm2/perm1,1/4)
        quar12 = numpy.power(perm1/perm2,1/4)

        return 0.28*numpy.sqrt(sqrt21+sqrt12)/(quar21+quar12)

    @staticmethod
    def get_harmonic_mean(perm1,perm2):
        return (2*perm1*perm2)/(perm1+perm2)