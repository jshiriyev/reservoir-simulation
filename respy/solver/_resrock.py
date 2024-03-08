from respy.utils._prop import Prop

from respy.solver._resinit import ResInit

class ResRock():

    _gravity = 9.807  # Gravitational acceleration in SI units

    def __init__(self,grid,fluids,rcomp=1e-5,wconds=None,bconds=None,**kwargs):
        """
        grid    : RecCuboid instance, rectangular cuboid grids

        fluids  : fluid properties and methods for calculations, tuple of (gas,oil,water)
        
        rcomp 	: rock compressibility

        wconds  : well condition, tuple of WellCond instance
        bconds  : boundary condition, tuple of BoundCond instance
        """

        self.grid = grid

        self.set_static()

        self.fluids = fluids

        self._rcomp = Prop(rcomp,0.000145038).set()

        self.wconds = () if wconds is None else wconds
        self.bconds = () if bconds is None else bconds

        self.resinit = ResInit(**kwargs)

    def rcomp(self,P):

    	return self._rcomp.get(P)

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

    def get_static_axis(self,axis):
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

    def get_static_well(self,cond):
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

    def get_static_face(self,cond):
        """Returns static transmissibility values for the given
        surface on the exterior boundary."""
        
        dims = getattr(cond.block,f"_{cond.face[0]}dims")
        area = getattr(cond.block,f"_{cond.face[0]}area")
        perm = getattr(cond.block,f"_{cond.face[0]}perm")

        return self.get_block_tranny(dims,area,perm)

    def init(self,pcow,pcog,pcgw):

    	Pw0 = self.resinit.waterpressure(self.grid._depth)

    	Po0 = self.resinit.oilpressure(self.grid._depth)
    	Pg0 = self.resinit.gaspressure(self.grid._depth)

    	Sw0,So0,Sg0 = self.resinit.saturations(pcow,pcog,pcgw)

    	cw0 = Sw0*self.fluids['water']._comp

    	co0 = So0*self.fluid['oil']._comp

    	cg0 = Sg0*self.fluid['gas']._comp

    	ct0 = cform+cw0+co0+cg0

    def __call__(self,press=None,temp=None):

        mats = {}

        for key,fluid in self.fluids.items():

            phase = fluid(press,temp)

            T = self.get_Tmatrix(phase)
            S = self.get_Smatrix()
            G = self.get_Gvector(phase)
            J = self.get_Jmatrix(phase)
            Q = self.get_Qvector(phase)

            mats[key] = Matrix(T,V,G,J,Q)

        return mats

    def get_Tmatrix(self,phase):
        """
        Sets transmissibility matrix with dynamic transmissibility values.
        """

        tmatrix = csr(self.matrix)

        tmatrix = self.tcharge('x',tmatrix,phase)

        if self.grid.flodim>1:
            tmatrix = self.tcharge('y',tmatrix,phase)

        if self.grid.flodim>2:
            tmatrix = self.tcharge('z',tmatrix,phase)

        return tmatrix

    def tcharge(self,axis:str,tmatrix:csr,phase:Fluid):
        """
        Returns updated transmissibility matrix:

        axis    : x, y, or z

        tmatrix : csr_matrix object defined for transmissibility matrix

        The transmissibility matrix is updated for diagonal and offset entries.
        """

        stat = getattr(self,f"_static{axis}")

        vals = stat/phase._viscosity

        dneg = getattr(self.grid,f"{axis}neg")
        dpos = getattr(self.grid,f"{axis}pos")

        tmatrix += csr((vals,(dneg,dneg)),shape=tmatrix.shape)
        tmatrix -= csr((vals,(dneg,dpos)),shape=tmatrix.shape)

        tmatrix += csr((vals,(dpos,dpos)),shape=tmatrix.shape)
        tmatrix -= csr((vals,(dpos,dneg)),shape=tmatrix.shape)

        return tmatrix

    def get_Smatrix(self):
        """
        Sets S matrix filled with pore volume values.
        """
        storage = self.grid._volume*self.grid._poro

        return diags(storage.flatten())

    def get_Gvector(self,phase):
        """
        Sets G vector filled with dynamic gravity coefficients.
        """
        return phase._density*self._gravity*self._Tmatrix.dot(self.grid._depth)

    def get_Jmatrix(self,phase):
        """
        Sets J matrix filled with dynamic constant pressure coefficients.
        """

        jmatrix = csr(self.matrix)

        for index,wcond in enumerate(self.wconds):
            jmatrix += self.jcharge(wcond,self._staticw[index],phase,jmatrix)

        for index,bcond in enumerate(self.bconds):
            jmatrix += self.jcharge(bcond,self._staticb[index],phase,jmatrix)

        return jmatrix

    def jcharge(self,cond,static,phase,jmatrix:csr):

        if cond.sort=="press":

            jvalues = 2*static/phase._viscosity
            jmatrix += csr((jvalues,(cond.block,cond.block)),shape=self.matrix)

        return jmatrix

    def get_Qvector(self,phase):
        """
        Sets Q vector filled with dynamic constant rate and pressure coefficients.
        """

        qvector = csr(self.vector)

        for index,wcond in enumerate(self.wconds):
            qvector += self.qcharge(wcond,self._staticw[index],phase,qvector)

        for index,bcond in enumerate(self.bconds):
            qvector += self.qcharge(bcond,self._staticb[index],phase,qvector)

        return qvector

    def qcharge(self,cond,static,phase,qvector:csr):

        tvalues = static/phase._viscosity

        if cond.sort=="press":
            qvalues = 2*tvalues*cond._cond
        else:
            qvalues = tvalues/tvalues.sum()*cond._cond

        indices = numpy.zeros(cond.block.size)

        qvector += csr((qvalues,(cond.block,indices)),shape=self.vector)

        return qvector

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

    @property
    def smatrix(self):
        """Shape of the matrices of transmissibility calculations"""
        return (self.grid.numtot,self.grid.numtot)

    @property
    def svector(self):
        """Shape of the vectors of transmissibility calculations"""
        return (self.grid.numtot,1)
