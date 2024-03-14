from scipy.sparse import csr_matrix as csr

from respy.fluid._fluid import Fluid

class Assembly:

    def Tmatrix(self,phase):
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

    def Smatrix(self):
        """
        Sets S matrix filled with pore volume values.
        """
        storage = self.grid._volume*self.grid._poro

        return diags(storage.flatten())

    def Gvector(self,phase):
        """
        Sets G vector filled with dynamic gravity coefficients.
        """
        return phase._density*self._gravity*self._Tmatrix.dot(self.grid._depth)

    def Jmatrix(self,phase):
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

    def Qvector(self,phase):
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