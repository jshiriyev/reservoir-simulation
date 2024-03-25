from scipy.sparse import csr_matrix as csr

from scipy.sparse import diags

from respy.solver._vector import Vector
from respy.solver._matrix import Matrix

class Build():

    def __init__(self,cube,wconds,bconds):

        self.cube = cube

        self.wconds = wconds
        self.bconds = bconds

    def __call__(self,vec:Vector):

        C = self.get_cmatrix(vec)
        T = self.get_tmatrix(vec)
        G = self.get_gvector(vec,T)
        J = self.get_jmatrix(vec)
        Q = self.get_qvector(vec)

        return Matrix(C,T,G,J,Q)

    def get_cmatrix(self,vec:Vector):
        """Returns C matrix filled with diagonal values."""
        cmatrix = diags(vec._C,shape=self.matrix)

        return cmatrix

    def get_tmatrix(self,vec:Vector):
        """Returns T matrix filled with diagonal and offset values."""
        tmatrix = csr(self.matrix)

        tmatrix = self.set_tmatrix_interblock(
            tmatrix,vec._X,self.cube.xneg,self.cube.xpos)

        tmatrix = self.set_tmatrix_interblock(
            tmatrix,vec._Y,self.cube.yneg,self.cube.ypos)

        tmatrix = self.set_tmatrix_interblock(
            tmatrix,vec._Z,self.cube.zneg,self.cube.zpos)

        return tmatrix

    def get_gvector(self,vec,tmatrix):
        """Returns G vector filled with gravity coefficients."""
        gvector = vec.fluid._rho*9.807*tmatrix.dot(vec.rrock._depth)
        
        return gvector

    def get_jmatrix(self,vec):
        """Returns J matrix filled with constant pressure coefficients."""

        jmatrix = csr(self.matrix)

        for vals,cond in zip(vec._W,self.wconds):
            jmatrix += self.set_jmatrix_constraint(jmatrix,vals,cond)

        for vals,cond in zip(vec._B,self.bconds):
            jmatrix += self.set_jmatrix_constraint(jmatrix,vals,cond)

        return jmatrix

    def get_qvector(self,vec):
        """Returns Q vector filled with constant rate and pressure coefficients."""

        qvector = csr(self.vector)

        for vals,cond in zip(vec._W,self.wconds):
            qvector += self.set_qvector_constraint(qvector,vals,cond)

        for vals,cond in zip(vec._B,self.bconds):
            qvector += self.set_qvector_constraint(qvector,vals,cond)

        return qvector

    @staticmethod
    def set_tmatrix_interblock(tmatrix:csr,vals,dneg,dpos):
        """Returns updated T matrix."""

        tmatrix += csr((vals,(dneg,dneg)),shape=tmatrix.shape)
        tmatrix -= csr((vals,(dneg,dpos)),shape=tmatrix.shape)

        tmatrix += csr((vals,(dpos,dpos)),shape=tmatrix.shape)
        tmatrix -= csr((vals,(dpos,dneg)),shape=tmatrix.shape)

        return tmatrix

    @staticmethod
    def set_jmatrix_constraint(jmatrix:csr,vals,cond):
        """Returns updated J matrix."""

        if cond.sort=="press":
            jmatrix += csr((vals,(cond.block,cond.block)),shape=jmatrix.shape)

        return jmatrix

    @staticmethod
    def set_qvector_constraint(qvector:csr,vals,cond):
        """Returns updated Q vector."""

        if cond.sort=="press":
            qvals = vals*cond._cond
        else:
            qvals = vals/vals.sum()*cond._cond

        qvector += csr((qvals,(cond.block,numpy.zeros(cond.block.size))),shape=qvector.shape)

        return qvector

    @property
    def matrix(self):
        """Returns the shape of matrices in the flow calculations."""
        return (self.edge.shape[0],self.edge.shape[0])

    @property
    def vector(self):
        """Returns the shape of vectors in the flow calculations."""
        return (self.edge.shape[0],1)
    
    