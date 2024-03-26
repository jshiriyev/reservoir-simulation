import numpy

from scipy.sparse import csr_matrix as csr

from scipy.sparse import diags

from respy.solver._vector import Vector
from respy.solver._matrix import Matrix

class Build():

    def __init__(self,grid,wconds,bconds):

        self.grid = grid

        self.wconds = wconds
        self.bconds = bconds

    def __call__(self,press,rrock,fluid,vec:Vector):

        P = self.get_pmatrix(press)
        C = self.get_cmatrix(vec)
        T = self.get_tmatrix(vec)
        G = self.get_gmatrix(rrock,fluid,T)
        J = self.get_jmatrix(vec)
        Q = self.get_qmatrix(vec)

        return Matrix(P,C,T,G,J,Q)

    def get_pmatrix(self,press):
        """Returns P column matrix filled with pressure values."""
        pmatrix = numpy.asarray(press).flatten().reshape((-1,1))

        return pmatrix

    def get_cmatrix(self,vec:Vector):
        """Returns C matrix filled with diagonal values."""
        cmatrix = diags(vec._C,shape=self.matrix)

        return cmatrix

    def get_tmatrix(self,vec:Vector):
        """Returns T matrix filled with diagonal and offset values."""
        tmatrix = csr(self.matrix)

        tmatrix = self.set_tmatrix_interblock(
            tmatrix,vec._X,self.grid.xneg,self.grid.xpos)

        tmatrix = self.set_tmatrix_interblock(
            tmatrix,vec._Y,self.grid.yneg,self.grid.ypos)

        tmatrix = self.set_tmatrix_interblock(
            tmatrix,vec._Z,self.grid.zneg,self.grid.zpos)

        return tmatrix

    def get_gmatrix(self,rrock,fluid,tmatrix):
        """Returns G column matrix filled with gravity coefficients."""
        gmatrix = fluid._rho*9.807*tmatrix.dot(rrock._depth)
        
        return gmatrix

    def get_jmatrix(self,vec):
        """Returns J matrix filled with constant pressure coefficients."""

        jmatrix = csr(self.matrix)

        for vals,cond in zip(vec._W,self.wconds):
            jmatrix += self.set_jmatrix_constraint(jmatrix,vals,cond)

        for vals,cond in zip(vec._B,self.bconds):
            jmatrix += self.set_jmatrix_constraint(jmatrix,vals,cond)

        return jmatrix

    def get_qmatrix(self,vec):
        """Returns Q column matrix filled with constant rate and pressure coefficients."""

        qmatrix = csr(self.column)

        for vals,cond in zip(vec._W,self.wconds):
            qmatrix += self.set_qmatrix_constraint(qmatrix,vals,cond)

        for vals,cond in zip(vec._B,self.bconds):
            qmatrix += self.set_qmatrix_constraint(qmatrix,vals,cond)

        return qmatrix

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
    def set_qmatrix_constraint(qmatrix:csr,vals,cond):
        """Returns updated Q column matrix."""

        if cond.sort=="press":
            qvals = vals*cond._cond
        else:
            qvals = vals/vals.sum()*cond._cond

        qmatrix += csr((qvals,(cond.block,numpy.zeros(cond.block.size))),shape=qmatrix.shape)

        return qmatrix

    @property
    def matrix(self):
        """Returns the shape of matrices in the flow calculations."""
        return (self.edge.shape[0],)*2

    @property
    def column(self):
        """Returns the shape of column matrices in the flow calculations."""
        return (self.edge.shape[0],1)
    
    