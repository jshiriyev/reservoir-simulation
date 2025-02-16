import sys

if __name__ == "__main__":
    # sys.path.append(r'C:\Users\javid.shiriyev\Documents\respy')
    sys.path.append(r'C:\Users\3876yl\Documents\respy')

import numpy

from scipy.sparse import csr_matrix as csr

from respy.solver._vector import Vector
from respy.solver._matrix import Matrix

class Build():

    def __init__(self,grid):

        self.grid = grid

    def __call__(self,vec:Vector):

        A = self.get_A(vec)
        T = self.get_T(vec)
        G = self.get_G(vec,T)
        J = self.get_J(vec)
        Q = self.get_Q(vec)

        return Matrix(A,T,G,J,Q)

    def get_A(self,vec:Vector):
        """Returns A matrix filled with diagonal values."""
        return csr((vec._A,(self.grid.rows,self.grid.rows)),shape=self.matrix)

    def get_T(self,vec:Vector):
        """Returns T matrix filled with diagonal and offset values."""
        tmatrix = csr(self.matrix)

        tmatrix = self.add_tmatrix_interblock(
            tmatrix,vec._X,self.grid.xneg,self.grid.xpos)

        tmatrix = self.add_tmatrix_interblock(
            tmatrix,vec._Y,self.grid.yneg,self.grid.ypos)

        tmatrix = self.add_tmatrix_interblock(
            tmatrix,vec._Z,self.grid.zneg,self.grid.zpos)

        return tmatrix

    def get_G(self,vec:Vector,tmatrix):
        """Returns G column matrix filled with gravity related terms."""
        return tmatrix.dot(vec._H.reshape((-1,1)))

    def get_J(self,vec):
        """Returns J matrix filled with constant pressure constraints on diagonal."""

        jmatrix = csr(self.matrix)

        for wcond in vec._W:
            jmatrix += self.add_jmatrix_constraint(jmatrix,wcond)

        for bcond in vec._B:
            jmatrix += self.add_jmatrix_constraint(jmatrix,bcond)

        return jmatrix

    def get_Q(self,vec):
        """Returns Q column matrix filled with constraints."""

        qmatrix = csr(self.column)

        for wcond in vec._W:
            qmatrix += self.add_qmatrix_constraint(qmatrix,wcond)

        for bcond in vec._B:
            qmatrix += self.add_qmatrix_constraint(qmatrix,bcond)

        return qmatrix

    @staticmethod
    def add_tmatrix_interblock(tmatrix:csr,vals,dneg,dpos):
        """Returns updated T matrix."""

        tmatrix += csr((vals,(dneg,dneg)),shape=tmatrix.shape)
        tmatrix -= csr((vals,(dneg,dpos)),shape=tmatrix.shape)

        tmatrix += csr((vals,(dpos,dpos)),shape=tmatrix.shape)
        tmatrix -= csr((vals,(dpos,dneg)),shape=tmatrix.shape)

        return tmatrix

    @staticmethod
    def add_jmatrix_constraint(jmatrix:csr,cond):
        """Returns updated J matrix."""

        if cond.sort=="press":
            jmatrix += csr((cond._prod,(cond.block,cond.block)),shape=jmatrix.shape)

        return jmatrix

    @staticmethod
    def add_qmatrix_constraint(qmatrix:csr,cond):
        """Returns updated Q column matrix."""

        qvector = cond._prod*cond._cond

        if cond.sort!="press":
            qvector /= cond._prod.sum()

        qmatrix += csr((qvector,(cond.block,numpy.zeros(cond.block.size))),shape=qmatrix.shape)

        return qmatrix

    @property
    def matrix(self):
        """Returns the shape of square matrices in the flow calculations."""
        return (self.grid.nums,)*2

    @property
    def column(self):
        """Returns the shape of column matrices in the flow calculations."""
        return (self.grid.nums,1)

    @staticmethod
    def get_diag(array):
        return csr((array,(self.grid.rows,self.grid.rows)),shape=self.matrix)

if __name__ == "__main__":

    pass