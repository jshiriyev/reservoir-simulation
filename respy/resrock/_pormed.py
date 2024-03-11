import sys

if __name__ == "__main__":
    sys.path.append(r'C:\Users\3876yl\Documents\respy')

import numpy

from respy.rock._resrock import ResRock

from respy.utils._prop import Prop

class PorMed():

    def __init__(self,grid,**kwargs):
        """
        grid    : GridDelta or GridRegular instance, rectangular cuboid grids

        """

        self.grid = grid

    @property
    def size(self):
        """Number of grids in the reservoir."""
        return self.grid.cube.shape[0]

    def get_property(self,quality,coeff=1.,dtype=None):
        """coeff    : conversion factor"""

        quality = numpy.asarray(quality)

        if dtype is not None:
            quality = quality.astype(dtype)

        quality = quality.flatten()*coeff

        if quality.size==1:
            quality = quality.repeat(self.numtot)

        return quality.reshape((-1,1))

    def set_depth(self,depth):
        """Assigns the depth values in ft to the grids."""

        self.grid.cube.cvol.set_prop(
            depth=self.get_property(depth,coeff=0.3048,dtype=numpy.float_))

    def set_poro(self,poro):
        """Assigns the porosity values in fractions to the grids."""

        self.grid.cube.cvol.set_prop(
            poro=self.get_property(poro,dtype=numpy.float_))

    def set_perm(self,xperm,yperm=None,zperm=None,yreduce=1.,zreduce=1.):
        """Assigns the permeability values in mD to the grids."""

        self.grid.cube.cvol.set_prop(
            xperm=self.get_property(xperm,coeff=9.869233e-16,dtype=numpy.float_))

        self.grid.cube.cvol.set_prop(
            yperm=self.xperm*yreduce if yperm is None else self.get_property(
            yperm,coeff=9.869233e-16,dtype=numpy.float_))

        self.grid.cube.cvol.set_prop(
            zperm=self.xperm*zreduce if zperm is None else self.get_property(
                zperm,coeff=9.869233e-16,dtype=numpy.float_))

    def set_comp(self,P):
        return self.__comp(P)

    def __call__(self,press):

        poro  = self._poro
        xperm = self._xperm
        yperm = self._yperm
        zperm = self._zperm
        comp  = self._comp
        depth = self._depth

        return ResRock(poro,xperm,yperm,zperm,comp,depth)

    @property
    def depth(self):
        return self._depth/0.3048
    
    @property
    def poro(self):
        return self._poro

    @property
    def xperm(self):
        return self._xperm/9.869233e-16

    @property
    def yperm(self):
        return self._yperm/9.869233e-16

    @property
    def zperm(self):
        return self._zperm/9.869233e-16

if __name__ == "__main__":

    import numpy as np

    from respy.grids._delta import GridDelta

    grid = GridDelta((750,1000,125),(750,1000,1250),(20,))

    rock = PorMed(grid)

    rock.set_poro(np.array([1,2,3,4,5,6,7,8,9]))

    print(rock.grid.cube.poro)
    print(rock.grid.xmin.poro)

