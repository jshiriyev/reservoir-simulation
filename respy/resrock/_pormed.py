import sys

if __name__ == "__main__":
    sys.path.append(r'C:\Users\3876yl\Documents\respy')

import numpy

from respy.rock._resrock import ResRock

from respy.utils._prop import Prop

class PorMed():

    def __init__(self,poro=None,xperm=None,yperm=None,zperm=None,comp=None):
        """
        """

        self.__poro  = Prop(poro)

        self.__xperm = Prop(xperm)
        self.__yperm = Prop(yperm)
        self.__zperm = Prop(zperm)

        self.__comp  = Prop(comp)
    
    @property
    def poro(self):
        return self.__poro()

    @property
    def xperm(self):
        return self.__xperm()

    @property
    def yperm(self):
        return self.__yperm()

    @property
    def zperm(self):
        return self.__zperm()

    @property
    def comp(self):
        return self.__comp()

    def __call__(self,press):

        poro  = self.__poro(press)
        xperm = self.__xperm(press)
        yperm = self.__yperm(press)
        zperm = self.__zperm(press)
        comp  = self.__comp(press)

        return ResRock(poro,xperm,yperm,zperm,comp)
    

if __name__ == "__main__":

    import numpy as np

    from respy.grids._delta import GridDelta

    grid = GridDelta((750,1000,125),(750,1000,1250),(20,))

    rock = PorMed(grid)

    rock.set_poro(np.array([1,2,3,4,5,6,7,8,9]))

    print(rock.grid.cube.poro)
    print(rock.grid.xmin.poro)

