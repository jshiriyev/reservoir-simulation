import sys

if __name__ == "__main__":
    sys.path.append(r'C:\Users\3876yl\Documents\respy')

import numpy

from respy.rock._resrock import ResRock

class PorMed():

    def __init__(self,poro=None,xperm=None,yperm=None,zperm=None,comp=None):
        """
        """

        self.__poro = poro
        self.__perm = perm
        self.__comp = comp
    
    def poro(self,press):
        if callable(self.__poro):
            return self.__poro(press)
        return self.__poro

    def xperm(self,press):
        if callable(self.__perm):
            return self.__perm(press)
        return self.__perm

    def yperm(self,press,yreduce=1.):
        if callable(self.__perm):
            return self.__perm(press)/yreduce
        return self.__perm/yreduce

    def zperm(self,press,zreduce=1):
        if callable(self.__perm):
            return self.__perm(press)/zreduce
        return self.__perm/zreduce

    def comp(self,press):
        if callable(self.__comp):
            return self.__comp(press)
        return self.__comp

    def __call__(self,press):

        poro  = self.__poro(press)
        xperm = self.__xperm(press)
        yperm = self.__yperm(press)
        zperm = self.__zperm(press)
        comp  = self.__comp(press)

        return ResRock(
            xperm,yperm=yperm,zperm=zperm,poro=poro,comp=comp,press=press)

if __name__ == "__main__":

    import numpy as np

    from respy.grids._delta import GridDelta

    grid = GridDelta((750,1000,125),(750,1000,1250),(20,))

    rock = PorMed(grid)

    rock.set_poro(np.array([1,2,3,4,5,6,7,8,9]))

    print(rock.grid.cube.poro)
    print(rock.grid.xmin.poro)

