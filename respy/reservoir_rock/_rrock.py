import sys

if __name__ == "__main__":
    sys.path.append(r'C:\Users\3876yl\Documents\respy')

import numpy

from respy.rrock._pmed import PorMed

class ResRock():

    def __init__(self,*args,**kwargs):
        """
        Reservoir Rock System
        """

        self.rrock = PorMed(*args,**kwargs)

    def __call__(self,press):

        if callable(self.rrock):
            return self.rrock(press)

        self.rrock._press = press

        return self.rrock

    @property
    def isstatic(self):
        return not callable(self.rrock)

    @property
    def isdynamic(self):
        return callable(self.rrock)

if __name__ == "__main__":

    import numpy as np

    from respy.grids._delta import GridDelta

    grid = GridDelta((750,1000,125),(750,1000,1250),(20,))

    rock = PorMed(grid)

    rock.set_poro(np.array([1,2,3,4,5,6,7,8,9]))

    print(rock.grid.cube.poro)
    print(rock.grid.xmin.poro)

