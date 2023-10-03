import matplotlib.pyplot as plt
import numpy as np

from scipy.sparse import csr_matrix as csr
from scipy.sparse import diags

from scipy.sparse.linalg import spsolve as sps

# import fluids

# from borepy.items import PorRock
# from borepy.items import Wells

from ._relperm import RelPerm

from ._cappres import BrooksCorey
from ._cappres import VanGenuchten
from ._cappres import JFunction
from ._cappres import ScanCurves

class ThreePhaseIMPES():

    def __init__(self):

        pass

class ThreePhaseSS():

    def __init__(self):

        pass

if __name__ == "__main__":

    import unittest

    from tests import test_porous_media

    unittest.main(test_porous_media)