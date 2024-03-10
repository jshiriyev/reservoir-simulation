# import os

# import json

# filedir = os.path.dirname(__file__)

# filepath = os.path.join(filedir,"fluids.json")

# with open(filepath,"r") as jsonfile:
#     library = json.load(jsonfile)

# I may need constitutive models to describe the relationship
# between shear stress and strain rate
# Newtonian and non newtonian fluids

from ._fluid import Fluid

from ._gas import Gas
from ._oil import Oil
from ._water import Water

from ._zfact import DAK
from ._zfact import DPR
from ._zfact import OE
from ._zfact import HY