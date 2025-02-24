# import os

# import json

# filedir = os.path.dirname(__file__)

# filepath = os.path.join(filedir,"fluids.json")

# with open(filepath,"r") as jsonfile:
#     library = json.load(jsonfile)

from . import phaseg
# from . import phaseo
# from . import phasew

# from . import cpress as CapPress

# from . import rperm as RelPerm

from ._fluid import Fluid

# from ._single import OnePhase

# from ._blackoil import BlackOil