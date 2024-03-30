# import os

# import json

# filedir = os.path.dirname(__file__)

# filepath = os.path.join(filedir,"fluids.json")

# with open(filepath,"r") as jsonfile:
#     library = json.load(jsonfile)

from . import phaseg as Gas
from . import phaseo as Oil
from . import phasew as Water

from . import cpress as CapPress

from . import rperm as RelPerm

from ._phase import Phase

from ._single import OnePhase

from ._blackoil import BlackOil