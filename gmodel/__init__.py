import os,sys,json

filedir = os.path.dirname(__file__)

from . import gstats
from . import mesh

# Geomodeling Elements
from ._items import Slot
from ._items import Zone
from ._items import Fault
from ._items import Fracture
from ._items import Segment


# Everything Associated with Well
wellspath = os.path.join(filedir,"_well.json")

with open(wellspath,"r") as wellsfile:
    json_wells = json.load(wellsfile)

from ._well import Well
from ._survey import Survey
from ._diagram import Diagram
from ._stock import Stock


# Everything Associated with Reservoir
formationpath = os.path.join(filedir,"_formation.json")

with open(formationpath,"r") as formationfile:
    json_zones = json.load(formationfile)

from ._surface import Surface
from ._formation import Formation
from ._formation import TopsView #must be depreciated later
from ._faults import Faults
from ._fracnet import Fractures
from ._reservoir import Reservoir