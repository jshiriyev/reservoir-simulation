# Analytical Solution
# from .analytical import linear
from .analytical import radial
# from .analytical import green

# 1. Input Data

# 2. Gridding
from .grids import Grids, GridDelta, GridRegular

# 3. Property Calculation
from .properties import Fluid
# from .properties import phaseg
# from .properties import cappres
# from .properties import relperm
from .properties import RRock

# 4. Configuration
from .configure import Well, Time

# 5. Simulation
from . import simulation

from .simulation import BaseSolver

# 6. Post-Processing