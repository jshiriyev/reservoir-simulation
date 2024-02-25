from dataclasses import dataclass

import numpy

@dataclass(frozen=True)
class Fluid:
    
    density		: float
    viscosity   : float

    def __post_init__(self):

        object.__setattr__(self,'_density',self.density*16.0185)
        object.__setattr__(self,'_viscosity',self.viscosity*0.001)
