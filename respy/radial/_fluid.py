from dataclasses import dataclass

@dataclass(frozen=True)
class Fluid:
    """It is a well property dictionary."""
    visc     : float
    fvf      : float = 1
    comp     : float = 0