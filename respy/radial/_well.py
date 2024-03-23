from dataclasses import dataclass

@dataclass(frozen=True)
class Well:
    """It is a well property dictionary."""
    rate     : float
    press    : float
    radius   : float