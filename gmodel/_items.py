from dataclasses import dataclass

@dataclass
class Slot:
    """It is a well item dictionary."""
    name: str
    index: int = None
    field: str = None
    platform: str = None
    xhead: float = 0.0
    yhead: float = 0.0
    datum: float = 0.0
    status: str = "prospect"

@dataclass
class Zone:
    """It is a formation surface item dictionary."""
    name: str
    field: str = None
    color: str = "white"
    hatch: str = ".."

@dataclass
class Fault:
    """"""
    name: str
    index: int = None
    field: str = None

@dataclass
class Fracture:
    """"""
    index: int
    field: str = None
    length: float = None
    aperture: float = None

    # % The fracture segment is defined as a plane joining two node points
    # % (point1 and point2). The heigth of fracture plane is taken the same
    # % as reservoir thickness (it is not diffcult to model shorter planes).
    # % z-coordinate of the points is given as the reservoir depth.

    # nodeCoord,     # map,     # thickness
    # fracID,     # nodeID,    # numAfrac,    # numAnode,    # conductivity
    # point1,     # point2,    # Length,    # areatoreservoir,    # areatofracture
    # volume,    # center,    # signX,    # signY,    # azimuth

@dataclass
class Segment:
    """It is a reservoir segment with its petrophysical characteristics"""
    name: str = None