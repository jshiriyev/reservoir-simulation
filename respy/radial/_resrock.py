from dataclasses import dataclass

@dataclass(frozen=True)
class ResRock:
    """It is a well property dictionary.

    perm    : reservoir permeability in mD
    height  : reservoir thickness in ft
    radius  : reservoir exterior boundary in ft

    """
    perm     : float # *9.869233e-16
    radius   : float # *0.3048
    height   : float = 1 # *0.3048

    @property
    def perm(self):
        return self._perm/9.869233e-16

    @property
    def height(self):
        return self._height/0.3048

    @property
    def radius(self):
        return self._redge/0.3048