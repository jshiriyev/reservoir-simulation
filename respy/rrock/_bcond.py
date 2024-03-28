import sys

if __name__ == "__main__":
    # sys.path.append(r'C:\Users\javid.shiriyev\Documents\respy')
    sys.path.append(r'C:\Users\3876yl\Documents\respy')

class BoundCond():
    """
    It is a boundary condition object used in the simulator.
    """

    def __init__(self,face,*,start:float=0,stop:float=None,**kwargs):
        """
        face    : boundary: xmin, xmax, ymin, ymax, zmin, or zmax

        start   : start time for implementing the boundary condition, days
        stop    : stop time for implementing the boundary condition, days

        Assign only one of the following conditions:

        press   : constant pressure values, psi
        orate   : constant flow boundary condition, 0 = no flow, bbl/day
        wrate   : constant flow boundary condition, 0 = no flow, bbl/day
        grate   : constant flow boundary condition, 0 = no flow, ft3/day
        """
        self._face = face

        self._start  = start*(24*60*60)
        self._stop   = None if stop is None else stop*(24*60*60)

        for key,value in kwargs.items():
            if value is not None:
                self._sort = key
                if key == "press":
                    self._cond = value*6894.76
                elif key in ("orate","wrate"):
                    self._cond = value*(0.3048**3)/(24*60*60)*5.615
                elif key == "grate":
                    self._cond = value*(0.3048**3)/(24*60*60)
                break

        self._prod = None

    @property
    def face(self):
        return self._face

    @property
    def axis(self):
        return self._face[0]

    @property
    def start(self):
        return self._start/(24*60*60)

    @property
    def stop(self):
        return None if self._stop is None else self._stop/(24*60*60)

    @property
    def sort(self):
        return self._sort
    
    @property
    def cond(self):
        if self._sort == "press":
            return self._cond/6894.76
        elif self._sort in ("orate","wrate"):
            return self._cond*(24*60*60)/(0.3048**3)/5.615
        elif self._sort == "grate":
            return self._cond*(24*60*60)/(0.3048**3)

    @property
    def prod(self):
        return self._prod*(3.28084**3)*(24*60*60)*6894.76
    

if __name__ == "__main__":

    bcond = BoundCond("xmin",press=500)

    print(bcond.cond)
    print(bcond.face)
    print(bcond.sort)