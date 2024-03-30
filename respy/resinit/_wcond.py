import sys

if __name__ == "__main__":
    # sys.path.append(r'C:\Users\javid.shiriyev\Documents\respy')
    sys.path.append(r'C:\Users\3876yl\Documents\respy')

class WellCond():
    """
    It is a well condition object used in the simulator.
    """

    def __init__(self,radius:float,*,block:tuple=None,axis:str="z",skin:float=0,start:float=0,stop:float=None,**kwargs):
        """
        radius  : well radius, ft

        block   : block indices containing the well

        axis    : (z) vertical or (x,y) horizontal well

        skin    : skin factor of the well, dimensionless

        start   : start time for implementing the well condition, days
        stop    : stop time for implementing the well condition, days

        Assign only one of the following conditions:

        press   : constant bottom hole pressure, psi
        orate   : constant oil rate
        wrate   : constant water rate
        grate   : constant gas rate
        """

        self._radius = radius*0.3048

        self._block  = block

        self._axis   = axis

        self._skin   = skin

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

        self._prod   = None
    
    @property
    def radius(self):
        return self._radius/0.3048

    @property
    def block(self):
        return self._block

    @property
    def axis(self):
        return self._axis

    @property
    def skin(self):
        return self._skin

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

    wcond = WellCond(0.5,(3,),"z",orate=500,start=3)

    print(wcond.sort,wcond.cond)

    print(wcond.start,wcond._start)