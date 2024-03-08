class WellCond:
    """
    It is a well condition object used in the simulator.
    """

    def __init__(self,block:tuple,axis:str,radius:float,skin:float,start:float,end:float,**kwargs):
        """
        block   : block indices containing the well
        axis    : (z) vertical or (x,y) horizontal well

        radius  : well radius, ft
        skin    : skin factor of the well, dimensionless

        start   : start time for implementing the condition, days
        end     : end time for implementing the condition, days

        Assign only one of the following conditions:

        bhp     : constant bottom hole pressure
        whp     : constant well head pressure

        orate   : constant oil rate
        wrate   : constant water rate
        grate   : constant gas rate
        """

        self._block     = block
        self._axis      = axis
        self._radius    = radius*0.3048
        self._skin      = skin
        self._start     = start*(24*60*60)
        self._end       = end*(24*60*60)

        for key,value in kwargs.items():
            if value is not None:
                self._sort = key
                if key == "press":
                    self._cond = value*6894.76
                elif key in ("orate","grate"):
                    self._cond = value*1.84013e-6
                elif key == "grate":
                    self._cond = value*3.27741e-7
                break

    @property
    def block(self):
        return self._block

    @property
    def axis(self):
        return self._axis

    @property
    def radius(self):
        return self._radius/0.3048

    @property
    def skin(self):
        return self._skin

    @property
    def start(self):
        return self._start/(24*60*60)

    @property
    def end(self):
        return self._end/(24*60*60)

    @property
    def sort(self):
        return self._sort
    
    @property
    def cond(self):
        if self._sort in ("bhp","whp"):
            return self._cond/6894.76
        elif self._sort in ("orate","wrate"):
            return self._cond/1.84013e-6
        elif self._sort == "grate":
            return self._cond/3.27741e-7