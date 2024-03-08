class BoundCond():

    def __init__(self,face,**kwargs):
        """
        face    : boundary, xmin, xmax, ymin, ymax, zmin, or zmax

        Assign only one of the following conditions:

        press   : constant pressure values
        orate   : constant flow boundary condition, 0 = no flow
        wrate   : constant flow boundary condition, 0 = no flow
        grate   : constant flow boundary condition, 0 = no flow
        """
        self._face = face

        for key,value in kwargs.items():
            if value is not None:
                self._sort = key
                if key == "press":
                    self._cond = value*6894.76
                elif key in ("orate","wrate"):
                    self._cond = value*1.84013e-6
                elif key == "grate":
                    self._cond = value*3.27741e-7
                break

    @property
    def face(self):
        return self._face

    @property
    def sort(self):
        return self._sort
    
    @property
    def cond(self):
        if self._sort == "press":
            return self._cond/6894.76
        elif self._sort in ("orate","wrate"):
            return self._cond/1.84013e-6
        elif self._sort == "grate":
            return self._cond/3.27741e-7