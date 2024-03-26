from numpy import ndarray

class Vector():

    def __init__(self,press,rrock,fluid,S:ndarray,X:ndarray,Y:ndarray,Z:ndarray,W:list,B:list):
        """
        Initialized vector interface with rrock and fluid properties,
        and the following items:
        
        S   : Storage Vector in SI units, (m**3)

        X   : Block transmissibility in x-direction
        Y   : Block transmissibility in y-direction
        Z   : Block transmissibility in z-direction

        W   : Well block transmissibility values
        B   : Boundary block transmissibility values

        Calculated parameters are:

        A   : Accumulation term (pore_volume)/(time_step)
              in SI units, (m**3)/(sec)
        C   : Accumulation times total compressibility term
              in SI units, (m**3)/(Pa*sec)
        """

        self._press = press # pressure values where properties are calculated

        self.rrock = rrock # pressure updated rrock properties
        self.fluid = fluid # pressure updated fluid properties

        self._S = S # storage

        self._X = X # x-transmissibility
        self._Y = Y # y-transmissibility
        self._Z = Z # z-transmissibility

        self._W = W # well transmissibility
        self._B = B # boundary transmissibility

    def set_A(self,tstep):
        """Sets accumulation vector for the class.
        tstep   : time step in the numerical calculations, sec"""
        self._A,self._tstep = self.get_A(self._S,tstep),tstep

    def set_C(self,tcomp=None):
        """Sets accumulation times total compressibility for the class.
        tcomp   : total compressibility in SI units, 1/Pa"""

        if tcomp is None:
            comp = self.rrock._comp+self.fluid._comp
        else:
            comp = numpy.asarray(tcomp).flatten()

        self._C,self._tcomp = self.get_C(self._A,comp),tcomp

    @property
    def S(self):
        """Converting from SI Units to Oil Field Units."""
        return self._S*(3.28084**3)

    @property
    def A(self):
        """Converting from SI Units to Oil Field Units."""
        return self._A*(3.28084**3)*(24*60*60)

    @property
    def C(self):
        """Converting from SI Units to Oil Field Units."""
        return self._C*(3.28084**3)*(24*60*60)*6894.76

    @property
    def X(self):
        """Converting from SI Units to Oil Field Units."""
        return self._X*(3.28084**3)*(24*60*60)*6894.76

    @property
    def Y(self):
        """Converting from SI Units to Oil Field Units."""
        return self._Y*(3.28084**3)*(24*60*60)*6894.76

    @property
    def Z(self):
        """Converting from SI Units to Oil Field Units."""
        return self._Z*(3.28084**3)*(24*60*60)*6894.76

    @staticmethod
    def get_A(storage,tstep):
        """Returns the accumulation vector for a given storage vector.
        tstep   : time step in the numerical calculations, sec"""
        return storage/tstep

    @staticmethod
    def get_C(accumulation,tcomp):
        """Returns accumulation times total compressibility vector.
        tcomp   : total compressibility in SI units, 1/Pa"""
        return accumulation*tcomp