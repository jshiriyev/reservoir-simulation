import numpy

class Vector():

    def __init__(self,S:numpy.ndarray,X:numpy.ndarray,Y:numpy.ndarray,Z:numpy.ndarray,W:list,B:list):
        """
        
        S   : Storage Vector in SI units, (m**3)

        X   : Block transmissibility in x-direction
        Y   : Block transmissibility in y-direction
        Z   : Block transmissibility in z-direction

        W   : Well block transmissibility values
        B   : Boundary block transmissibility values

        Calculated parameters are:

        A   : Accumulation term (pore_volume)/(time_step)
              in SI units, (m**3)/(sec)
        C   : A.dot(c) in SI units, (m**3)/(Pa*sec)
        """

        self._S = S

        self._X = X
        self._Y = Y
        self._Z = Z

        self._W = W
        self._B = B

    def set_A(self,tstep):
        """Sets (pore_volume)/(time_step).
        tstep   : time step in the numerical calculations, sec"""
        self._A = self.get_A(self._S,tstep)

    def set_C(self,tcomp):
        """Sets A.dot(c) in sparse matrix form.
        tcomp   : total compressibility in SI units, 1/Pa"""
        self._C = self.get_C(self._A,tcomp)

    @staticmethod
    def get_A(S,t):
        """Returns (pore_volume)/(time_step).
        t   : time step in the numerical calculations, sec"""
        return S/t

    @staticmethod
    def get_C(A,c):
        """Returns A.dot(c) in sparse matrix form.
        c   : total compressibility in SI units, 1/Pa"""
        return A*c

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