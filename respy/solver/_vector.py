class Vector():

    def __init__(self,X,Y,Z,S,W,B,t=None,c=None):
        """
        X   : Block transmissibility in x-direction
        Y   : Block transmissibility in y-direction
        Z   : Block transmissibility in z-direction

        S   : Storage Matrix in SI units, (m**3)

        t   : time step in the numerical calculations, sec
        c   : total compressibility in SI units, 1/Pa

        Calculated parameters are:

        A   : Accumulation term (pore_volume)/(time_step)
              in SI units, (m**3)/(sec)
        C   : A.dot(c) in SI units, (m**3)/(Pa*sec)
        """

        self._X = X
        self._Y = Y
        self._Z = Z

        self._S = S

        if t is not None:
            self._A = self.get_A(self._S,t)

        if c is not None:
            self._C = self.get_C(self._A,c)

        self._W = W
        self._B = B

    @staticmethod
    def get_A(S,t):
        """Returns (pore_volume)/(time_step)."""
        return S/t

    @staticmethod
    def get_C(A,c):
        """Returns A.dot(c) in sparse matrix form."""
        return A.dot(c)

    @property
    def S(self):
        """Converting from SI Units to Oil Field Units."""
        return self._S*(3.28084**3)

    @property
    def A(self):
        """Converting from SI Units to Oil Field Units."""
        return self._A*(3.28084**3)*(24*60*60)