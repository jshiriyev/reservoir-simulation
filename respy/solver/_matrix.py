from scipy.sparse import diags

class Matrix():

    def __init__(self,T,S,G,J,Q,t=None,c=None):
        """Inputs should be in SI units.

        T   : Transmissibility Matrix in SI units, (m**3)/(Pa*sec)
        S   : Storage Matrix in SI units, (m**3)
        G   : Gravity Vector in SI units, (m**3)/(sec)
        J   : Constant Pressure Matrix in SI units, (m**3)/(Pa*sec)
        Q   : Constant Flow Rate Vector in SI units, (m**3)/(sec)
        
        t   : time step in the numerical calculations, sec
        c   : total compressibility in SI units, 1/Pa

        Calculated parameters are:

        A   : (pore_volume)/(time_step) in SI units, (m**3)/(sec)
        B   : A.dot(c) in SI units, (m**3)/(Pa*sec)
        """

        self._T = T
        self._S = S
        self._G = G
        self._J = J
        self._Q = Q

        if t is not None:
            self._A = self.get_A(self._S,t)

        if c is not None:
            self._B = self.get_B(self._A,c)

    @staticmethod
    def get_A(S,t):
        """Returns (pore_volume)/(time_step)."""
        return S/t

    @staticmethod
    def get_B(A,c):
        """Returns A.dot(c) in sparse matrix form."""
        return A.dot(c)

    @property
    def T(self):
        """Converting from SI Units to Oil Field Units."""
        return self._T*(3.28084**3)*(24*60*60)*6894.76
    
    @property
    def S(self):
        """Converting from SI Units to Oil Field Units."""
        return self._S*(3.28084**3)

    @property
    def A(self):
        """Converting from SI Units to Oil Field Units."""
        return self._A*(3.28084**3)*(24*60*60)

    @property
    def B(self):
        return self._B*(3.28084**3)*(24*60*60)*6894.76
    
    @property
    def G(self):
        """Converting from SI Units to Oil Field Units."""
        return self._G*(3.28084**3)*(24*60*60)

    @property
    def J(self):
        """Converting from SI Units to Oil Field Units."""
        return self._J*(3.28084**3)*(24*60*60)*6894.76
    
    @property
    def Q(self):
        """Converting from SI Units to Oil Field Units."""
        return self._Q*(3.28084**3)*(24*60*60)

    def __call__(self,*args,**kwargs):

        return self

    def residual(self,P,Pprev):
        """Returns residual vector for the given n+1 values of P and (T,J,A,Q,G),
        and n values of P (Pprev) in SI units, (m**3)/(sec)"""
        return -csr.dot(self._T+self._J+self._B,P)+csr.dot(self._B,Pprev)+self._Q+self._G

    @property
    def pa2psi(self):
        """Pressure conversion factor from SI Units to Oil Field Units."""
        return 1/6894.76

    @property
    def psi2pa(self):
        """Pressure conversion factor from Oil Field Units to SI Units."""
        return 6894.76

    @property
    def shape(self):
        """Shape of the matrices of transmissibility calculations"""
        return self.T.shape

    @property
    def size(self):
        """Shape of the vectors of transmissibility calculations"""
        return self.T.size

    @staticmethod
    def implicit(Mnext,Pprev):
        """Implicit solution: matrices is defined at Pnext, and
        returned pressure is Pnext."""

        LHS = Mnext._T+Mnext._J+Mnext._B
        RHS = csr.dot(Mnext._B,Pprev)+Mnext._Q+Mnext._G

        return linalg.spsolve(LHS,RHS)

    @staticmethod
    def mixed(Mprev,Mnext,Pprev,Theta=0.5):
        """Mixed solution: matrices is defined at Pprev and Pnext, and
        returned pressure is Pnext."""

        LHS = (1-Theta)(Mnext._T+Mnext._J)+Mnext._B
        RHS = csr.dot(Mprev._B-Theta*(Mprev._T+Mprev._J),Pprev)+Mprev._Q+Mprev._G

        return linalg.spsolve(LHS,RHS)

    @staticmethod
    def explicit(Mprev,Pprev):
        """Explicit solution: matrices are defined at Pprev, and
        returned pressure is Pnext."""
        LHS = Mprev._B
        RHS = csr.dot(-(Mprev._T+Mprev._J),Pprev)+Mprev._Q+Mprev._G
        return Pprev+linalg.spsolve(LHS,RHS)

if __name__ == "__main__":

    pass