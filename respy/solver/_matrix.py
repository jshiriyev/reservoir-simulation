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
    def implicit(m,P):
        """Implicit solution of one-phase flow."""
        return linalg.spsolve(m.T+m.J+m.A,csr.dot(m.A,P)+m.Q+m.G)

    @staticmethod
    def mixed(m,P,theta):
        """Mixed solution of one-phase flow."""

        LHS = (1-theta)(m.T+m.J)+m.A
        RHS = csr.dot(m.A-theta*(m.T+m.J),P)+m.Q+m.G

        return linalg.spsolve(LHS,RHS)

    @staticmethod
    def explicit(m,P):
        """Explicit solution of one-phase flow."""
        return P+linalg.spsolve(m.A,csr.dot(-(m.T+m.J),P)+m.Q+m.G)

if __name__ == "__main__":

    pass