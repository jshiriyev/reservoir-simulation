from scipy.sparse import csr_matrix as csr

class Matrix():

    def __init__(self,P,C,T,G,J,Q):
        """Inputs should be in SI units.

        P   : Pressure values at which class matrices are
              calculated, (Pa)

        C   : Accumulation and compressibility multiplication
              in SI units, (m**3)/(Pa*sec)
        T   : Transmissibility Matrix in SI units, (m**3)/(Pa*sec)
        G   : Gravity Vector in SI units, (m**3)/(sec)
        J   : Constant Pressure Matrix in SI units, (m**3)/(Pa*sec)
        Q   : Constant Flow Rate Vector in SI units, (m**3)/(sec)
        
        """

        self._P = self.get_matpress(P)

        self._C = C
        self._T = T
        self._G = G
        self._J = J
        self._Q = Q

    def implicit_pressure(self,Pprev):
        """Implicit solution: self is defined at Pnext, and
        returned pressure is Pnext."""

        LHS = self._T+self._J+self._C
        RHS = self._C
        RHS = csr.dot(RHS,Pprev)+self._Q+self._G

        return linalg.spsolve(LHS,RHS)

    def implicit_residual(self,Pprev):
        """Returns residual vector for the self defined at Pnext in
        SI units, (m**3)/(sec)"""

        LHS = self._T+self._J+self._C
        RHS = self._C
        RHS = csr.dot(RHS,Pprev)+self._Q+self._G

        return -csr.dot(LHS,self._P)+RHS

    def explicit_pressure(self):
        """Explicit solution: self is defined at Pprev, and
        returned pressure is Pnext."""

        LHS = self._C
        RHS = self._C-(self._T+self._J)
        RHS = csr.dot(RHS,self._P)+self._Q+self._G

        return linalg.spsolve(LHS,RHS)

    def explicit_residual(self,Pnext):
        """Returns residual vector for the self defined at Pprev in
        SI units, (m**3)/(sec)"""

        LHS = self._C
        RHS = self._C-(self._T+self._J)
        RHS = csr.dot(RHS,self._P)+self._Q+self._G

        return -csr.dot(LHS,Pnext)+RHS

    @property
    def P(self):
        """Converting from SI Units to Oil Field Units."""
        return self._P/6894.76
    
    @property
    def T(self):
        """Converting from SI Units to Oil Field Units."""
        return self._T*(3.28084**3)*(24*60*60)*6894.76

    @property
    def C(self):
        return self._C*(3.28084**3)*(24*60*60)*6894.76
    
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

    @property
    def shape(self):
        """Shape of the matrices of transmissibility calculations"""
        return self.T.shape

if __name__ == "__main__":

    pass