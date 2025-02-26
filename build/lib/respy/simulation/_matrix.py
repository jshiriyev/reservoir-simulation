from scipy.sparse import csr_matrix as csr

class Matrix():

    def __init__(self,A,T,G,J,Q,P=None):
        """Inputs should be in SI units.

        A   : Accumulation-Compressibility Matrix,
              (m**3)/(Pa*sec)

        T   : Inter-Block Transmissibility Matrix,
              (m**3)/(Pa*sec)

        G   : Gravity Related Column Matrix,
              (m**3)/(sec)

        J   : Constant Pressure Constraint Matrix,
              (m**3)/(Pa*sec)

        Q   : Flow Rate Constraint Column Matrix,
              (m**3)/(sec)

        P   : Pressure Column Matrix Defining the Instance, (Pa)
        
        """

        self._A = A
        self._T = T
        self._G = G
        self._J = J
        self._Q = Q

        self._P = P

    def imppress(self,Pprev):
        """Implicit solution: self is defined at Pnext, and
        returned pressure is Pnext."""

        LHS = self._T+self._J+self._A
        RHS = self._A
        RHS = csr.dot(RHS,Pprev)+self._Q+self._G

        return linalg.spsolve(LHS,RHS)

    def impresid(self,Pprev):
        """Returns residual vector for the self defined at Pnext in
        SI units, (m**3)/(sec)"""

        LHS = self._T+self._J+self._A
        RHS = self._A
        RHS = csr.dot(RHS,Pprev)+self._Q+self._G

        return -csr.dot(LHS,self._P)+RHS

    def exppress(self):
        """Explicit solution: self is defined at Pprev, and
        returned pressure is Pnext."""

        LHS = self._A
        RHS = self._A-(self._T+self._J)
        RHS = csr.dot(RHS,self._P)+self._Q+self._G

        return linalg.spsolve(LHS,RHS)

    def expresid(self,Pnext):
        """Returns residual vector for the self defined at Pprev in
        SI units, (m**3)/(sec)"""

        LHS = self._A
        RHS = self._A-(self._T+self._J)
        RHS = csr.dot(RHS,self._P)+self._Q+self._G

        return -csr.dot(LHS,Pnext)+RHS
    
    @property
    def A(self):
        return self._A*(3.28084**3)*(24*60*60)*6894.76

    @property
    def T(self):
        """Converting from SI Units to Oil Field Units."""
        return self._T*(3.28084**3)*(24*60*60)*6894.76
    
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
    def P(self):
        """Converting from SI Units to Oil Field Units."""
        return self._P/6894.76
    

if __name__ == "__main__":

    pass