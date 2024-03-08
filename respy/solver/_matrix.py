from scipy.sparse import diags

class Matrix():

    def __init__(self,T,S,G,J,Q,tcomp=None,tstep=None):
        """Inputs should be in SI units.

        T  		: Transmissibility Matrix
        S 		: Storage Matrix
        G 		: Gravity Vector
        J 		: Constant Pressure Matrix
        Q 		: Constant Flow Rate Vector
		
		tcomp 	: total compressibility, 1/Pa
		tstep 	: time step in the numerical calculations, sec

        """

        self._T = T
        self._S = S
        self._G = G
        self._J = J
        self._Q = Q

        tcomp = numpy.asarray(tcomp).flatten()

        if tcomp.size == 1:
        	tcomp = tcomp.repeat(self._T.shape[0])

        self._tcomp = diags(tcomp)
        self._tstep = tstep

    @property
    def T(self):
        """Converting from SI Units to Oil Field Units."""
        return self._T*(3.28084**3*(24*60*60)*6894.76)
    
    @property
    def S(self):
        """Converting from SI Units to Oil Field Units."""
        return self._S*(3.28084**3)

    @property
    def tcomp(self):
        """Converting from SI Units to Oil Field Units."""
        return self._tcomp/6894.76

    @property
    def tstep(self):
        """Converting from SI Units to Oil Field Units."""
        return self._tstep/(24*60*60)

    @property
    def A(self):
        """Returning in Oil Field Units."""
        return self.S.dot(self.tcomp)/self.tstep

    @property
    def _A(self):
        """Returning in SI Units."""
        return self._S.dot(self._tcomp)/self._tstep

    @property
    def G(self):
        """Converting from SI Units to Oil Field Units."""
        return self._G*(3.28084**3*(24*60*60))

    @property
    def J(self):
        """Converting from SI Units to Oil Field Units."""
        return self._J*(3.28084**3*(24*60*60)*6894.76)
    
    @property
    def Q(self):
        """Converting from SI Units to Oil Field Units."""
        return self._Q*(3.28084**3*(24*60*60))

    @property
    def all(self):
        return (self.T,self.A,self.G,self.J,self.Q)

    @property
    def _all(self):
        return (self._T,self._A,self._G,self._J,self._Q)

    def residual(self,P,Pn):
        """Returns residual vector for the given n+1 or n-1 values of
        P and (T,J,A,Q,G), and n values of (P,) in Oil Field Units"""
        T,A,G,J,Q = self.all
        return -csr.dot(T+J+Act,P)+csr.dot(Act,Pn)+Q+G

    def _residual(self,P,Pn):
        """Returns residual vector for the given n+1 or n-1 values of
        P and (T,J,A,Q,G), and n values of (P,) in SI units"""
        T,A,G,J,Q = self._all
        return -csr.dot(T+J+Act,P)+csr.dot(Act,Pn)+Q+G

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