class VanGenuchten():

    def __init__(self,irreducible=0,residual=0,n=2,m=2,gamma=0.02):

        self.irreducible = irreducible
        self.residual = residual

        self.n,self.m = n,m

        self.gamma = gamma

    def imbibition(self,sw):

        saturation = (sw-self.irreducible)/(1-self.irreducible-self.residual)

        return 1/self.gamma*(saturation**(-1/self.m)-1)**(1/self.n)

    def iimbibition(self,cappres):
        """Inverse imbibition calculation, saturation from capillary pressure."""

        saturation = ((self.gamma*cappres)**self.n+1)**(-self.m)

        return self.irreducible+saturation*(1-self.irreducible-self.residual)