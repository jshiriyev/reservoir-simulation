class Fluid():
    """
    Base Class that defines fluid properties at the
    given pressure and temperature.
    """

    def __init__(self,visc=None,rho=None,comp=None,fvf=None):
        
        self._visc  = visc*0.001
        self._rho   = rho*16.0185
        self._comp  = comp/6894.76
        self._fvf   = fvf

    def __call__(self,*args,**kwargs):

        return (self._visc,self._rho,self._comp,self._fvf)

    @property
    def visc(self):
        if self._visc is not None:
            return self._visc/0.001

    @property
    def rho(self):
        if self._rho is not None:
            return self._rho/16.0185

    @property
    def comp(self):
        if self._comp is not None:
            return self._comp*6894.76

    @property
    def fvf(self):
        return self._fvf
        



    