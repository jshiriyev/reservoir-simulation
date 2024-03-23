class Fluid():
    """
    Base Class that defines constant fluid properties at the
    given pressure and temperature.
    """

    def __init__(self,visc=None,rho=None,comp=None,fvf=None):
        """
        Initializes a fluid with certain viscosity, density,
        compressibility and formation volume factor:

        visc    : viscosity of fluid, cp
        rho     : density of fluid, lb/ft3
        comp    : compressibility of fluid, 1/psi
        fvf     : formation volume factor, ft3/scf

        """
        
        self._visc = None if visc is None else visc*0.001
        self._rho  = None if rho is None else rho*16.0185
        self._comp = None if comp is None else comp/6894.75729
        self._fvf  = fvf

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
            return self._comp*6894.75729

    @property
    def fvf(self):
        return self._fvf

if __name__ == "__main__":

    f = Fluid(5,5)

    print(f)

    print(f(45,47))

    print(f.visc,f._visc)
        



    