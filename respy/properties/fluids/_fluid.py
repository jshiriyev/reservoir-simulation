import numpy

class Fluid():
    """
    A class that defines constant fluid properties at the
    given pressure and temperature.
    """

    def __init__(self,*,visc=None,rho=None,comp=None,fvf=None,press=None,satur=1.0,rperm=1.0):
        """
        Initializes a fluid with certain viscosity, density,
        compressibility and formation volume factor:

        visc    : viscosity of fluid, cp
        rho     : density of fluid, lb/ft3
        comp    : compressibility of fluid, 1/psi
        fvf     : formation volume factor, ft3/scf

        press   : pressure at which properties are defined; if it is None,
                  properties are pressure independent, psi
        
        satur   : saturation value, default is 1, dimensionless
        rperm   : relative permeability value, default is 1, dimensionless

        For any given pressure and temperature values.

        """
        
        self._visc  = self.get_prop(visc,0.001)
        self._rho   = self.get_prop(rho,16.0185)

        self._grad  = self._rho*9.807
        
        self._comp  = self.get_prop(comp,1/6894.76)
        self._fvf   = self.get_prop(fvf)

        self._press = self.get_prop(press,6894.76)
        self._satur = self.get_prop(satur)
        self._rperm = self.get_prop(rperm)

        self._mobil = self.get_mobil()

    @property
    def visc(self):
        if self._visc is not None:
            return self._visc/0.001

    @property
    def rho(self):
        if self._rho is not None:
            return self._rho/16.0185

    @property
    def grad(self):
        return self._grad

    @property
    def comp(self):
        if self._comp is not None:
            return self._comp*6894.76

    @property
    def fvf(self):
        return self._fvf

    @property
    def press(self):
        if self._press is not None:
            return self._press/6894.76

    @property
    def mobil(self):
        if self._mobil is not None:
            return self._mobil*6894.76*(24*60*60)

    def get_mobil(self):
        try:
            return (self._rperm)/(self._visc*self._fvf)
        except TypeError:
            return

    @staticmethod
    def get_prop(prop,conv=1.):
        if prop is not None:
            return numpy.asarray(prop).astype(numpy.float_)*conv

if __name__ == "__main__":

    f = Fluid(5,5)

    print(f)

    print(f.visc,f._visc)
    
    print(f.visc.dtype)

    print(f.comp)

    print(callable(f))


    