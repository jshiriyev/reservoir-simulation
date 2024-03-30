import sys

if __name__ == "__main__":
    # sys.path.append(r'C:\Users\javid.shiriyev\Documents\respy')
    sys.path.append(r'C:\Users\3876yl\Documents\respy')

import numpy

class Phase():
    """
    A class that defines constant fluid properties at the
    given pressure and temperature.
    """

    def __init__(self,*,visc=None,rho=None,comp=None,fvf=None,rperm=1.0,press=None):
        """
        Initializes a fluid with certain viscosity, density,
        compressibility and formation volume factor:

        visc    : viscosity of fluid, cp
        rho     : density of fluid, lb/ft3
        comp    : compressibility of fluid, 1/psi
        fvf     : formation volume factor, ft3/scf

        rperm   : relative permeability value, default is 1, dimensionless
        press   : pressure at which properties are defined; if it is None,
                  properties are pressure independent, psi

        For any given pressure and temperature values.

        """
        
        self._visc  = self.get_prop(visc,0.001)
        self._rho   = self.get_prop(rho,16.0185)

        self._grad  = self._rho*9.807
        
        self._comp  = self.get_prop(comp,1/6894.76)
        self._fvf   = self.get_prop(fvf)

        self._rperm = self.get_prop(rperm)
        self._press = self.get_prop(press,6894.76)        

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


    