import sys

if __name__ == "__main__":
    # sys.path.append(r'C:\Users\javid.shiriyev\Documents\respy')
    sys.path.append(r'C:\Users\3876yl\Documents\respy')

import numpy

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

        For any given pressure and temperature values.

        """
        
        self._visc = self.set_prop(visc,0.001)
        self._rho  = self.set_prop(rho,16.0185)
        self._comp = self.set_prop(comp,1/6894.76)
        self._fvf  = self.set_prop(fvf)

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

    @staticmethod
    def set_prop(prop,conv=1.):
        if prop is not None:
            return numpy.asarray(prop).astype(numpy.float_)*conv

if __name__ == "__main__":

    f = Fluid(5,5)

    print(f)

    print(f.visc,f._visc)
    
    print(f.visc.dtype)

    print(f.comp)

    print(callable(f))


    