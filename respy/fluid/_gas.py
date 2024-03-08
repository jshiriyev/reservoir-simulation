from respy.fluid import Fluid

from respy.utils._prop import Prop

class Gas():
    """
    The hydrocarbon gases that are normally found in a natural gas are
    methanes, ethanes, propanes, butanes, pentanes, and small amounts of
    hexanes and heavier. The nonhydrocarbon gases include carbon dioxide,
    hydrogen sulfide, and nitrogen.
    """

    def __init__(self,spgr:float=None,*,crit=None,zfact=1.0,visc=None):
        """
        Real gas defined based on specific gravity and critical parameters:
        
        spgr    : specific gravity of the gas at
                  standard conditions.

        crit    : tuple of (pcrit,tcrit)

        zfact   : either constant value, or method to calculate z-factor
                  from pressure and temperature.

        visc    : either constant value, or method to calculate viscosity
                  from pressure and temperature.
        """

        self._spgr = spgr

        self._mw = 28.964*self._spgr/1000 # Molecular Mass, kg/mol

        crit = (None,None) if crit is None else crit

        pcrit,tcrit = crit

        self._pcrit = pcrit*6894.76
        self._tcrit = tcrit*5.0/9.0

        self.__zfact = Prop(zfact)
        self.__visc  = Prop(visc,0.001)

    @property
    def spgr(self):
        return self._spgr
    
    @property
    def mw(self):
        """Molecular Mass, lb/lbmol"""
        return self._mw*1000

    @property
    def pcrit(self):
        """Critical Pressure in psi"""
        return self._pcrit/6894.76

    @property
    def tcrit(self):
        """Critical Temperature in Rankine"""
        return self._tcrit/(5.0/9.0)
    
    def __call__(self,P,T):
        """
        P   : Pressure, psia
        T   : Temperature, °R
        """

        Z,Zprime = self.get_zfact(P)

        visc = self.get_visc(P,Z)
        rho  = self.get_rho(P,Z)
        comp = self.get_comp(P,Z,Zprime)
        fvf  = self.get_fvf(P,Z)

        return Fluid(visc,rho,comp,fvf)

    def pred(self,P):
        return P/self.pcrit

    def tred(self,T):
        return T/self.tcrit

    def get_zfact(self,P):

        # Pr = self.pred(P)
        # Tr = self.tred(T)
        
        # return self._zfact(Pr,Tr,derivative=True)

        if self.__zfact.const:
            zprime = 0

        return self.__zfact(P),zprime

    def get_visc(self,P,Z=1,**kwargs):
        """Calculates pressure and temperature dependent viscosity in cp"""
        # return self.__visc(P)
        return self.__visc(P)

    def get_rho(self,P,T,Z):
        """Returns density value in lb/ft3"""
        return (P*self.mw)/(Z*10.731577089016*T)

    def get_spvol(self,P,T,Z):
        """Returns specific volume in ft3/lb"""
        return (Z*10.731577089016*T)/(P*self.mw)

    def get_comp(self,P,T,Z,Zprime):
        """Returns cpr=cg*Pc. This method should not be used
        at Tr<1.4 for 0.4<Pr<3.0."""
        if zprime==0:
            return 1/P

        Pr,Tr = self.pred(P),self.tred(T)

        return 1/(1+(0.27*Pr)/(Z**2*Tr)*Zprime)/P

    def get_fvf(self,P,T,Z):
        """Function to Calculate Gas Formation Volume Factor in ft3/scf"""
        return 0.02827*Z*T/P

    def get_fef(self,P,T,Z):
        """Returns gas expansion factor."""
        return P/(0.02827*Z*T)

if __name__ == "__main__":

    gas = Gas(5,visc=5)

    print(gas.visc)

    # fluids = gas(
    #     component=['methane','ethane'],
    #     molefraction=[0.2,0.4],
    #     moleweight=[12.5,None])

    # print(fluids.standard['pressure']['SI'])

    # print(fluids.composition['H2S'])

    # print(fluids.composition.molefraction)
    # print(fluids.composition.moleweight)
    # print(fluids.composition)