from respy.fluid import Fluid

from respy.fluid._critic import natural_gas
from respy.fluid._critic import gas_condensate
from respy.fluid._critic import composition

from respy.fluid._zfact import Dranchuk_Abu_Kassem
from respy.fluid._zfact import Dranchuk_Purvis_Robinson
from respy.fluid._zfact import Orkahub_Energy
from respy.fluid._zfact import Hall_Yarborough

from respy.fluid._visc import Carr_Kobayashi_Burrows
from respy.fluid._visc import Lee_Gonzalez_Eakin

class Gas():
    """
    The hydrocarbon gases that are normally found in a natural gas are
    methanes, ethanes, propanes, butanes, pentanes, and small amounts of
    hexanes and heavier. The nonhydrocarbon gases include carbon dioxide,
    hydrogen sulfide, and nitrogen.
    """

    def __init__(self,spgr:float=None,**kwargs):
        """
        Gas can be defined based on specific gravity or
            molecular composition.
        
        spgr    : specific gravity of the gas at
                  standard conditions.

        """

        self._spgr = spgr

        self._mw = 28.964*self._spgr/1000 #kg/mol

    @property
    def spgr(self):
        return self._spgr
    
    @property
    def mw(self):
        """Molecular Mass, lb/lbmol"""
        return 28.964*self._spgr

    @property
    def pcrit(self):
        """Function to Calculate Gas Critical Pressure in psia
        for high molecular weigth reservoir gases"""
        return 756.8-131*self._spgr-3.6*self._spgr**2

    @property
    def tcrit(self):
        """Function to Calculate Gas Critical Temperature in °R
        for high molecular weigth reservoir gases"""
        return 169.2+349.5*self._spgr-74*self._spgr**2

    @staticmethod
    def pred(P,Pc):
        return P/Pc

    @staticmethod
    def tred(T,Tc):
        return T/Tc

    def __call__(self,P,T,method=None):
        """
        P   : Pressure, psia
        T   : Temperature, °R
        """

        Pr = self.pred(P,self.pcrit)
        Tr = self.tred(T,self.tcrit)

        Z,Zprime = self.get_zfact(Pr,Tr,method=method)

        visc = self.get_visc(P,T,Z,M)
        rho  = self.get_rho(P,T,Z,M)
        comp = self.get_comp(Pr,Tr,Z,Zprime)/self.pcrit
        fvf  = self.get_fvf(P,T,Z)

        return Fluid(visc,rho,comp,fvf)

    @staticmethod
    def get_zfactor(Pr,Tr,method=None,derivative=False):

        if method=="HY":
            return Hall_Yarborough(Pr,Tr,derivative=derivative)
        elif method=="DAK":
            return Dranchuk_Abu_Kassem(Pr,Tr,derivative=derivative)
        elif method=="DPR":
            return Dranchuk_Purvis_Robinson(Pr,Tr,derivative=derivative)
        else:
            return Orkahub_Energy(Pr,Tr,derivative=derivative)

    @staticmethod
    def get_rho(P,T,M,Z):
        """Returns density value in lb/ft3"""
        return (P*M)/(Z*10.731577089016*T)

    @staticmethod
    def get_spvol(P,T,M,Z):
        """Returns specific volume in ft3/lb"""
        return (Z*10.731577089016*T)/(P*M)
    
    @staticmethod
    def get_comp(Pr,Tr,Z,Zprime):
        """Returns cpr=cg*Pc. This method should not be used
        at Tr<1.4 for 0.4<Pr<3.0."""
        return 1/(1+(0.27*Pr)/(Z**2*Tr)*Zprime)/Pr

    @staticmethod
    def get_fvf(P,T,Z):
        """Function to Calculate Gas Formation Volume Factor in ft3/scf"""
        return 0.02827*Z*T/P

    @staticmethod
    def get_fef(P,T,Z):
        """Returns gas expansion factor."""
        return P/(0.02827*Z*T)

    @staticmethod
    def get_visc(P,T,M,Z=None,method="LGE",**kwargs):
        """Returns Gas Viscosity in cp"""

        if method=="CKB":
            return Carr_Kobayashi_Burrows(M/28.964,T,**kwargs)(P,T)
        elif method=="LGE":
            return Lee_Gonzalez_Eakin(P,T,Z,M)

if __name__ == "__main__":

    fluids = gas(
        component=['methane','ethane'],
        molefraction=[0.2,0.4],
        moleweight=[12.5,None])

    print(fluids.standard['pressure']['SI'])

    print(fluids.composition['H2S'])

    print(fluids.composition.molefraction)
    print(fluids.composition.moleweight)
    print(fluids.composition)