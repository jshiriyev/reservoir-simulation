import numpy as np

from ._brooks_corey import BrooksCorey

class StonesI():
    """
    This Model Provides IMBIBITION Relative Permeability MODELS for a THREE-PHASE system.
    
    """

    def __init__(self,sorow=0.4,sorgo=0.4,swr=0.1,sgr=0.05,k0ro_wr=0.8,k0rw_or=0.3,k0ro_gr=0.8,k0rg_lr=0.3,no=2,nw=2,ng=2,som=None):
        """Initializes the class.

        Parameters:
        ----------
        swr     = residual water saturation
        sorow   = residual oil saturtaion in oil-water system
        sorgo   = residual oil saturation in gas-oil system
        sgr     = criticial gas saturation

        k0rw_or = water relative permeability at residual oil saturation in oil-water system
        k0ro_wr = oil relative permeability at residual water saturation in oil-water system
        k0ro_gr = oil relative permeability at residual gas saturation in gas-oil system
        k0rg_lr = gas relative permeability at residual liquid saturation in gas-oil system

        nw      = water exponent on relative permeability curve
        no      = oil exponent on relative permeability curve
        ng      = gas exponent on relative permeability curve

        som     = minimum oil saturation in three-phase system

        """
        self.swr     = swr
        self.sorow   = sorow
        self.sorgo   = sorgo
        self.sgr     = sgr

        self.k0rw_or = k0rw_or
        self.k0ro_wr = k0ro_wr
        self.k0ro_gr = k0ro_gr
        self.k0rg_lr = k0rg_lr

        self.nw      = nw
        self.no      = no
        self.ng      = ng
 
        self.som     = som

    @property
    def slrgo(self):
        """Getter for the residual liquid saturation in gas-oil system."""
        return self.swr+self.sorgo

    @property
    def som(self):
        """Getter for the minimum oil saturation in three-phase system."""
        return self._som
    
    @som.setter
    def som(self,value:float=None):
        """Setter for the minimum oil saturation in three-phase system."""
        self._som = np.min((self.sorow,self.sorgo)) if value is None else value

    def imbibition(self,Sw,So,Sg):
        """Computes three-phase relative permeabilities krw, kro and krg.

        Parameters:
        ----------
        Sw : water saturation
        So : oil saturation
        Sg : gas saturation

        """
        ow = BrooksCorey(self.swr,self.sorow,self.k0rw_or,self.k0ro_wr,self.nw,self.no)
        go = BrooksCorey(self.slrgo,self.sgr,self.k0ro_gr,self.k0rg_lr,self.no,self.ng)

        movable_o = So-self.som
        movable_w = Sw-self.swr
        movable_g = Sg

        movable_f = 1-self.swr-self.som

        So_star = movable_o/movable_f
        Sw_star = movable_w/movable_f
        Sg_star = movable_g/movable_f

        kroow,krw = self.water_oil(Sw,So)
        krogo,krg = self.gas_oil(Sw,So,Sg)

        beta_w = (kroow)/(1-Sw_star)
        beta_g = (krogo)/(1-Sg_star)

        kro = So_star*beta_w*beta_g

        return kro,krw,krg

if __name__ == "__main__":

    import unittest

    from fluidflow.pormed.tests.conrelation import TestRelativePermeability

    unittest.main()

    # rp = relative_permeability_balhoff(
    #     Swi=0.2,
    #     Swr=0.2,
    #     krwo=0.2,
    #     kroo=1.0,
    #     nw=3,
    #     no=3,
    #     )

    # Sw = 0.2001

    # rp.system2phase(Sw=Sw,model="oil-water")

    # print(rp.krw)