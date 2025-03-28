import numpy as np

class StonesII():
    """
    This Model Provides IMBIBITION Relative Permeability MODELS for a THREE-PHASE system.

    Som     = minimum oil saturation in three-phase system
    Sorow   = residual oil saturtaion in oil-water system
    Sorgo   = residual oil saturation in gas-oil system
    Swc     = connate water saturation
    Sgc     = criticial gas saturation
    krowc   = oil relative permeability at connate water saturaton in oil-water system
    krwor   = water relative permeability at the residual oil saturation in oil-water system
    krogc   = oil relative permeability at critical gas saturation in gas-oil system
    krglc   = gas relative permeability at critical liquid saturation in gas-oil system
    no      = oil exponent on relative permeability curve
    nw      = water exponent on relative permeability curve
    ng      = gas exponent on relative permeability curve
    
    """

    def __init__(self,Sorow=0.4,Sorgo=0.4,Swc=0.1,Sgc=0.05,krowc=0.8,krwor=0.3,krogc=0.8,krglc=0.3,no=2,nw=2,ng=2,Som=None):

        self.Sorow  = Sorow
        self.Sorgo  = Sorgo
        self.Swc    = Swc
        self.Sgc    = Sgc
        self.krowc  = krowc
        self.krwor  = krwor
        self.krogc  = krogc
        self.krglc  = krglc
        self.no     = no
        self.nw     = nw
        self.ng     = ng

        self.Som    = Som

    def _stones_model_II(self,Sw,So,Sg):
        """
        Sw      = water saturation
        So      = oil saturation
        Sg      = gas saturation

        """
        movable_o = So-self.Som
        movable_w = Sw-self.Swc
        movable_g = Sg

        movable_f = 1-self.Swc-self.Som

        So_star = movable_o/movable_f
        Sw_star = movable_w/movable_f
        Sg_star = movable_g/movable_f

        kroow,krw = self.water_oil(Sw,So)
        krogo,krg = self.gas_oil(Sw,So,Sg)

        kro = self.krowc*((kroow/self.krowc+krw)*(krogo/self.krowc+krg)-(krw+krg))

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