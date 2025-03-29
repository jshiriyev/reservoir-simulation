import numpy as np

from ._stones_I import StonesI

class StonesII(StonesI):
    """
    This Model Provides IMBIBITION Relative Permeability MODELS for a THREE-PHASE system.

    """

    def __init__(self,**kwargs):

        super().__init__(**kwargs)

    def kro(self,krw,kro_ow,kro_go,krg):
        """Returns relative oil permeability based on Stone's Model II."""
        krw = np.ravel(krw)
        krg = np.ravel(krg)

        beta_w = kro_ow/self.k0ro_ow+krw
        beta_g = kro_go/self.k0ro_ow+krg

        return self.k0ro_ow*(beta_w*beta_g-(krw+krg))

    def get(self,sw:np.ndarray,so:np.ndarray,sg:np.ndarray):
        """Computes three-phase relative permeabilities krw, kro and krg.

        Parameters:
        ----------
        sw : water saturation
        so : oil saturation
        sg : gas saturation

        """
        krw,kro_ow = self.ow.get(sw)
        kro_go,krg = self.go.get(sw+so)

        kro = self.kro(krw,kro_ow,kro_go,krg)

        sod = self.sod(so)

        kro[sod<0] = 0
        kro[sod>1] = self.k0ro_ow

        return krw,kro,krg