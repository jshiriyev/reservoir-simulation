import numpy as np

from _brooks_corey import BrooksCorey

class StonesI():
    """
    This Model Provides IMBIBITION Relative Permeability MODELS for a THREE-PHASE system.
    
    """

    def __init__(self,
        swr:float = 0.1, sor_ow:float = 0.4, sor_go:float = 0.4, sgr:float = 0.05,
        k0rw_ow:float = 0.3, k0ro_ow:float = 0.8, k0ro_go:float = 0.8, k0rg_go:float = 0.3,
        nw:float = 2, no:float = 2, ng:float = 2
        ):
        """Initializes the class.

        Parameters:
        ----------
        swr     = residual water saturation
        sor_ow  = residual oil saturtaion in oil-water system
        sor_go  = residual oil saturation in gas-oil system
        sgr     = criticial gas saturation

        k0rw_ow = water relative permeability at residual oil saturation in oil-water system
        k0ro_ow = oil relative permeability at residual water saturation in oil-water system
        k0ro_go = oil relative permeability at residual gas saturation in gas-oil system
        k0rg_go = gas relative permeability at residual liquid saturation in gas-oil system

        nw      = water exponent on relative permeability curve
        no      = oil exponent on relative permeability curve
        ng      = gas exponent on relative permeability curve

        """
        self.swr     = swr
        self.sor_ow  = sor_ow
        self.sor_go  = sor_go
        self.sgr     = sgr

        self.k0rw_ow = k0rw_ow
        self.k0ro_ow = k0ro_ow
        self.k0ro_go = k0ro_go
        self.k0rg_go = k0rg_go

        self.nw      = nw
        self.no      = no
        self.ng      = ng

    @property
    def slr_go(self):
        """Getter for the residual liquid saturation in gas-oil system."""
        return self.swr+self.sor_go

    def swd(self,sw:np.ndarray):
        """Returns dimensionless water saturation."""
        values = (np.ravel(sw)-self.swr)/(1-self.swr-self.sor_ow)

        return np.clip(values,0,1)

    def sod(self,so:np.ndarray):
        """Returns dimensionless oil saturation."""
        values = (np.ravel(so)-self.sor_ow)/(1-self.swr-self.sor_ow-self.sgr)

        return np.clip(values,0,1)

    def sgd(self,sg:np.ndarray):
        """Returns dimensionless gas saturation."""
        values = (np.ravel(sg)-self.sgr)/(1-self.sgr-self.slr_go)

        return np.clip(values,0,1)

    def som(self,sg:np.ndarray):
        """Returns minimum oil saturation in three-phase system."""
        a = np.ravel(sg)/(1-self.swr-self.sor_go)

        return (1-a)*self.sor_ow+a*self.sor_go

    def sw_star(self,sw,som):
        """Returns dimensionless water saturation (starred)."""
        values = (np.ravel(sw)-self.swr)/(1-self.swr-som)

        return np.clip(values,0,1)

    def so_star(self,so,som):
        """Returns dimensionless oil saturation (starred)."""
        values = (np.ravel(so)-som)/(1-self.swr-som)

        return np.clip(values,0,1)

    def sg_star(self,sg,som):
        """Returns dimensionless gas saturation (starred)."""
        values = (np.ravel(sg))/(1-self.swr-som)

        return np.clip(values,0,1)

    def imbibition(self,sw:np.ndarray,so:np.ndarray,sg:np.ndarray):
        """Computes three-phase relative permeabilities krw, kro and krg.

        Parameters:
        ----------
        sw : water saturation
        so : oil saturation
        sg : gas saturation

        """
        ow = BrooksCorey(self.swr,self.sor_ow,self.k0rw_ow,self.k0ro_ow,self.nw,self.no)
        go = BrooksCorey(self.slr_go,self.sgr,self.k0ro_go,self.k0rg_go,self.no,self.ng)

        krw,kro_ow = ow.imbibition(sw)
        kro_go,krg = go.imbibition(sw+so)

        sod = self.sod(so)
        som = self.som(sg)

        beta_w = (kro_ow)/(1-self.sw_star(sw,som))
        beta_o = (self.k0ro_ow)/self.so_star(so,som)
        beta_g = (kro_go)/(1-self.sg_star(sg,som))

        kro = (beta_w*beta_g)/beta_o

        kro[sod<0] = 0
        kro[sod>1] = self.k0ro_ow

        return krw,kro,krg

if __name__ == "__main__":

    rp = StonesI(swr=0.1,sor_ow=0.4,sor_go=0.2,sgr=0.05,k0rw_ow=0.3,k0ro_ow=0.8,k0ro_go=0.8,k0rg_go=0.3)

    print(rp.imbibition(0.3,0.5,0.2))