from inspect import isfunction

if __name__ == "__main__":
    from _fluid import Fluid
else:
    from respy.fluid import Fluid

class Gas():
    """
    The hydrocarbon gases that are normally found in a natural gas are
    methanes, ethanes, propanes, butanes, pentanes, and small amounts of
    hexanes and heavier. The nonhydrocarbon gases include carbon dioxide,
    hydrogen sulfide, and nitrogen.
    """

    def __init__(self,spgr:float=None,*,crit=None,zfact=None,visc=None):
        """
        Real gas defined based on specific gravity and critical parameters:
        
        spgr    : specific gravity of the gas at
                  standard conditions.

        crit    : either tuple of (pcrit,tcrit) or method name to calculate (pcrit,tcrit)

        zfact   : either constant value, or method to calculate z-factor
                  from pressure and temperature.

        visc    : either constant value, or method to calculate viscosity
                  from pressure and temperature.
        """

        self._spgr = spgr

        if isinstance(crit,tuple):
            self._pcrit,self._tcrit = crit
        elif isinstance(crit,str):
            self._pcrit,self._tcrit = getattr(self,crit)(spgr)

        self._zfact = zfact

        if isfunction(visc) or visc is None:
            self._visc = visc
        else:
            self._visc = visc*0.001

    @property
    def spgr(self):
        return self._spgr
    
    @property
    def mw(self):
        """Molecular Mass, lb/lbmol"""
        return 28.964*self._spgr

    @property
    def _mw(self):
        """Molecular Mass, kg/mol"""
        return 28.964*self._spgr/1000

    @property
    def pcrit(self):
        """Critical Pressure in psi"""
        return self._pcrit

    @property
    def tcrit(self):
        """Critical Temperature in Rankine"""
        return self._tcrit
    
    @staticmethod
    def natural_gas(spgr):
        """It calculates pseudo-critical temperature and pressure for
        natural-gas systems based on specific gravity."""
        ppc = 677+15*spgr-37.5*spgr**2
        tpc = 168+325*spgr-12.5*spgr**2
        return ppc,tpc

    @staticmethod
    def gas_condensate(spgr):
        """It calculates pseudo-critical temperature and pressure for
        gas-condensate systems based on specific gravity."""
        ppc = 706-51.7*spgr-11.1*spgr**2
        tpc = 187+330.0*spgr-71.5*spgr**2
        return ppc,tpc

    @staticmethod
    def high_mole_weight(spgr):
        """It calculates critical properties for high molecular
        weight reservoir gases, and returns:

        ppc     : Critical Pressure in psia
        tpc     : Critical Temperature in °R
        """
        ppc = 756.8-131*spgr-3.6*spgr**2
        tpc = 169.2+349.5*spgr-74*spgr**2
        return ppc,tpc

    def pred(self,P):
        return P/self.pcrit

    def tred(self,T):
        return T/self.tcrit

    @property
    def zfact(self):
        return self._zfact

    @property
    def zprime(self):
        if isfunction(self._zprime):
            return
        return 0

    @property
    def visc(self):
        if isfunction(self._visc) or self._visc is None:
            return self._visc
        return self._visc/0.001
    
    def __call__(self,P,T):
        """
        P   : Pressure, psia
        T   : Temperature, °R
        """

        Z,Zprime = self.get_zfact(P,T)

        visc = self.get_visc(P,T,Z)
        rho  = self.get_rho(P,T,Z)
        comp = self.get_comp(P,T,Z,Zprime)
        fvf  = self.get_fvf(P,T,Z)

        return Fluid(visc,rho,comp,fvf)

    def get_zfact(self,P,T):

        if not isfunction(self._zfact):
            return self._zfact,0

        Pr = self.pred(P)
        Tr = self.tred(T)
        
        return self._zfact(Pr,Tr,derivative=True)

    def get_visc(self,P,T,Z=1,**kwargs):
        """Calculates pressure and temperature dependent viscosity in cp"""
        if isfunction(self._visc):
            return self._visc(P,T,Z,**kwargs)

        return self._visc

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