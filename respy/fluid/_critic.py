def natural_gas(spgr):
    """It calculates pseudo-critical temperature and pressure for
    natural-gas systems based on specific gravity."""
    tpc = 168+325*spgr-12.5*spgr**2
    ppc = 677+15*spgr-37.5*spgr**2
    return tpc,ppc

def gas_condensate(spgr):
    """It calculates pseudo-critical temperature and pressure for
    gas-condensate systems based on specific gravity."""
    tpc = 1887+33*spgr-71.5*spgr**2
    ppc = 706-51.7*spgr-11.1*spgr**2
    return tpc,ppc

def composition(mfrac,tcrit,pcrit):
    """It calculates pseudo-critical temperature and pressure based on
    mole fraction and pseudo properties of each component."""
    tpcs = [frac*T for frac,T in zip(mfrac,tcrit)]
    ppcs = [frac*P for frac,P in zip(mfrac,tcrit)]
    return sum(tpcs),sum(ppcs)