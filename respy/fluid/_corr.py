def wichert_aziz(tpc,ppc,H2S=0.0,CO2=0.0):
    """It corrects pseudo-critical temperature and pressure based on
    the non-hydrocarbon mole fraction."""

    A,B = H2S+CO2,H2S
    
    epsilon = 120*(A**0.9-A**1.6)+15*(B**0.5-B**4.0)

    tpc_dash = tpc-epsilon
    ppc_dash = ppc*tpc_dash/(tpc+B*(1-B)*epsilon)

    tpc,ppc = tpc_dash,ppc_dash

    return tpc,ppc

def carr_kobayashi_burrows(tpc,ppc,H2S=0.0,CO2=0.0,N2=0.0):
    """It corrects pseudo-critical temperature and pressure based on
    the non-hydrocarbon mole fraction."""

    tpc = tpc-80*CO2+130*H2S-250*N2
    ppc = ppc+440*CO2+600*H2S-170*N2

    return tpc,ppc

def high_molecular_weight():

    pass