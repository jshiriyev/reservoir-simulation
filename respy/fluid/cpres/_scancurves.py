import numpy

def ScanCurves(sw,hyst,cappres,epspc=0,residual=0):
    """Initialization of a reservoir is assumed to occur by primary drainage (oleic
    displaces aqueous phase) by migration from a source rock.

    During waterflooding, the capillary pressure transitions from a drainage curve to
    an imbibition curve.

    Capillary scanning curves are utilized in simulators to model this
    transition.

    A capillary scanning curve is a weighted average of the imbibition
    and drainage capillary pressure curves."""

    conds = sw<=hyst

    cappres = numpy.zeros(sw.shape)

    drainage = cappres.drainage(sw)

    imbibition = cappres.imbibition(sw)

    Av = 1-residual-hyst
    Bv = sw-hyst

    fv = (Av+epspc)/(Bv+epspc)*Bv/Av

    cappres[conds] = drainage

    cappres[~conds] = fv*imbibition+(1-fv)*drainage

    return cappres

if __name__ == "__main__":

    cappres = BrooksCorey(0.1,0.4,lamda=2,entry=3.5)

    pcDm = cappres.drainage(0.3)
    pcIm = cappres.imbibition(0.3)

    print(pcDm,pcIm)

    jfunctionD = JFunction.direct(pcDm,100,0.2,20,0)
    jfunctionI = JFunction.direct(pcIm,100,0.2,20,0)

    pcDc = JFunction.inverse(jfunctionD,50,0.15,25,0)
    pcIc = JFunction.inverse(jfunctionI,50,0.15,25,0)

    print(pcDc,pcIc)