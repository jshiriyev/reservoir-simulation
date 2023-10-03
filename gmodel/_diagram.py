from ._survey import Survey

class Diagram(Survey):
    """It is a well diagram"""

    def __init__(self,**kwargs):

        super().__init__(**kwargs)

    def set_indiameter(self,*args,csaFlag=True):

        """
        indiameter  : Inner Diameter
        csa         : Cross Sectional Are
        H_diameter  : Hydraulic Diameter; 4*Hydraulic_Radius
        H_radius    : Hydraulic Radius; the ratio of the cross-sectional area of
                      a channel or pipe in which a fluid is flowing to the 
                      wetted perimeter of the conduit.
        """
        
        [self.indiameter.append(indiameter) for indiameter in args]

        if csaFlag:
            [self.csa.append(np.pi*indiameter**2/4) for indiameter in args]

        [self.H_diameter.append(indiameter) for indiameter in args]
        [self.H_radius.append(indiameter/4) for indiameter in args]

    def set_roughness(self,*args):

        [self.roughness.append(roughness) for roughness in args]

        [self.roughness_R.append(arg/indiameter) for (arg,indiameter) in zip(args,self.indiameter)]

    def set_nodes(self,zloc=None,elevation=[0,0]):

        """
        Nodes are the locations where the measurements are available, and
        coordinates are selected in such a way that:
        - r-axis shows radial direction
        - theta-axis shows angular direction 
        - z-axis shows lengthwise direction
        """

        if zloc is None:
            self.zloc = [0,self.length]

        self.elevation = elevation

    def plot(self):

        pass