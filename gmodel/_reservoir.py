from ._items import Segment
from ._faults import Faults
from ._fracnet import Fractures
from ._stock import Stock

class Reservoir():
    """It is a geometry to grid and assign petrophysical properties.
    It is a collection of segments"""

    def __init__(self):

        self.segments = []

    def add_segment(self,**kwargs):

        self.segments.append(Segment(kwargs))

    def add_faultsystem(self,**kwargs):

        self.faults = FaultSystem(kwargs)

    def add_fracturenetwork(self,**kwargs):

        self.fractures = FractureNetwork(kwargs)

    def add_stock(self,**kwargs):

        self.stock = Stock(**kwargs)

    def set_porosity(self,porosity,homogeneous=True,X=None,Y=None,Z=None):

        if isinstance(porosity,float):
            porosity = numpy.array([porosity])
        elif isinstance(porosity,int):
            porosity = numpy.array([porosity])
        elif isinstance(porosity,tuple):
            porosity = numpy.array(porosity)

        if homogeneous and self.gridFlag:
            self.porosity = numpy.ones(self.grid_numtot)*porosity
        else:
            self.porosity = porosity

    def set_permeability(self,permeability,isotropy=True,homogeneous=True,X=None,Y=None,Z=None):

        if isinstance(permeability,float):
            permeability = numpy.array([permeability])
        elif isinstance(permeability,int):
            permeability = numpy.array([permeability])
        elif isinstance(permeability,tuple):
            permeability = numpy.array(permeability)

        if homogeneous and isotropy:

            if self.gridFlag:
                self.permeability = numpy.ones((self.grid_numtot,3))*permeability
            else:
                self.permeability = permeability
            
        elif not homogeneous and isotropy:
            
            self.permeability = numpy.repeat(permeability,3).reshape((-1,3))
            
        elif homogeneous and not isotropy:
            
            if self.gridFlag:
                self.permeability = numpy.repeat(permeability,self.grid_numtot).reshape((3,-1)).T
            else:
                self.permeability = permeability
            
        elif not homogeneous and not isotropy:
            
            self.permeability = permeability

    def set_compressibility(self,compressibility):

        self.compressibility = compressibility

    def vtkwrite(res,frac,well,time,sol):
        
        pass