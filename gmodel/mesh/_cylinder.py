from matplotlib import pyplot

import numpy

class Ellipse():

    # This class is supposed to create 2-D surface in 3-D domain
    # thickness: thickness of the ellipse
    # rinner: inner radius
    # dip_angle: 

    # lamda: node spacing, radius ratio

    def __init__(self,radii=None,inner_radii=(0,0),thickness=1,**kwargs):

        super().__init__(**kwargs)

        if radii is not None:
            self.set_radii(radii)

        self.set_thickness(thickness)

        self.gridFlag = False

    def set_origin(self,origin=(0,0,0)):

        # origin: location of the center of ellipse

        self.origin = origin # origing has not been implemented yet

    def set_radii(self,radii=(1,1)):

        # lengths: (major_radius,minor_radius)

        if type(radii) == int:
            self.radii = (float(radii),float(radii))
        elif type(radii) == float:
            self.radii = (radii,radii)
        elif type(radii) == str:
            raise ValueError('Radii must be defined as a number or tuple.')
        elif len(radii) == 2:
            self.radii = radii

        self.radiusMajor = max(self.radii)
        self.radiusMinor = min(self.radii)

        self.set_ellipse()

    def set_inradii(self,inradii=(0,0)):

        if type(inradii) == int:
            self.inradii = (float(inradii),float(inradii))
        elif type(inradii) == float:
            self.inradii = (inradii,inradii)
        elif type(inradii) == str:
            raise ValueError('Inner radii must be defined as a number or tuple.')
        elif len(inradii) == 2:
            self.inradii = inradii

        self.inradiusMajor = max(inradii)
        self.inradiusMinor = min(inradii)

    def set_area(self,area,aspect=1):
        
        self.area = area

        radius1 = numpy.sqrt(area*aspect/numpy.pi)
        radius2 = radius1/aspect

        self.set_radii(radii=(radius1,radius2))

    def set_thickness(self,thickness):

        self.thickness = thickness

        if hasattr(self,"radii"):
            self.set_ellipse()

    def set_ellipse(self):

        numverts = 50

        thetas = numpy.linspace(0,2*numpy.pi,numverts+1)[:-1]

        self.edge_vertices = numpy.zeros((2*numverts,3))

        self.edge_vertices[:,0] = numpy.tile(self.radii[0]/2*numpy.cos(thetas),2)+self.radii[0]/2
        self.edge_vertices[:,1] = numpy.tile(self.radii[1]/2*numpy.sin(thetas),2)+self.radii[1]/2
        self.edge_vertices[:,2] = numpy.append(numpy.zeros(numverts),self.thickness*numpy.ones(numverts))

        indices = numpy.empty((2*numverts,2),dtype=int)

        vertices_0 = numpy.arange(numverts)
        vertices_1 = numpy.append(numpy.arange(numverts)[1:],0)

        indices[:,0] = numpy.append(vertices_0,vertices_0+numverts)
        indices[:,1] = numpy.append(vertices_1,vertices_1+numverts)

        x_aspects = self.edge_vertices[:,0][indices]
        y_aspects = self.edge_vertices[:,1][indices]
        z_aspects = self.edge_vertices[:,2][indices]

        self.boundaries = []

        for x_aspect,y_aspect,z_aspect in zip(x_aspects,y_aspects,z_aspects):
            self.boundaries.append(numpy.array([x_aspect,y_aspect,z_aspect]))

        self.gridFlag = False

    def grid(self,lamda):

        self.lamda = lamda

    def plot(self,axis,showVertices=False,showBounds=True,showGridEdges=False,showGridCenters=False):

        if showVertices:
            axis.scatter(*self.edge_vertices.T)

        if showBounds:
            for line in self.boundaries:
                axis.plot(*line,color='grey')

        if showGridEdges:
            for node in self.grid_xnodes[1:-1]:
                axis.vlines(x=node,ymin=0,ymax=self.lengths[1],linestyle="--")
            for node in self.grid_ynodes[1:-1]:
                axis.hlines(y=node,xmin=0,xmax=self.lengths[0],linestyle="--")

        # axis.set_box_aspect(self.radii[1]/self.radii[0])

        if showGridCenters:
            axis.scatter(*self.grid_centers.T)

class Cylinder():

    """
    For cylindrical disk, dimensions is a tuple with two entries for sizes in r,z direction
    """
    
    def __init__(self,lengths,**kwargs):

        super().__init__(**kwargs)

        self.lengths = lengths

        numverts = 50

        thetas = numpy.linspace(0,2*numpy.pi,numverts+1)[:-1]

        self.edge_vertices = numpy.zeros((2*numverts,3))

        self.edge_vertices[:,0] = numpy.tile(self.lengths[0]/2*numpy.cos(thetas),2)+self.lengths[0]/2
        self.edge_vertices[:,1] = numpy.tile(self.lengths[1]/2*numpy.sin(thetas),2)+self.lengths[1]/2
        self.edge_vertices[:,2] = numpy.append(numpy.zeros(numverts),self.lengths[2]*numpy.ones(numverts))

        indices = numpy.empty((2*numverts,2),dtype=int)

        vertices_0 = numpy.arange(numverts)
        vertices_1 = numpy.append(numpy.arange(numverts)[1:],0)

        indices[:,0] = numpy.append(vertices_0,vertices_0+numverts)
        indices[:,1] = numpy.append(vertices_1,vertices_1+numverts)

        x_aspects = self.edge_vertices[:,0][indices]
        y_aspects = self.edge_vertices[:,1][indices]
        z_aspects = self.edge_vertices[:,2][indices]

        self.boundaries = []

        for x_aspect,y_aspect,z_aspect in zip(x_aspects,y_aspects,z_aspects):
            self.boundaries.append(numpy.array([x_aspect,y_aspect,z_aspect]))

    def grid(self):
        pass

    def plot(self):
        pass