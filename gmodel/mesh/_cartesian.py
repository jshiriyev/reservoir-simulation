import numpy

from dataclasses import dataclass

class Voxel():

    """
    For rectangular parallelepiped, dimensions is a tuple
    with three entries for sizes in x,y,z direction.
    """

    def __init__(self,lengths,**kwargs):

        super().__init__(**kwargs)

        self.lengths = lengths

        self.edge_vertices = numpy.zeros((8,3))

        self.edge_vertices[0,:] = (0,0,0)
        self.edge_vertices[1,:] = (self.lengths[0],0,0)
        self.edge_vertices[2,:] = (self.lengths[0],self.lengths[1],0)
        self.edge_vertices[3,:] = (0,self.lengths[1],0)

        self.edge_vertices[4,:] = (0,0,self.lengths[2])
        self.edge_vertices[5,:] = (self.lengths[0],0,self.lengths[2])
        self.edge_vertices[6,:] = (self.lengths[0],self.lengths[1],self.lengths[2])
        self.edge_vertices[7,:] = (0,self.lengths[1],self.lengths[2])

        indices = numpy.empty((12,2),dtype=int)

        indices[:,0] = (0,1,2,3,0,1,2,3,4,5,6,7)
        indices[:,1] = (1,2,3,0,4,5,6,7,5,6,7,4)
        
        x_aspects = self.edge_vertices[:,0][indices]
        y_aspects = self.edge_vertices[:,1][indices]
        z_aspects = self.edge_vertices[:,2][indices]

        self.boundaries = []

        for x_aspect,y_aspect,z_aspect in zip(x_aspects,y_aspects,z_aspects):
            self.boundaries.append(numpy.array([x_aspect,y_aspect,z_aspect]))

        self.gridFlag = False

    def grid(self,grid_num):

        """
        self.grid_num        : number of grids in x, y, z directions
        self.grid_numtot     : number of totla grids 
        self.grid_indices    : connectivity map containing index of all grids and their neighbours.
        self.grid_sizes      : size of grids in all directions.
        self.grid_areas      : area of all faces
        self.grid_volumes    : volume of grids
        self.grid_centers    : coordinates of the center of grids
        """
        
        self.grid_num = grid_num

        self.grid_numtot = numpy.prod(self.grid_num)

        idx = numpy.arange(self.grid_numtot)
        
        self.grid_indices = numpy.tile(idx,(7,1)).T

        self.grid_indices[idx.reshape(-1,self.grid_num[0])[:,1:].ravel(),1] -= 1
        self.grid_indices[idx.reshape(-1,self.grid_num[0])[:,:-1].ravel(),2] += 1
        self.grid_indices[idx.reshape(self.grid_num[2],-1)[:,self.grid_num[0]:],3] -= self.grid_num[0]
        self.grid_indices[idx.reshape(self.grid_num[2],-1)[:,:-self.grid_num[0]],4] += self.grid_num[0]
        self.grid_indices[idx.reshape(self.grid_num[2],-1)[1:,:],5] -= self.grid_num[0]*self.grid_num[1]
        self.grid_indices[idx.reshape(self.grid_num[2],-1)[:-1,:],6] += self.grid_num[0]*self.grid_num[1]

        self.grid_hasxmin = ~(self.grid_indices[:,0]==self.grid_indices[:,1])
        self.grid_hasxmax = ~(self.grid_indices[:,0]==self.grid_indices[:,2])
        self.grid_hasymin = ~(self.grid_indices[:,0]==self.grid_indices[:,3])
        self.grid_hasymax = ~(self.grid_indices[:,0]==self.grid_indices[:,4])
        self.grid_haszmin = ~(self.grid_indices[:,0]==self.grid_indices[:,5])
        self.grid_haszmax = ~(self.grid_indices[:,0]==self.grid_indices[:,6])

        self.grid_xnodes = numpy.linspace(0,self.lengths[0],self.grid_num[0]+1)
        self.grid_ynodes = numpy.linspace(0,self.lengths[1],self.grid_num[1]+1)
        self.grid_znodes = numpy.linspace(0,self.lengths[2],self.grid_num[2]+1)
        
        xsize = self.grid_xnodes[1:]-self.grid_xnodes[:-1]
        ysize = self.grid_ynodes[1:]-self.grid_ynodes[:-1]
        zsize = self.grid_znodes[1:]-self.grid_znodes[:-1]
        
        self.grid_sizes = numpy.zeros((self.grid_numtot,3))
        self.grid_sizes[:,0] = numpy.tile(xsize,self.grid_num[1]*self.grid_num[2])
        self.grid_sizes[:,1] = numpy.tile(ysize.repeat(self.grid_num[0]),self.grid_num[2])
        self.grid_sizes[:,2] = zsize.repeat(self.grid_num[0]*self.grid_num[1])

        self.grid_areas = numpy.zeros((self.grid_numtot,3))
        self.grid_areas[:,0] = self.grid_sizes[:,1]*self.grid_sizes[:,2]
        self.grid_areas[:,1] = self.grid_sizes[:,2]*self.grid_sizes[:,0]
        self.grid_areas[:,2] = self.grid_sizes[:,0]*self.grid_sizes[:,1]

        self.grid_volumes = numpy.prod(self.grid_sizes,axis=1)

        xcenter = self.grid_xnodes[:-1]+xsize/2
        ycenter = self.grid_ynodes[:-1]+ysize/2
        zcenter = self.grid_znodes[:-1]+zsize/2
        
        self.grid_centers = numpy.zeros((self.grid_numtot,3))
        self.grid_centers[:,0] = numpy.tile(xcenter,self.grid_num[1]*self.grid_num[2])
        self.grid_centers[:,1] = numpy.tile(ycenter.repeat(self.grid_num[0]),self.grid_num[2])
        self.grid_centers[:,2] = zcenter.repeat(self.grid_num[0]*self.grid_num[1])

        self.gridFlag = True

    def plot(self):

        pass

class OneDimCuboid():

    def __init__(self,length,numtot,csa):

        self.length = length

        self.numtot = numtot

        self.csa = csa

        self.num = (numtot,1,1)

        self._index()

        self._size()

        self._area()

    def _index(self):

        idx = numpy.arange(self.numtot)

        self.index = numpy.tile(idx,(7,1)).T

        self.index[idx.reshape(-1,self.num[0])[:,1:].ravel(),1] -= 1
        self.index[idx.reshape(-1,self.num[0])[:,:-1].ravel(),2] += 1
        self.index[idx.reshape(self.num[2],-1)[:,self.num[0]:],3] -= self.num[0]
        self.index[idx.reshape(self.num[2],-1)[:,:-self.num[0]],4] += self.num[0]
        self.index[idx.reshape(self.num[2],-1)[1:,:],5] -= self.num[0]*self.num[1]
        self.index[idx.reshape(self.num[2],-1)[:-1,:],6] += self.num[0]*self.num[1]

    def _size(self):

        self.size = numpy.zeros((self.numtot,3))

        self.size[:,0] = self.length/self.num[0]
        self.size[:,1] = self.csa
        self.size[:,2] = 1

        self.xaxis = numpy.linspace(0,self.length,self.numtot,endpoint=False)+self.size[:,0]/2

    def _area(self):

        self.area = numpy.zeros((self.numtot,6))

        self.area[:,0] = self.size[:,1]*self.size[:,2]
        self.area[:,1] = self.size[:,1]*self.size[:,2]
        self.area[:,2] = self.size[:,2]*self.size[:,0]
        self.area[:,3] = self.size[:,2]*self.size[:,0]
        self.area[:,4] = self.size[:,0]*self.size[:,1]
        self.area[:,5] = self.size[:,0]*self.size[:,1]

    def set_permeability(self,permeability,homogeneous=True,isotropic=True):

        self.permeability = numpy.zeros((self.numtot,3))

        if homogeneous and isotropic:
            
            self.permeability[:] = permeability

        elif homogeneous and not isotropic:

            self.permeability[:] = permeability

        elif not homogeneous and isotropic:

            self.permeability[:,0] = permeability
            self.permeability[:,1] = permeability
            self.permeability[:,2] = permeability

        else:

            self.permeability[:] = permeability

class RectRectGrid():
    """It is a 2D object in 3D space."""

    def __init__(self,length=1,width=1,height=1,centroid=None):
        """Initialization of rectangle in 3D domain.
        If centroid is not defined, left-bottom vertex will be assigned to (0,0) point."""

        super().__setattr__('length',length)
        super().__setattr__('width',width)
        super().__setattr__('height',height)

        if centroid is None:
            centroid = (length/2,width/2)

        super().__setattr__('centroid',numpy.array(centroid))

        self._vertices()

    def _vertices(self):

        vertices = numpy.zeros((4,2),dtype='float64')

        vertices[0,:] = self.centroid+(-self.length/2,-self.width/2)
        vertices[1,:] = self.centroid+(-self.length/2,+self.width/2)
        vertices[2,:] = self.centroid+(+self.length/2,+self.width/2)
        vertices[3,:] = self.centroid+(+self.length/2,-self.width/2)

        super().__setattr__('vertices',vertices)

        indices = numpy.array((0,1,2,3,0),dtype='int32')

        super().__setattr__('indices',indices)

    def mesh(self,number,**kwargs):

        self._initialize(number)

        xdelta = kwargs.get('xdelta')
        ydelta = kwargs.get('ydelta')
        zdelta = kwargs.get('zdelta')

        self._delta(xdelta,ydelta,zdelta)

        self._center()

        self._area()

        self._volume()

    def _initialize(self,number):
        """It initializes grids:

        number: must be a tuple or list (number of x grids, number of y grids)"""

        super().__setattr__('number',number)

        super().__setattr__('count',numpy.prod(self.number))

        self._cmap()

        hasxmin = ~(self.cmap[:,0]==self.cmap[:,1])
        hasxmax = ~(self.cmap[:,0]==self.cmap[:,2])

        hasymin = ~(self.cmap[:,0]==self.cmap[:,3])
        hasymax = ~(self.cmap[:,0]==self.cmap[:,4])

        super().__setattr__('hasxmin',hasxmin)
        super().__setattr__('hasxmax',hasxmax)
        super().__setattr__('hasymin',hasymin)
        super().__setattr__('hasymax',hasymax)

    def _cmap(self):

        indices = numpy.arange(self.count)
        
        cmap = numpy.tile(indices,(5,1)).T

        cmap[indices.reshape(-1,self.xnumber)[:,1:].ravel(),1] -= 1
        cmap[indices.reshape(-1,self.xnumber)[:,:-1].ravel(),2] += 1

        cmap[indices.reshape(1,-1)[:,self.xnumber:],3] -= self.xnumber
        cmap[indices.reshape(1,-1)[:,:-self.xnumber],4] += self.xnumber

        super().__setattr__('cmap',cmap)

    def _delta(self,xdelta=None,ydelta=None,zdelta=None):

        if xdelta is None:
            xdelta = self.length/self.xnumber

        if ydelta is None:
            ydelta = self.width/self.ynumber

        if zdelta is None:
            zdelta = self.height

        delta = numpy.zeros((self.count,3),dtype='float64')

        xdelta = numpy.array(xdelta).flatten()
        ydelta = numpy.array(ydelta).flatten()

        if xdelta.size == 1:
            delta[:,0] = xdelta
        elif xdelta.size == self.xnumber:
            delta[:,0] = numpy.tile(xdelta,self.ynumber)
        else:
            raise ValueError

        if ydelta.size == 1:
            delta[:,1] = ydelta
        elif ydelta.size == self.ynumber:
            delta[:,1] = ydelta.repeat(self.xnumber)
        else:
            raise ValueError

        delta[:,2] = zdelta
        
        super().__setattr__('delta',delta)

    def _center(self):

        xdelta = self.xdelta[:self.xnumber]
        ydelta = self.ydelta[::self.xnumber]

        xcenter = numpy.cumsum(xdelta)-xdelta/2
        ycenter = numpy.cumsum(ydelta)-ydelta/2

        zcenter = self.height/2

        center = numpy.zeros((self.count,3),dtype='float64')

        xcenter = numpy.array(xcenter).flatten()
        ycenter = numpy.array(ycenter).flatten()

        if xcenter.size==1:
            center[:,0] = xcenter
        elif xcenter.size==self.xnumber:
            center[:,0] = numpy.tile(xcenter,self.ynumber)
        else:
            raise ValueError("xcenter")

        if ycenter.size==1:
            center[:,1] = ycenter
        elif ycenter.size == self.ynumber:
            center[:,1] = ycenter.repeat(self.xnumber)
        else:
            raise ValueError("ycenter")

        center[:,2] = zcenter

        super().__setattr__('center',center)

    def _area(self):

        area = numpy.zeros((self.count,3))

        area[:,0] = self.zdelta*self.ydelta
        area[:,1] = self.xdelta*self.zdelta
        area[:,2] = self.ydelta*self.xdelta

        super().__setattr__('area',area)

    def _volume(self):

        volume = numpy.prod(self.delta,axis=1)

        super().__setattr__('volume',volume)

    def view(self,axis=None,vertices=True,bounds=True,edges=True,centers=True):

        show = True if axis is None else False

        if axis is None:
            axis = pyplot.figure().add_subplot()

        if vertices:
            axis.scatter(*self.numvertices.T)

        if bounds:
            axis.plot(*self.numvertices[self.indices].T,color='grey')

        if edges:
            xdelta = self.xdelta[:self.xnumber]
            ydelta = self.ydelta[::self.xnumber]
            
            xinner = numpy.cumsum(xdelta)[:-1]
            yinner = numpy.cumsum(ydelta)[:-1]

            axis.vlines(x=xinner,ymin=0,ymax=self.width,linestyle="--")
            axis.hlines(y=yinner,xmin=0,xmax=self.length,linestyle="--")

        if centers:
            axis.scatter(*self.center.T[:2,:])

        axis.set_box_aspect(self.width/self.length)

        if show:
            pyplot.show()

    def __setattr__(self,key,value):

        raise AttributeError(f'RectGrid does not have attribute {key}')

    @property
    def xnumber(self):

        return self.number[0]

    @property
    def ynumber(self):

        return self.number[1]

    @property
    def znumber(self):

        return 1
    
    @property
    def xcenter(self):

        return self.center[:,0]

    @property
    def ycenter(self):

        return self.center[:,1]

    @property
    def zcenter(self):

        return self.center[:,2]

    @property
    def xdelta(self):

        return self.delta[:,0]

    @property
    def ydelta(self):

        return self.delta[:,1]

    @property
    def zdelta(self):

        return self.delta[:,2]

class Hexahedron():

    def __init__(self,origin=(0,0,0),dimens=(1,1,1),azimuth=0.,dx=10.,dy=10.,dz=1.0,zcorn=None):
        """
        Initializes Hexahedral grid model with:
        
        origin  : (Ox,Oy,Oz)
        dimens  : (Nx,Ny,Nz)

        azimuth : azimuth of the collective grid orientation
                  0 degree is East direction, 90 degree is North direction
        """

        self.origin = numpy.array(origin)
        self.dimens = numpy.array(dimens,dtype="int32")

        self.azimuth = azimuth

        self.set_dx(dx)
        self.set_dy(dy)

        if zcorn is None:
            self.set_dz(dz)
        else:
            self.set_zcorn(zcorn)

        self.scalars = {}

    def set_dx(self,dx):

        Nx = self.dimens[0]

        dx = numpy.array(dx).flatten()

        if dx.size==1:
            self.dx = numpy.repeat(dx,Nx)
        elif dx.size==Nx:
            self.dx = dx
        else:
            raise ValueError()

    def set_dy(self,dy):

        Ny = self.dimens[1]

        dy = numpy.array(dy).flatten()

        if dy.size==1:
            self.dy = numpy.repeat(dy,Ny)
        elif dy.size==Ny:
            self.dy = dy
        else:
            raise ValueError()

    def set_dz(self,dz):

        Nz = self.dimens[2]

        dz = numpy.array(dz).flatten()

        if dz.size==1:
            self.dz = numpy.repeat(dz,Nz)
        elif dz.size==Nz:
            self.dz = dz
        else:
            raise ValueError()

        zcorn = numpy.zeros((Nz+1,))+self.origin[2]

        zcorn[1:] += numpy.cumsum(self.dz)

        self.set_zcorn(zcorn.repeat(self.numverticesxy))

    def set_zcorn(self,zcorn=None):
        """The size of zcorn must be equal to (Nx+1)*(Ny+1)*(Nz+1)."""

        if zcorn.size==self.numvertices:
            self.zcorn = zcorn
        else:
            raise ValueError()

    def set_scalars(self,**kwargs):

        if len(kwargs)==0:
            raise AttributeError
        elif len(kwargs)>1:
            raise AttributeError

        key,array = kwargs.popitem()

        self.scalars[key] = array

    @property
    def numgrids(self):
        """Returns number of total grids."""
        return numpy.prod(self.dimens)

    @property
    def numgridsxy(self):
        """Returns number of total grids in one layer."""
        return numpy.prod(self.dimens[:2])

    @property
    def numvertices(self):
        """Returns number of total vertices."""
        return numpy.prod(self.dimens+1)

    @property
    def numverticesxy(self):
        """Returns number of total vertices on the surface."""
        return numpy.prod(self.dimens[:2]+1)
    
    @property
    def vertices(self):
        """Returns ID of nodes on the surface."""
        Nx,Ny,Nz = self.dimens

        surface_nodes = numpy.zeros((self.numgridsxy,4),dtype="int32")

        indices = numpy.arange(self.numverticesxy)

        indices = indices.reshape((Ny+1,Nx+1))

        surface_nodes[:,0] = indices[:Ny,:Nx].flatten()
        surface_nodes[:,1] = indices[:Ny,1:].flatten()
        surface_nodes[:,2] = indices[1:,1:].flatten()
        surface_nodes[:,3] = indices[1:,:Nx].flatten()

        return surface_nodes

    @property
    def connectivity(self):
        """Returns connectivity map of the grids"""
        
        indices = numpy.arange(self.numgrids)

        Nx,Ny,Nz = self.dimens
        
        connectivity_map = numpy.tile(indices,(6,1)).T

        connectivity_map[indices.reshape(-1,Nx)[:,1:].ravel(),0] -= 1
        connectivity_map[indices.reshape(-1,Nx)[:,:-1].ravel(),1] += 1

        connectivity_map[indices.reshape(Nz,-1)[:,Nx:],2] -= Nx
        connectivity_map[indices.reshape(Nz,-1)[:,:-Nx],3] += Nx

        connectivity_map[indices.reshape(Nz,-1)[1:,:],4] -= Nx*Ny
        connectivity_map[indices.reshape(Nz,-1)[:-1,:],5] += Nx*Ny

        return connectivity_map

    @property
    def coords(self):
        
        Nx,Ny,Nz = self.dimens

        coordinates = numpy.zeros((self.numvertices,3))

        xcoord = numpy.zeros((Nx+1,))+self.origin[0]
        ycoord = numpy.zeros((Ny+1,))+self.origin[1]

        xcoord[1:] += numpy.cumsum(self.dx)
        ycoord[1:] += numpy.cumsum(self.dy)

        coordinates[:,0] = numpy.tile(xcoord,(Ny+1)*(Nz+1))
        coordinates[:,1] = numpy.tile(numpy.repeat(ycoord,Nx+1),(Nz+1))
        coordinates[:,2] = self.zcorn

        return coordinates

    @property
    def coordsxy(self):

        Nx,Ny,_ = self.dimens

        coordinates = numpy.zeros((self.numverticesxy,3))

        xcoord = numpy.zeros((Nx+1,))+self.origin[0]
        ycoord = numpy.zeros((Ny+1,))+self.origin[1]

        xcoord[1:] += numpy.cumsum(self.dx)
        ycoord[1:] += numpy.cumsum(self.dy)

        coordinates[:,0] = numpy.tile(xcoord,(Ny+1))
        coordinates[:,1] = numpy.repeat(ycoord,Nx+1)

        return coordinates
    
    @property
    def tops(self):
        """Returns top of surface grids, (Ngrids @surface,)."""
        return self.zcorn[self.vertices].mean(axis=1)

    @property
    def centers(self):
        """Returns center of all grids, (Ngrids,3)."""
        Nx,Ny,Nz = self.dimens

        grid_centers = numpy.zeros((self.numgrids,3))

        xcoords = self.coordsxy[:,0][self.vertices]
        ycoords = self.coordsxy[:,1][self.vertices]

        xcenter = xcoords.mean(axis=1)
        ycenter = ycoords.mean(axis=1)

        zcoords = self.coords[:,2].reshape((-1,self.numverticesxy))

        layers = [zcoords[k][self.vertices].mean(axis=1) for k in range(Nz+1)]

        zcenters = [(layers[k]+layers[k+1])/2 for k in range(Nz)]

        grid_centers[:,0] = numpy.tile(xcenter,Nz)
        grid_centers[:,1] = numpy.tile(ycenter,Nz)
        grid_centers[:,2] = numpy.array(zcenters).flatten()

        return grid_centers

    @property
    def areas(self):
        """Returns area of all grids, (Ngrids,6)"""
        pass

    @property
    def volumes(self):
        """Returns volume of all grids, (Ngrids,)."""
        pass

    def vtk(self,filename):

        with open(filename,"w") as vtkfile:

            vtkfile.write("# vtk DataFile Version 3.0\n")
            vtkfile.write("3D scalar data\n")
            vtkfile.write("ASCII\n")
            vtkfile.write("DATASET UNSTRUCTURED_GRID\n")

            vtkfile.write(f"POINTS {self.numvertices} float\n")

            numpy.savetxt(vtkfile,self.coords,fmt="%.1f")

            vtkfile.write("\n")

            vtkfile.write(f"CELLS {self.numgrids} {self.numgrids*9}\n")

            vertices = self.vertices

            for k in range(self.dimens[2]):

                array = numpy.concatenate((
                    vertices+k*self.numverticesxy,
                    vertices+(k+1)*self.numverticesxy),
                    axis=1)

                types = numpy.array([8],dtype='int32').repeat(self.numgridsxy)

                types = types.reshape((-1,1))

                array = numpy.concatenate((types,array),axis=1)

                numpy.savetxt(vtkfile,array,fmt="%d")

            vtkfile.write("\n")
            vtkfile.write(f"CELL_TYPES {self.numgrids}\n")

            types = numpy.array([12],dtype='int32').repeat(self.numgrids)

            numpy.savetxt(vtkfile,types,fmt="%d")

            vtkfile.write("\n")

            vtkfile.write(f"CELL_DATA {self.numgrids}\n")

            for key,array in self.scalars.items():

                dtype = ''.join([i for i in str(array.dtype) if not i.isdigit()])

                vtkfile.write(f"SCALARS {key} {dtype} 1")

                vtkfile.write("\n")
                vtkfile.write("LOOKUP_TABLE default\n")

                numpy.savetxt(vtkfile,array,fmt="%.1f")

    def grdecl(self,filename):

        pass

@dataclass
class Scalar:
    """It is a scalar data type dictionary."""
    name: str
    array: numpy.ndarray
    vtkdatatype: str = None
