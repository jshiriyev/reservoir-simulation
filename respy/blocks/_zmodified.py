import numpy

class ZModified():

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