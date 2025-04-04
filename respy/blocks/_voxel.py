import numpy as np

class Voxel():
    """
    For rectangular parallelepiped, dimensions is a tuple
    with three entries for sizes in x,y,z direction.
    """

    def __init__(self,dims):

        self._dims = dims

        edge = np.empty((12,2),dtype=int)

        edge[:,0] = (0,1,2,3,0,1,2,3,4,5,6,7)
        edge[:,1] = (1,2,3,0,4,5,6,7,5,6,7,4)

        self._edge = edge

    @property
    def dims(self):
        return self._dims
    
    @property
    def vertex(self):

        vertices = np.zeros((8,3))

        x,y,z = self._dims

        vertices[0,:] = (0,0,0)
        vertices[1,:] = (x,0,0)
        vertices[2,:] = (x,y,0)
        vertices[3,:] = (0,y,0)
        vertices[4,:] = (0,0,z)
        vertices[5,:] = (x,0,z)
        vertices[6,:] = (x,y,z)
        vertices[7,:] = (0,y,z)

        return vertices

    @property
    def edge(self):
        return self._edge

    @property
    def edgex(self):
        return self.vertex[:,0][self.edge]

    @property
    def edgey(self):
        return self.vertex[:,1][self.edge]

    @property
    def edgez(self):
        return self.vertex[:,2][self.edge]

if __name__ == "__main__":

    vox = Voxel((3,3,1))

    print(vox.edge)

    # print(vox.boundaries)