from matplotlib import pyplot

import numpy

class Line():

    def __init__(self,**kwargs):

        super().__init__(**kwargs)

    def set_tracks(self,tracks):
        
        self.tracks = numpy.array(tracks)

    def set2Dview(self,axis,view='zy'):

        if view == 'yx':
            axis.plot(self.tracks[:,0],self.tracks[:,1])
            axis.set_xlabel("x-axis")
            axis.set_ylabel("y-axis")

        elif view == 'zy':
            axis.plot(self.tracks[:,1],self.tracks[:,2])
            axis.set_xlabel("y-axis")
            axis.set_ylabel("z-axis")
            axis.invert_yaxis()

        elif view == 'xz':
            axis.plot(self.tracks[:,2],self.tracks[:,0])
            axis.set_xlabel("z-axis")
            axis.set_ylabel("x-axis")

        elif view == 'xy':
            axis.plot(self.tracks[:,1],self.tracks[:,0])
            axis.set_xlabel("y-axis")
            axis.set_ylabel("x-axis")

        elif view == 'yz':
            axis.plot(self.tracks[:,2],self.tracks[:,1])
            axis.set_xlabel("z-axis")
            axis.set_ylabel("y-axis")

        elif view == 'zx':
            axis.plot(self.tracks[:,0],self.tracks[:,2])
            axis.set_xlabel("x-axis")
            axis.set_ylabel("z-axis")
            axis.invert_yaxis()

    def set3Dview(self,axis):

        axis.plot3D(*self.tracks.T)
        axis.set_xlabel("x-axis")
        axis.set_ylabel("y-axis")
        axis.set_zlabel("z-axis")
        axis.invert_zaxis()