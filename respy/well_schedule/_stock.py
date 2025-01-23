import numpy

# from .directory._browser import Browser

class WellStock():

    def __iter__(self):

        wells = zip(self.Wells.Trajectory.tracks,
            self.Wells.consbhp,
            self.Wells.limits,
            )

        for index,(track,flag,limit) in enumerate(wells):
            pass

    def get_wconds(self):

        for index,(track,flag,limit) in enumerate(wells):

            ttrack = np.transpose(track[:,:,np.newaxis],(2,1,0))

            vector = self.PorRock.grid_centers[:,:,np.newaxis]-ttrack

            distance = np.sqrt(np.sum(vector**2,axis=1))

            wgrid = np.unique(np.argmin(distance,axis=0))

            well_grids = np.append(self.well_grids,wgrid)
            well_indices = np.append(self.well_indices,np.full(wgrid.size,index,dtype=int))
            well_bhpflags = np.append(self.well_bhpflags,np.full(wgrid.size,flag,dtype=bool))
            well_limits = np.append(self.well_limits,np.full(wgrid.size,limit,dtype=float))

        return well_grids, well_indices, well_bhpflags, well_limits