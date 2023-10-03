import os
import sys

import numpy as np

from scipy.stats import norm

# from geomodel.connectivity import SpatProp

class kriging():

    def __init__(self,obsSpatProp,estSpatProp):
        
        # super(kriging,self).__init__(None,**kwargs)
        
        self.obs = obsSpatProp
        self.est = estSpatProp

    def simple(self,mean=None):

        if mean is None:
            self.mean = self.obs.mean()
        else:
            self.mean = mean
        
        "perc -> percentile, perc=0.5 gives mean values"

        self.distance = SpatProp.get_distance(self.est,self.obs)
        
        _,self.covariance = SpatProp.get_varmodel(
            self.distance,self.obs.type,
            self.obs.sill,self.obs.range,
            self.obs.nugget)

        self.lambdas = np.linalg.solve(self.obs.covariance,self.covariance)
        
        calc_diff_arr = self.obs.reshape((-1,1))-self.mean
        calc_prop_mat = self.lambdas*(calc_diff_arr)
        calc_vars_mat = self.lambdas*self.covariance
        
        self.property = self.mean+calc_prop_mat.sum(axis=0)
        self.variance = self.obs.sill-calc_vars_mat.sum(axis=0)

    def ordinary(self):
        
        "perc -> percentile, perc=0.5 gives mean values"

        self.distance = SpatProp.get_distance(self.est,self.obs)
        
        _,self.covariance = SpatProp.get_varmodel(
            self.distance,self.obs.type,
            self.obs.sill,self.obs.range,
            self.obs.nugget)
        
        Am = self.var.covariance
        Ar = np.ones(Am.shape[0]).reshape((-1,1))
        Ab = np.ones(Am.shape[0]).reshape((1,-1))
        Ab = np.append(Ab,np.array([[0]]),axis=1)
        Am = np.append(Am,Ar,axis=1)
        Am = np.append(Am,Ab,axis=0)

        bm = self.covariance
        bb = np.ones(bm.shape[1]).reshape((1,-1))
        bm = np.append(bm,bb,axis=0)

        xm = np.linalg.solve(Am,bm)
        
        self.lambdas = xm[:-1,:]
        self.beta = xm[-1,:]

        calc_prop_mat = self.lambdas*self.obs.reshape((-1,1))
        calc_vars_mat = self.lambdas*self.covariance

        self.property = calc_prop_mat.sum(axis=0)
        self.variance = self.sill-self.beta-calc_vars_mat.sum(axis=0)

    def get_percentile(self,perc=0.5):

        return self.property+norm.ppf(perc)*np.sqrt(self.variance)

    def gaussian_simulation(self):

        perc = np.random.rand(self.x.size)

        self.get_percentile(perc=perc)