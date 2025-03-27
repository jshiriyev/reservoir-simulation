import numpy as np

def k_model(self):

    N = 1000

    self.Sw = np.linspace(self.Swr,1-self.Sor,N)

    self.kro = 2*(1-self.Sw-self.Sor)**2
    self.krw = (self.Sw-self.Swr)**3

def coreymodel(self,koro,korw,m,n):

    N = 1000

    self.Sw = np.linspace(self.Swr,1-self.Sor,N)

    S = (self.Sw-self.Swr)/(1-self.Swr-self.Sor)

    self.kro = koro*(1-S)**m
    self.krw = korw*S**n

    ## end-point mobility ratio calculation
    self.Mo = (korw/self.muw)/(koro/self.muo)