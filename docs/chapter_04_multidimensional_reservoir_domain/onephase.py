import numpy as np

class GridBase():

	def __init__(self,xdelta,ydelta,zdelta):
		"""Initialize grid deltas in feet. Internally stored in meters."""
		self.xdelta = xdelta # ft
		self.ydelta = ydelta # ft
		self.zdelta = zdelta # ft
		
	@property
	def xdelta(self):
		return self._xdelta/0.3048

	@xdelta.setter
	def xdelta(self,value):
		self._xdelta = np.asarray(value).flatten().astype(np.float64)*0.3048

	@property
	def ydelta(self):
		return self._ydelta/0.3048

	@ydelta.setter
	def ydelta(self,value):
		self._ydelta = np.asarray(value).flatten().astype(np.float64)*0.3048

	@property
	def zdelta(self):
		return self._zdelta/0.3048

	@zdelta.setter
	def zdelta(self,value):
		self._zdelta = np.asarray(value).flatten().astype(np.float64)*0.3048

class Grid(GridBase):

    def __init__(self,plat,xdelta,ydelta,zdelta):        
        super().__init__(xdelta,ydelta,zdelta)
        self.plat = plat
        self.volume = None

    @property
    def nums(self):
        """Returns total number of grids."""
        return self.plat.shape[0]

    @property
    def dims(self):
        return int(self.plat.shape[1]/2)

    @property
    def index(self):
        return np.arange(self.nums)

    @property
    def volume(self):
        return self._volume*35.3147

    @volume.setter
    def volume(self,value):
        self._volume = np.prod((self._xdelta,self._ydelta,self._zdelta),axis=0)

class GridDelta(GridBase):
	"""Gridding based on delta values."""
	def __init__(self,xdelta,ydelta,zdelta):
		super().__init__(xdelta,ydelta,zdelta)	

	@property
	def nums(self):
		"""Returns tuple of (xnums,ynums,znums)."""
		return (self.xnums,self.ynums,self.znums)

	@property
	def xnums(self):
		return self._xdelta.size

	@property
	def ynums(self):
		return self._ydelta.size

	@property
	def znums(self):
		return self._zdelta.size

	@property
	def dims(self) -> int:
		"""Returns the number of dimensions based on self.nums."""
		return next((i+1 for i in range(2,-1,-1) if self.nums[i]>1),1)

	@property
	def index(self):
		return np.arange(np.prod(self.nums))

	@property
	def grid(self):

		xdelta = np.tile(self.xdelta,self.ynums*self.znums)
		ydelta = np.tile(np.repeat(self.ydelta,self.xnums),self.znums)
		zdelta = np.repeat(self.zdelta,self.xnums*self.ynums)

		return Grid(self.plat,xdelta,ydelta,zdelta)

	@property
	def plat(self):
		"""Plat of grids that locates neighbor index information."""
		map_ = np.tile(self.index,(self.dims*2,1)).T

		map_[self.index.reshape(-1,self.xnums)[:,1:].ravel(),0] -= 1
		map_[self.index.reshape(-1,self.xnums)[:,:-1].ravel(),1] += 1

		if self.dims>1:
			map_[self.index.reshape(self.znums,-1)[:,self.xnums:],2] -= self.xnums
			map_[self.index.reshape(self.znums,-1)[:,:-self.xnums],3] += self.xnums

		if self.dims>2:
			map_[self.index.reshape(self.znums,-1)[1:,:],4] -= self.xnums*self.ynums
			map_[self.index.reshape(self.znums,-1)[:-1,:],5] += self.xnums*self.ynums

		return map_

class PorMed():
    
    def __init__(self,xperm,poro=None,**kwargs):

        self.set_perm(xperm,**kwargs)
        
        self.poro = poro # fraction
        
    @property
    def perm(self):
        return self._perm/9.869233e-16
        
    def set_perm(self,xperm,*,yperm=None,zperm=None,yreduce:float=1.,zreduce:float=10.):
        """Assigns the permeability values in mD to the grids.
    
        xperm 	: permeability in x-direction, mD
        yperm   : permeability in y direction, mD
        zperm   : permeability in z direction, mD
    
        yreduce : xperm to yperm ratio, dimensionless
        zreduce : xperm to zperm ratio, dimensionless
    
        """
        self.xperm = xperm
        self.yperm = (yperm,yreduce)
        self.zperm = (zperm,zreduce)

        self._perm = np.column_stack((self._xperm,self._yperm,self._zperm))

    @property
    def xperm(self):
        return self._xperm/9.869233e-16

    @xperm.setter
    def xperm(self,value):
        self._xperm = value*9.869233e-16

    @property
    def yperm(self):
        return self._yperm/9.869233e-16

    @yperm.setter
    def yperm(self,value):
        self._yperm = self._xperm/value[1] if value[0] is None else value[0]*9.869233e-16
    
    @property
    def zperm(self):
        return self._zperm/9.869233e-16

    @zperm.setter
    def zperm(self,value):
        self._zperm = self._xperm/value[1] if value[0] is None else value[0]*9.869233e-16

    @property
    def poro(self):
        return self._poro
        
    @poro.setter
    def poro(self,value):
        self._poro = value

class Fluid():

	def __init__(self,visc):
		self.visc = visc # cp

	@property
	def visc(self):
		return self._visc*1000

	@visc.setter
	def visc(self,value):
		self._visc = value/1000

class Time():
	"""Time class for the reservoir simulation"""
	def __init__(self,delta:float,total:float):
		"""Defining the time settings for the simulator

		delta   : first time step defined in days
		total   : total simulation time defined in days
		"""

		self.delta = delta
		self.total = total

		self.times = None
		self.steps = None

	@property
	def delta(self):
		return self._delta/(24*60*60)

	@delta.setter
	def delta(self,value):
		self._delta = value*24*60*60

	@property
	def total(self):
		return self._total/(24*60*60)

	@total.setter
	def total(self,value):
		self._total = value*24*60*60

	@property
	def times(self):
		return self._times/(24*60*60)

	@times.setter
	def times(self,value):
		"""Sets linearly spaced time data."""
		self._times = np.arange(
			self._total+self._delta/2,step=self._delta)

	@property
	def steps(self):
		return self._steps/(24*60*60)

	@steps.setter
	def steps(self,value):
		"""Sets the time steps in the given time array."""
		self._steps = self._times[1:]-self._times[:-1]

	def __iter__(self):
		"""It starts from 0 time and iterates till the last time step."""

		zipped = zip(self._times,self._steps)

		for index,(time,step) in enumerate(zipped):

			yield index,time,step

	@property
	def nums(self):
		return self._steps.size

class Matrix():

    def __init__(self,T:np.ndarray,J:np.ndarray,A:np.ndarray,Q:np.ndarray,G:np.ndarray):
        """
        T   : Inter-block transmissibility, matrix
        J   : Productivity index, diagonal matrix
        A   : Block accumulation, diagonal matrix
        Q   : Source term, vector
        G   : Gravity term, vector
        """
        self.T = T # transmissibility
        self.J = J # productivity
        self.A = A # block accumulation
        self.Q = Q # source
        self.G = G # gravity

    @property
    def T(self):
        """Converting from SI Units to Oil Field Units."""
        return self._T*(3.28084**3)*(24*60*60)*6894.76

    @T.setter
    def T(self,value):
        self._T = value/((3.28084**3)*(24*60*60)*6894.76)

    @property
    def J(self):
        """Converting from SI Units to Oil Field Units."""
        return self._J*(3.28084**3)*(24*60*60)*6894.76

    @J.setter
    def J(self,value):
        self._J = value/((3.28084**3)*(24*60*60)*6894.76)
        
    @property
    def A(self):
        """Converting from SI Units to Oil Field Units."""
        return self._A*(3.28084**3)*(24*60*60)

    @A.setter
    def A(self,value):
        self._A = value/((3.28084**3)*(24*60*60))

    @property
    def Q(self):
        """Converting from SI Units to Oil Field Units."""
        return self._Q*(3.28084**3)*(24*60*60)

    @Q.setter
    def Q(self,value):
        self._Q = value/((3.28084**3)*(24*60*60))

    @property
    def G(self):
        """Converting from SI Units to Oil Field Units."""
        return self._G*(3.28084**3)*(24*60*60)

    @G.setter
    def G(self,value):
        self._G = value/((3.28084**3)*(24*60*60))

class Block():

    def __init__(self,grid,res,fluid,*args):

        self.grid = grid
        self.res = res
        self.fluid = fluid
        self.wells = args

    def transmissibility(self):
        self._Tx = self.res._xperm*self.grid._ydelta*self.grid._zdelta/self.grid._xdelta/self.fluid._visc
        self._Ty = self.res._yperm*self.grid._zdelta*self.grid._xdelta/self.grid._ydelta/self.fluid._visc
        self._Tz = self.res._zperm*self.grid._xdelta*self.grid._ydelta/self.grid._zdelta/self.fluid._visc

    @property
    def Tx(self):
        return self._Tx*(3.28084**3)*(24*60*60)*6894.76
    
    @property
    def Ty(self):
        return self._Ty*(3.28084**3)*(24*60*60)*6894.76
    
    @property
    def Tz(self):
        return self._Tz*(3.28084**3)*(24*60*60)*6894.76

    def accumulation(self,tdelta:float):
        self._A = (self.grid._volume*self.res._poro)/(tdelta*86400)
        
    @property
    def A(self):
        return self._A*(3.28084**3)*(24*60*60)

    def productivity(self):

        self._J = []

        for well in self.wells:
            
            if well.axis=="x":
                k1 = self.res._yperm[list(well.block)]
                k2 = self.res._zperm[list(well.block)]
                w1 = self.grid._ydelta[list(well.block)]
                w2 = self.grid._zdelta[list(well.block)]
                w3 = self.grid._xdelta[list(well.block)]
            elif well.axis=='y':
                k1 = self.res._xperm[list(well.block)]
                k2 = self.res._zperm[list(well.block)]
                w1 = self.grid._xdelta[list(well.block)]
                w2 = self.grid._zdelta[list(well.block)]
                w3 = self.grid._ydelta[list(well.block)]
            elif well.axis=='z':
                k1 = self.res._xperm[list(well.block)]
                k2 = self.res._yperm[list(well.block)]
                w1 = self.grid._xdelta[list(well.block)]
                w2 = self.grid._ydelta[list(well.block)]
                w3 = self.grid._zdelta[list(well.block)]
    
            req = self.requivalent(k1,k2,w1,w2)
        
            self._J.append(2*np.pi*w3*np.sqrt(k1*k2)/(self.fluid._visc*(np.log(req/well._radius)+well.skin)))

    @property
    def J(self):
        return [j*(3.28084**3)*(24*60*60)*6894.76 for j in self._J]

    @staticmethod
    def requivalent(perm1,perm2,edge1,edge2):
        """Returns equivalent radius of grids that contains well."""
        sqrt21 = np.power(perm2/perm1,1/2)*np.power(edge1,2)
        sqrt12 = np.power(perm1/perm2,1/2)*np.power(edge2,2)

        quar21 = np.power(perm2/perm1,1/4)
        quar12 = np.power(perm1/perm2,1/4)

        return 0.28*np.sqrt(sqrt21+sqrt12)/(quar21+quar12)

class Build():

    def __init__(self,plat):
        self.plat = plat

    @property
    def nums(self):
        return self.plat.shape[0]

    @property
    def dims(self):
        return int(self.plat.shape[1]/2)

    @property
    def index(self):
        return np.arange(self.nums)

    def T(self,Tx,Ty,Tz):

        index = self.index
        
        T = np.zeros((self.nums,self.nums))

        Sxneg = index[self.plat[:,0]!=index]
        Nxneg = self.plat[self.plat[:,0]!=index,0]

        txneg = Build.harmonic_mean(Tx[Sxneg],Tx[Nxneg])
        
        T[Sxneg,Nxneg] -= txneg
        T[Sxneg,Sxneg] += txneg
        
        Sxpos = index[self.plat[:,1]!=index]
        Nxpos = self.plat[self.plat[:,1]!=index,1]

        txpos = Build.harmonic_mean(Tx[Sxpos],Tx[Nxpos])
        
        T[Sxpos,Nxpos] -= txpos
        T[Sxpos,Sxpos] += txpos
        
        Syneg = index[self.plat[:,2]!=index]
        Nyneg = self.plat[self.plat[:,2]!=index,2]

        tyneg = Build.harmonic_mean(Ty[Syneg],Ty[Nyneg])

        T[Syneg,Nyneg] -= tyneg
        T[Syneg,Syneg] += tyneg
        
        Sypos = index[self.plat[:,3]!=index]
        Nypos = self.plat[self.plat[:,3]!=index,3]

        typos = Build.harmonic_mean(Ty[Sypos],Ty[Nypos])

        T[Sypos,Nypos] -= typos
        T[Sypos,Sypos] += typos

        return T

    def A(self,accumulation):

        A = np.zeros((self.nums,self.nums))

        A[self.index,self.index] = accumulation

        return A

    @staticmethod
    def harmonic_mean(term1,term2):
        """Returns harmonic mean of two terms, term1 and term 2."""
        return (2*term1*term2)/(term1+term2)