import numpy

class CornerPoint():

	def __init__(self,grdecl):

		self.SPECGRID = grdecl['SPECGRID']

		self.COORD = grdecl['COORD']
		self.ZCORN = grdecl['ZCORN']

		self.ACTNUM = grdecl['ACTNUM'].astype(bool)

	@property
	def dimensions(self):

		return tuple([int(x) for x in self.SPECGRID[0].split()[:3]])

	@property
	def ncells(self):
		return numpy.prod(self.dimensions)

	@property
	def pillars2d(self):
		
		Nx,Ny,_ = self.dimensions

		pillars = numpy.zeros((Nx*Ny,4),dtype=int)

		column0 = numpy.arange(Nx,dtype=int)
		column0 = numpy.tile(column0,Ny)

		column0 += numpy.arange(Ny,dtype=int).repeat(Nx)*(Nx+1)

		pillars[:,0] = column0
		pillars[:,1] = pillars[:,0]+1
		pillars[:,2] = pillars[:,1]+Nx
		pillars[:,3] = pillars[:,2]+1

		return pillars

	@property
	def pillars3d(self):

		_,_,Nz = self.dimensions
		
		pillars = self.pillars2d

		return numpy.tile(pillars,(Nz,1))

	@property
	def xcorn(self):

		Ncells = self.ncells

		zcorn = self.zcorn

		pillars = self.pillars3d

		XCORN = numpy.zeros((Ncells,8))

		for cell in range(Ncells):

			cell_zs = zcorn[cell] #8 values

			cell_ps = pillars[cell] #4 values

			cell_xs = numpy.zeros((8,))

			cell_xs[0::4] = self.getxcorn(cell_zs[0::4],self.COORD[cell_ps[0]])
			cell_xs[1::4] = self.getxcorn(cell_zs[1::4],self.COORD[cell_ps[1]])
			cell_xs[2::4] = self.getxcorn(cell_zs[2::4],self.COORD[cell_ps[2]])
			cell_xs[3::4] = self.getxcorn(cell_zs[3::4],self.COORD[cell_ps[3]])

			XCORN[cell,:] = cell_xs

		return XCORN

	@property
	def ycorn(self):

		Ncells = self.ncells

		zcorn = self.zcorn

		pillars = self.pillars3d

		YCORN = numpy.zeros((Ncells,8))

		for cell in range(Ncells):

			cell_zs = zcorn[cell] #8 values

			cell_ps = pillars[cell] #4 values

			cell_ys = numpy.zeros((8,))

			cell_ys[0::4] = self.getycorn(cell_zs[0::4],self.COORD[cell_ps[0]])
			cell_ys[1::4] = self.getycorn(cell_zs[1::4],self.COORD[cell_ps[1]])
			cell_ys[2::4] = self.getycorn(cell_zs[2::4],self.COORD[cell_ps[2]])
			cell_ys[3::4] = self.getycorn(cell_zs[3::4],self.COORD[cell_ps[3]])

			YCORN[cell,:] = cell_ys

		return YCORN

	@property
	def zcorn(self):

		Nx,Ny,Nz = self.dimensions

		indices = numpy.zeros((Nx*Ny,8),dtype=int)

		column0 = numpy.arange(Nx*2,step=2,dtype=int)
		column0 = numpy.tile(column0,Ny)

		column0 += numpy.arange(Ny,dtype=int).repeat(Nx)*Nx*4

		indices[:,0] = column0
		indices[:,1] = indices[:,0]+1
		indices[:,2] = indices[:,0]+Nx*2
		indices[:,3] = indices[:,2]+1
		indices[:,4] = indices[:,0]+Nx*Ny*4
		indices[:,5] = indices[:,4]+1
		indices[:,6] = indices[:,2]+Nx*Ny*4
		indices[:,7] = indices[:,6]+1

		indices = numpy.tile(indices,(Nz,1))

		kvalues = numpy.arange(Nz,dtype=int).repeat(Nx*Ny)*Nx*Ny*8

		indices += kvalues.reshape((-1,1))

		return self.ZCORN[indices]

	def hexahedron(self):

		Ncells = self.ncells
		
		XCORN = self.xcorn
		YCORN = self.ycorn
		ZCORN = self.zcorn

		p1 = numpy.array((XCORN[:,0],YCORN[:,0],ZCORN[:,0])).T
		p2 = numpy.array((XCORN[:,1],YCORN[:,1],ZCORN[:,1])).T
		p3 = numpy.array((XCORN[:,3],YCORN[:,3],ZCORN[:,3])).T
		p4 = numpy.array((XCORN[:,2],YCORN[:,2],ZCORN[:,2])).T
		p5 = numpy.array((XCORN[:,4],YCORN[:,4],ZCORN[:,4])).T
		p6 = numpy.array((XCORN[:,5],YCORN[:,5],ZCORN[:,5])).T
		p7 = numpy.array((XCORN[:,7],YCORN[:,7],ZCORN[:,7])).T
		p8 = numpy.array((XCORN[:,6],YCORN[:,6],ZCORN[:,6])).T

		return HexaHedron(p1,p2,p3,p4,p5,p6,p7,p8)
	
	@staticmethod
	def getxcorn(zcorn,pillar):

		x1,_,z1,x2,_,z2 = pillar

		if z2==z1:
			return x1

		return x1+(x2-x1)/(z2-z1)*(zcorn-z1)

	@staticmethod
	def getycorn(zcorn,pillar):

		_,y1,z1,_,y2,z2 = pillar

		if z2==z1:
			return y1

		return y1+(y2-y1)/(z2-z1)*(zcorn-z1)

class HexaHedron():

	def __init__(self,r1,r2,r3,r4,r5,r6,r7,r8):

		self.r1 = numpy.array(r1)
		self.r2 = numpy.array(r2)
		self.r3 = numpy.array(r3)
		self.r4 = numpy.array(r4)
		self.r5 = numpy.array(r5)
		self.r6 = numpy.array(r6)
		self.r7 = numpy.array(r7)
		self.r8 = numpy.array(r8)

	@property
	def normal1(self):
		v1 = self.r4-self.r1
		v2 = self.r2-self.r1
		return numpy.cross(v1,v2)

	@property
	def normal2(self):
		v1 = self.r2-self.r1
		v2 = self.r5-self.r1
		return numpy.cross(v1,v2)

	@property
	def normal3(self):
		v1 = self.r3-self.r2
		v2 = self.r6-self.r2
		return numpy.cross(v1,v2)

	@property
	def normal4(self):
		v1 = self.r8-self.r4
		v2 = self.r3-self.r4
		return numpy.cross(v1,v2)

	@property
	def normal5(self):
		v1 = self.r5-self.r1
		v2 = self.r4-self.r1
		return numpy.cross(v1,v2)

	@property
	def normal6(self):
		v1 = self.r6-self.r5
		v2 = self.r8-self.r5
		return numpy.cross(v1,v2)

	def unormal(self,face):

		normal = getattr(self,f"normal{face}")

		magnitude = numpy.linalg.norm(normal,2,axis=1)

		magnitude[magnitude==0] = 1

		return normal/magnitude.reshape((-1,1))

	@property
	def center(self):

		return 1/8*(
			self.r1+self.r2+self.r3+self.r4+
			self.r5+self.r6+self.r7+self.r8)

	@property
	def center1(self):
		return 1/4*(self.r1+self.r2+self.r3+self.r4)

	@property
	def center2(self):
		return 1/4*(self.r1+self.r2+self.r5+self.r6)

	@property
	def center3(self):
		return 1/4*(self.r2+self.r3+self.r6+self.r7)

	@property
	def center4(self):
		return 1/4*(self.r3+self.r4+self.r7+self.r8)

	@property
	def center5(self):
		return 1/4*(self.r1+self.r4+self.r5+self.r8)

	@property
	def center6(self):
		return 1/4*(self.r5+self.r6+self.r7+self.r8)

	@property
	def area1(self):
		
		A1 = self.areaTriangle(self.r1,self.r4,self.r2)
		A2 = self.areaTriangle(self.r3,self.r2,self.r4)
		
		return A1+A2

	@property
	def area2(self):

		A1 = self.areaTriangle(self.r1,self.r2,self.r5)
		A2 = self.areaTriangle(self.r6,self.r5,self.r2)
		
		return A1+A2

	@property
	def area3(self):

		A1 = self.areaTriangle(self.r2,self.r3,self.r6)
		A2 = self.areaTriangle(self.r7,self.r6,self.r3)
		
		return A1+A2

	@property
	def area4(self):
		
		A1 = self.areaTriangle(self.r4,self.r8,self.r3)
		A2 = self.areaTriangle(self.r7,self.r3,self.r8)
		
		return A1+A2

	@property
	def area5(self):
		
		A1 = self.areaTriangle(self.r4,self.r1,self.r8)
		A2 = self.areaTriangle(self.r5,self.r8,self.r1)
		
		return A1+A2

	@property
	def area6(self):
		
		A1 = self.areaTriangle(self.r5,self.r6,self.r8)
		A2 = self.areaTriangle(self.r7,self.r8,self.r6)
		
		return A1+A2
	
	@property
	def volume(self):

		center = self.center

		v1 = numpy.sum((self.center1-center)*self.unormal(1),axis=1)*self.area1
		v2 = numpy.sum((self.center2-center)*self.unormal(2),axis=1)*self.area2
		v3 = numpy.sum((self.center3-center)*self.unormal(3),axis=1)*self.area3
		v4 = numpy.sum((self.center4-center)*self.unormal(4),axis=1)*self.area4
		v5 = numpy.sum((self.center5-center)*self.unormal(5),axis=1)*self.area5
		v6 = numpy.sum((self.center6-center)*self.unormal(6),axis=1)*self.area6

		return 1/3*(v1+v2+v3+v4+v5+v6)

	@staticmethod
	def areaTriangle(r1,r2,r3):

		area = numpy.zeros((r1.shape[0],))

		AB,AC = r2-r1,r3-r1

		ABmag = numpy.linalg.norm(AB,2,axis=1)
		ACmag = numpy.linalg.norm(AC,2,axis=1)

		lengths = ABmag*ACmag

		scalar = numpy.sum(AB*AC,axis=1)

		cos_theta = scalar[lengths!=0]/lengths[lengths!=0]

		sin_theta = numpy.sqrt(1-cos_theta**2)

		area[lengths!=0] = 1/2*lengths[lengths!=0]*sin_theta

		return area

class RecCuboid():

	def __init__(self,xdelta:tuple=None,ydelta:tuple=None,zdelta:tuple=None,length:float=None,width:float=None,height:float=None,num:tuple=None,flodim:int=None):
		"""Three-dimensional reservoir model can be initialized in two different ways:
		
		Method 1:

		xdelta	: length of grids in ft, (Nlength,)
		ydelta	: width of grids in ft, (Nwidth,)
		zdelta	: height of grids in ft, (Nheight,)

		Method 2:

		length 	: length of reservoir in ft, (x direction)
		width 	: width of reservoir in ft, (y direction)
		height 	: height of reservoir in ft, (z direction)

		num 	: number of grids, (Nlength, Nwidth, Nheight)

		Other Inputs:

		flodim 	: flow dimension it can be 1, 2, or 3

		The object has methods to calculate the control volume implementation parameters.

		"""

		if isinstance(xdelta,tuple): # Method 1:
			self.init1(xdelta,ydelta,zdelta)
		else:						 # Method 2:
			self.init2(legnth,width,height,num)

		# The parameters starting with underscore are defined in SI units.
		# The same parameters without underscore are in Oil Field units.

		if flodim is not None:
			self.flodim = flodim
		elif self._num[2]>1:
			self.flodim = 3
		elif self._num[1]>1:
			self.flodim = 2
		else:
			self.flodim = 1

		self.set_gplat()

		self.set_dims()
		self.set_area()
		self.set_volume()

	def init1(self,xdelta,ydelta,zdelta):
		"""Initialization with the first method."""

		self._xdelta = numpy.asarray(xdelta).astype(numpy.float_)*0.3048
		self._ydelta = numpy.asarray(ydelta).astype(numpy.float_)*0.3048
		self._zdelta = numpy.asarray(zdelta).astype(numpy.float_)*0.3048

		self._length = self._xdelta.sum()
		self._width  = self._ydelta.sum()
		self._height = self._zdelta.sum()

		self._num 	 = (len(xdelta),len(ydelta),len(zdelta))

	def init2(self,length,width,height,num):
		"""Initialization with the second method."""

		self._length = length*0.3048
		self._width  = width*0.3048
		self._height = height*0.3048

		self._xdelta = numpy.repeat(self._length/num[0],num[0])
		self._ydelta = numpy.repeat(self._width/num[1],num[1])
		self._zdelta = numpy.repeat(self._height/num[2],num[2])

		self._num 	 = num

	def set_gplat(self):

		num = self._num

		idx = numpy.arange(self.numtot,dtype=numpy.int_)

		self.gplat = numpy.tile(idx,(1+self.flodim*2,1)).T

		self.gplat[idx.reshape(-1,num[0])[:,1:].ravel(),1] -= 1
		self.gplat[idx.reshape(-1,num[0])[:,:-1].ravel(),2] += 1

		if self.flodim>1:
			self.gplat[idx.reshape(num[2],-1)[:,num[0]:],3] -= num[0]
			self.gplat[idx.reshape(num[2],-1)[:,:-num[0]],4] += num[0]

		if self.flodim>2:
			self.gplat[idx.reshape(num[2],-1)[1:,:],5] -= num[0]*num[1]
			self.gplat[idx.reshape(num[2],-1)[:-1,:],6] += num[0]*num[1]

	def set_dims(self):
		"""should be depreciated later"""

		self._xdims = numpy.tile(
			self._xdelta,self._num[1]*self._num[2]).reshape((-1,1))

		self._ydims = numpy.tile(
			numpy.repeat(self._ydelta,self._num[0]),self._num[2]).reshape((-1,1))

		self._zdims = numpy.repeat(
			self._zdelta,self._num[0]*self._num[1]).reshape((-1,1))

	def set_area(self):

		self._xarea = self._ydims*self._zdims
		self._yarea = self._zdims*self._xdims
		self._zarea = self._xdims*self._ydims

	def set_volume(self):

		self._volume = self._xdims*self._ydims*self._zdims

	def get_property(self,quality,conversion_factor=1.,dtype=None):

		quality = numpy.asarray(quality)

		if dtype is not None:
			quality = quality.astype(dtype)

		quality = quality.flatten()*conversion_factor

		if quality.size==1:
			quality = quality.repeat(self.numtot)

		return quality.reshape((-1,1))

	def set_depth(self,depth):
		"""Assigns the depth values in ft to the grids."""
		self._depth = self.get_property(depth,conversion_factor=0.3048,dtype=numpy.float_)

	def set_poro(self,poro):
		"""Assigns the porosity values in fractions to the grids."""

		self._poro = self.get_property(poro,dtype=numpy.float_)

	def set_perm(self,xperm,yperm=None,zperm=None,yreduce=1.,zreduce=1.):
		"""Assigns the permeability values in mD to the grids."""

		self._xperm = self.get_property(
			xperm,conversion_factor=9.869233e-16,dtype=numpy.float_)

		self._yperm = self._xperm*yreduce if yperm is None else self.get_property(
			yperm,conversion_factor=9.869233e-16,dtype=numpy.float_)

		self._zperm = self._xperm*zreduce if zperm is None else self.get_property(
				zperm,conversion_factor=9.869233e-16,dtype=numpy.float_)

	@property
	def num(self):
		return self._num

	@property
	def numtot(self):
		return numpy.prod(self._num).item()

	@property
	def xdelta(self):
		return self._xdelta/0.3048

	@property
	def ydelta(self):
		return self._ydelta/0.3048

	@property
	def zdelta(self):
		return self._zdelta/0.3048
	
	@property
	def length(self):
		return self._length/0.3048

	@property
	def width(self):
		return self._width/0.3048

	@property
	def height(self):
		return self._height/0.3048

	@property
	def _xaxis(self):

		side2 = numpy.cumsum(self._xdelta)
		side1 = numpy.roll(side2,1)

		side1[0] = 0

		return (side1+side2)/2

	@property
	def _yaxis(self):

		side2 = numpy.cumsum(self._ydelta)
		side1 = numpy.roll(side2,1)

		side1[0] = 0

		return (side1+side2)/2

	@property
	def _zaxis(self):

		side2 = numpy.cumsum(self._zdelta)
		side1 = numpy.roll(side2,1)

		side1[0] = 0

		return (side1+side2)/2

	@property
	def xaxis(self):
		return self._xaxis/0.3048

	@property
	def yaxis(self):
		return self._yaxis/0.3048

	@property
	def zaxis(self):
		return self._zaxis/0.3048

	@property
	def xdims(self):
		return self._xdims/0.3048

	@property
	def ydims(self):
		return self._ydims/0.3048

	@property
	def zdims(self):
		return self._zdims/0.3048

	@property
	def xarea(self):
		return self._xarea/0.3048**2

	@property
	def yarea(self):
		return self._yarea/0.3048**2

	@property
	def zarea(self):
		return self._zarea/0.3048**2

	@property
	def volume(self):
		return self._volume/0.3048**3

	@property
	def depth(self):
		return self._depth/0.3048
	
	@property
	def poro(self):
		return self._poro

	@property
	def xperm(self):
		return self._xperm/9.869233e-16

	@property
	def yperm(self):
		return self._yperm/9.869233e-16

	@property
	def zperm(self):
		return self._zperm/9.869233e-16

	@property
	def xmin(self):
		"""Properties of grids on x-minimum boundary"""
		return RecCuboidGrid(self,1,edge=True)

	@property
	def xpos(self):
		"""Properties of grids that has x-positive neighbors"""
		return RecCuboidGrid(self,2)

	@property
	def xneg(self):
		"""Properties of grids that has x-negative neighbors"""
		return RecCuboidGrid(self,1)

	@property
	def xmax(self):
		"""Properties of grids on x-maximum boundary"""
		return RecCuboidGrid(self,2,edge=True)

	@property
	def ymin(self):
		"""Properties of grids on y-minimum boundary"""
		if self.flodim>1:
			return RecCuboidGrid(self,3,edge=True)

	@property
	def ypos(self):
		"""Properties of grids that has y-positive neighbors"""
		if self.flodim>1:
			return RecCuboidGrid(self,4)

	@property
	def yneg(self):
		"""Properties of grids that has y-negative neighbors"""
		if self.flodim>1:
			return RecCuboidGrid(self,3)

	@property
	def ymax(self):
		"""Properties of grids on y-maximum boundary"""
		if self.flodim>1:
			return RecCuboidGrid(self,4,edge=True)

	@property
	def zmin(self):
		"""Properties of grids on z-minimum boundary"""
		if self.flodim>2:
			return RecCuboidGrid(self,5,edge=True)

	@property
	def zpos(self):
		"""Properties of grids that has z-positive neighbors"""
		if self.flodim>2:
			return RecCuboidGrid(self,6)

	@property
	def zneg(self):
		"""Properties of grids that has z-negative neighbors"""
		if self.flodim>2:
			return RecCuboidGrid(self,5)

	@property
	def zmax(self):
		"""Properties of grids on z-maximum boundary"""
		if self.flodim>2:
			return RecCuboidGrid(self,6,edge=True)

class RecCuboidGrid(numpy.ndarray):

	def __new__(cls,grid,path,edge=False):
		"""
		grid 	: RecCuboid instance
		path 	: direction (1: west, 2: east, 3: south, 4: north, 5: down, 6: up)
		edge 	: boundary (True) or inner (False)
		"""

		item = (grid.gplat[:,0]==grid.gplat[:,path])

		if not edge:
			item = ~item

		obj = numpy.asarray(grid.gplat[item,0],dtype=numpy.int_).view(cls)

		obj.grid = grid

		obj.axis = int((path-1)/2)

		return obj

	def __array_finalize__(self,obj):

		if obj is None: return

		self.grid = getattr(obj,'grid',None)
		self.axis = getattr(obj,'axis',None)

	def __getattr__(self,key):

		("dims","area","volume","depth","poro","perm")

		unit = key[1:] if key.startswith('_') else key	

		if unit in ("volume","depth","poro"):
			return getattr(self.grid,key)[self,0]

		if unit[0] not in ("x","y","z"):
			return

		if unit[1:] in ("dims","area","perm"):
			return getattr(self.grid,key)[self,0]		

if __name__ == "__main__":

	# cells = Hexahedron(
	# 	((-5,-5,5),(5,-5,5)),
	# 	((5,-5,5),(15,-5,5)),
	# 	((5,-5,-5),(15,-5,-5)),
	# 	((-5,-5,-5),(5,-5,-5)),
	# 	((-5,5,5),(5,5,5)),
	# 	((5,5,5),(15,5,5)),
	# 	((5,5,-5),(15,5,-5)),
	# 	((-5,5,-5),(5,5,-5)),
	# )

	# print(cells.unormal1)

	# print(cells.center)

	# print(cells.center1)
	# print(cells.center2)
	# print(cells.center3)
	# print(cells.center4)
	# print(cells.center5)
	# print(cells.center6)

	# print(cells.volume)

	# grids = RecCuboid((4,1,1))

	# print(grids.gplat)
	# print(grids.area)
	# print(grids.zaxis)
	# print(grids.volume)

	grid = RecCuboid((750,1000,1250),(750,1000,1250),(20,))

	# print(grid.size)

	# print(grid.size_test)

	# print(grid.size_test.ypos)

	print(grid.xneg1._area)

	# print(area.xmin)