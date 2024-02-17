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

	def __init__(self,
		xdelta:tuple=None,
		ydelta:tuple=None,
		zdelta:tuple=None,
		length:float=None,
		width:float=None,
		height:float=None,
		num:tuple=None,
		dim:int=None):
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

		dim 	: flow dimension it can be 1, 2, or 3

		The object has methods to calculate the control volume implementation parameters.

		"""

		if isinstance(xdelta,tuple): # Method 1:
			self.init1(xdelta,ydelta,zdelta)
		else:						 # Method 2:
			self.init2(legnth,width,height,num)

		# The parameters starting with underscore are defined in SI units.
		# The same parameters without underscore are in Oil Field units.

		if dim is not None:
			self.dimension = dim
		elif self._num[2]>1:
			self.dimension = 3
		elif self._num[1]>1:
			self.dimension = 2
		else:
			self.dimension = 1

		self.set_index()
		self.set_size()
		self.set_area()
		self.set_volume()

	def init1(self,xdelta,ydelta,zdelta):
		"""Initialization with the first method."""

		self._xdelta = numpy.asarray(xdelta).astype(numpy.float64)*0.3048
		self._ydelta = numpy.asarray(ydelta).astype(numpy.float64)*0.3048
		self._zdelta = numpy.asarray(zdelta).astype(numpy.float64)*0.3048

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

	def set_index(self):

		idx = numpy.arange(self.numtot)

		num = self._num

		self.index = numpy.tile(idx,(1+self.dimension*2,1)).T

		self.index[idx.reshape(-1,num[0])[:,1:].ravel(),1] -= 1
		self.index[idx.reshape(-1,num[0])[:,:-1].ravel(),2] += 1

		if self.dimension>1:
			self.index[idx.reshape(num[2],-1)[:,num[0]:],3] -= num[0]
			self.index[idx.reshape(num[2],-1)[:,:-num[0]],4] += num[0]

		if self.dimension>2:
			self.index[idx.reshape(num[2],-1)[1:,:],5] -= num[0]*num[1]
			self.index[idx.reshape(num[2],-1)[:-1,:],6] += num[0]*num[1]

	def set_xaxis(self):

		side2 = numpy.cumsum(self._xdelta)
		side1 = numpy.roll(side2,1)

		side1[0] = 0

		self._xaxis = (side1+side2)/2

		# self._xaxis = numpy.arange(
		# 	self._length/self._num[0]/2,
		# 	self._length,
		# 	self._length/self._num[0])

	def set_yaxis(self):

		side2 = numpy.cumsum(self._ydelta)
		side1 = numpy.roll(side2,1)

		side1[0] = 0

		self._yaxis = (side1+side2)/2

		# self._yaxis = numpy.arange(
		# 	self._width/self._num[1]/2,
		# 	self._width,
		# 	self._width/self._num[1])

	def set_zaxis(self):

		side2 = numpy.cumsum(self._zdelta)
		side1 = numpy.roll(side2,1)

		side1[0] = 0

		self._zaxis = (side1+side2)/2

		# self._zaxis = numpy.arange(
		# 	self._height/self._num[2]/2,
		# 	self._height,
		# 	self._height/self._num[2])

	def set_size(self):

		self._size = numpy.zeros((self.numtot,3))

		# self._size[:,0] = self._length/self._num[0]
		# self._size[:,1] = self._width/self._num[1]
		# self._size[:,2] = self._height/self._num[2]

		self._size[:,0] = numpy.tile(self._xdelta,self._num[1]*self._num[2])

		self._size[:,1] = numpy.tile(
			numpy.repeat(self._ydelta,self._num[0]),self._num[2])

		self._size[:,2] = numpy.repeat(self._zdelta,self._num[0]*self._num[1])

	def set_area(self):

		self._area = numpy.zeros((self.numtot,self.dimension*2))

		self._area[:,0] = self._size[:,1]*self._size[:,2]
		self._area[:,1] = self._size[:,1]*self._size[:,2]

		if self.dimension>1:
			self._area[:,2] = self._size[:,2]*self._size[:,0]
			self._area[:,3] = self._size[:,2]*self._size[:,0]

		if self.dimension>2:
			self._area[:,4] = self._size[:,0]*self._size[:,1]
			self._area[:,5] = self._size[:,0]*self._size[:,1]

	def set_volume(self):

		self._volume = numpy.prod(self._size,axis=1).reshape((-1,1))

	def get_property(self,quality,conversion_factor=1.,dtype=None):

		quality = numpy.asarray(quality)

		if dtype is not None:
			quality = quality.astype(dtype)

		quality = quality.flatten()*conversion_factor

		if quality.size==1:
			quality = quality.repeat(self._numtot)

		return quality.reshape((-1,1))

	def set_porosity(self,porosity):
		"""porosity in fractions"""

		self._poro = self.get_property(porosity,dtype=numpy.float64)

	def set_permeability(self,xperm,yperm=None,zperm=None,yreduce=1.,zreduce=1.):
		"""xperm in mD"""

		_xperm = self.get_property(
			xperm,conversion_factor=9.869233e-16,dtype=numpy.float64)

		if yperm is not None:
			_yperm = self.get_property(
				yperm,conversion_factor=9.869233e-16,dtype=numpy.float64)
		else:
			_yperm = _xperm*yreduce

		if zperm is not None:
			_zperm = self.get_property(
				zperm,conversion_factor=9.869233e-16,dtype=numpy.float64)
		else:
			_zperm = _xperm*zreduce

		self._perm = numpy.concatenate((_xperm,_yperm,_zperm),axis=1)

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
	def xaxis(self):
		if not hasattr(self,"_xaxis"):
			self.set_xaxis()
		return self._xaxis/0.3048

	@property
	def yaxis(self):
		if not hasattr(self,"_yaxis"):
			self.set_yaxis()
		return self._yaxis/0.3048

	@property
	def zaxis(self):
		if not hasattr(self,"_zaxis"):
			self.set_zaxis()
		return self._zaxis/0.3048

	@property
	def size(self):
		return self._size/0.3048

	@property
	def area(self):
		return self._area/0.3048**2

	@property
	def volume(self):
		return self._volume/0.3048**3
	
	@property
	def perm(self):
		return self._perm/9.869233e-16

	@property
	def poro(self):
		return self._poro
	
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

	grids = RecCuboid((4,1,1))

	print(grids.index)
	# print(grids.area)
	# print(grids.zaxis)
	# print(grids.volume)