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

	def __init__(self,num:tuple,length:float=None,width:float=None,height:float=None,dimension:int=None):
		"""Three-dimensional reservoir model defined by:
		
		num 	: number of grids, (Nlength, Nwidth, Nheight)

		length 	: length of reservoir in ft, (x direction
		width 	: width of reservoir in ft, (y direction)
		height 	: height of reservoir in ft, (z direction)

		dimension 	: it can be 1, 2, or 3

		to calculate the control volume implementation parameters.

		"""

		self._num 	 = num

		self._length = num[0]*1000*0.3048 if length is None else length*0.3048
		self._width  = num[1]*1000*0.3048 if width is None else width*0.3048
		self._height = num[2]*1000*0.3048 if height is None else height*0.3048

		# The parameters starting with underscore are defined in SI units.
		# The same parameters without underscore are in Oil Field units.

		if dimension is not None:
			pass
		elif num[2]>1:
			dimension = 3
		elif num[1]>1:
			dimension = 2
		else:
			dimension = 1

		self.set_index(dimension)
		self.set_size()
		self.set_area(dimension)
		self.set_volume()

	def set_index(self,dimension:int=3):
		"""dimension can be 1, 2, or 3"""

		idx = numpy.arange(self.numtot)

		num = self._num

		self.index = numpy.tile(idx,(1+dimension*2,1)).T

		self.index[idx.reshape(-1,num[0])[:,1:].ravel(),1] -= 1
		self.index[idx.reshape(-1,num[0])[:,:-1].ravel(),2] += 1

		if dimension>1:
			self.index[idx.reshape(num[2],-1)[:,num[0]:],3] -= num[0]
			self.index[idx.reshape(num[2],-1)[:,:-num[0]],4] += num[0]

		if dimension>2:
			self.index[idx.reshape(num[2],-1)[1:,:],5] -= num[0]*num[1]
			self.index[idx.reshape(num[2],-1)[:-1,:],6] += num[0]*num[1]

	def set_xaxis(self):

		self._xaxis = numpy.arange(
			self._length/self._num[0]/2,
			self._length,
			self._length/self._num[0])

	def set_yaxis(self):

		self._yaxis = numpy.arange(
			self._width/self._num[1]/2,
			self._width,
			self._width/self._num[1])

	def set_zaxis(self):

		self._zaxis = numpy.arange(
			self._height/self._num[2]/2,
			self._height,
			self._height/self._num[2])

	def set_size(self):

		self._size = numpy.zeros((self.numtot,3))

		self._size[:,0] = self._length/self._num[0]
		self._size[:,1] = self._width/self._num[1]
		self._size[:,2] = self._height/self._num[2]

	def set_area(self,dimension:int=3):
		"""dimension can be 1, 2, or 3"""

		self._area = numpy.zeros((self.numtot,dimension*2))

		self._area[:,0] = self._size[:,1]*self._size[:,2]
		self._area[:,1] = self._size[:,1]*self._size[:,2]

		if dimension>1:
			self._area[:,2] = self._size[:,2]*self._size[:,0]
			self._area[:,3] = self._size[:,2]*self._size[:,0]

		if dimension>2:
			self._area[:,4] = self._size[:,0]*self._size[:,1]
			self._area[:,5] = self._size[:,0]*self._size[:,1]

	def set_volume(self):

		self._volume = numpy.prod(self._size,axis=1).reshape((-1,1))

	def set_permeability(self,permeability):
		"""permeability in mD"""

		self._perm = numpy.asarray(permeability).flatten()
		self._perm *= 9.869233e-16

		if self._perm.size==1:
			self._perm = self._perm.repeat(self._numtot)

	def set_porosity(self,porosity):
		"""porosity in fractions"""

		self._poro = numpy.asarray(porosity).flatten()

		if self._poro.size==1:
			self._poro = self._poro.repeat(self._numtot)

	@property
	def num(self):
		return self._num

	@property
	def numtot(self):
		return numpy.prod(self._num).item()

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