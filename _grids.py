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

class OneDimRecCuboid():

	def __init__(self,length,csarea,numtot):
		"""One dimensional reservoir model defined by:

		length 	: length of reservoir in ft
		csarea 	: cross sectional area perpendicular to the flow in ft2
		numtot 	: number of grids in the direction of flow

		"""

		self._length = length*0.3048

		self._csarea = csarea*0.3048**2

		# The parameters starting with underscore are defined in SI units.
		# The same parameters without underscore are in Oil Field units.

		self._numtot = numtot

		self._num = (numtot,1,1)

		# The following parameters calculated for the Finite Difference implementation.

		self.set_index()

		self.set_xaxis()

		self.set_size()

		self.set_area()

	def set_index(self):

		idx = numpy.arange(self.numtot)

		self.index = numpy.tile(idx,(3,1)).T

		self.index[idx.reshape(-1,self.numtot)[:,1:].ravel(),1] -= 1
		self.index[idx.reshape(-1,self.numtot)[:,:-1].ravel(),2] += 1
		# self.index[idx.reshape(1,-1)[:,self.numtot:],3] -= self.numtot
		# self.index[idx.reshape(1,-1)[:,:-self.numtot],4] += self.numtot
		# self.index[idx.reshape(1,-1)[1:,:],5] -= self.numtot
		# self.index[idx.reshape(1,-1)[:-1,:],6] += self.numtot

	def set_xaxis(self):

		self._xaxis = numpy.arange(
			self._length/self._numtot/2,
			self._length,
			self._length/self._numtot)

	def set_size(self):

		self._size = numpy.zeros((self.numtot,1))

		self._size[:,0] = self._length/self._numtot
		# self._size[:,1] = self._area
		# self._size[:,2] = 1

	def set_area(self):

		self._area = numpy.zeros((self.numtot,2))

		self._area[:,0] = self._csarea # self._size[:,1]*self._size[:,2]
		self._area[:,1] = self._csarea # self._size[:,1]*self._size[:,2]
		# self._area[:,2] = self._size[:,2]*self._size[:,0]
		# self._area[:,3] = self._size[:,2]*self._size[:,0]
		# self._area[:,4] = self._size[:,0]*self._size[:,1]
		# self._area[:,5] = self._size[:,0]*self._size[:,1]

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
	def length(self):
		return self._length/0.3048

	@property
	def csarea(self):
		return self._csarea/0.3048**2

	@property
	def numtot(self):
		return self._numtot

	@property
	def num(self):
		return self._num

	@property
	def xaxis(self):
		return self._xaxis/0.3048

	@property
	def size(self):
		return self._size/0.3048
	
	@property
	def area(self):
		return self._area/0.3048**2

	@property
	def perm(self):
		return self._perm/9.869233e-16

	@property
	def poro(self):
		return self._poro
	
if __name__ == "__main__":

	cells = Hexahedron(
		((-5,-5,5),(5,-5,5)),
		((5,-5,5),(15,-5,5)),
		((5,-5,-5),(15,-5,-5)),
		((-5,-5,-5),(5,-5,-5)),
		((-5,5,5),(5,5,5)),
		((5,5,5),(15,5,5)),
		((5,5,-5),(15,5,-5)),
		((-5,5,-5),(5,5,-5)),
	)

	# print(cells.unormal1)

	# print(cells.center)

	# print(cells.center1)
	# print(cells.center2)
	# print(cells.center3)
	# print(cells.center4)
	# print(cells.center5)
	# print(cells.center6)

	print(cells.volume)