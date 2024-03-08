class CornPoint():

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

		return Hexahedron(p1,p2,p3,p4,p5,p6,p7,p8)
	
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