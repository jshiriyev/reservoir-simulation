from respy.utils._prop import Prop

from respy.solver._resinit import ResInit

class ResRock():

    _gravity = 9.807  # Gravitational acceleration in SI units

    def __init__(self,grid,rcomp=1e-5,**kwargs):
        """
        grid    : RecCuboid instance, rectangular cuboid grids

        """

        self.grid = grid

        self.__rcomp = Prop(rcomp,0.000145038)

        self.resinit = ResInit(**kwargs)

    def init(self,pcow,pcog,pcgw):

    	Pw0 = self.resinit.waterpressure(self.grid._depth)

    	Po0 = self.resinit.oilpressure(self.grid._depth)
    	Pg0 = self.resinit.gaspressure(self.grid._depth)

    	Sw0,So0,Sg0 = self.resinit.saturations(pcow,pcog,pcgw)

    	cw0 = Sw0*self.fluids['water']._comp

    	co0 = So0*self.fluids['oil']._comp

    	cg0 = Sg0*self.fluids['gas']._comp

    	ct0 = self.get_rcomp(None)+cw0+co0+cg0

    def __call__(self,press=None,temp=None):

        mats = {}

        for key,fluid in self.fluids.items():

            phase = fluid(press,temp)

            T = self.get_Tmatrix(phase)
            S = self.get_Smatrix()
            G = self.get_Gvector(phase)
            J = self.get_Jmatrix(phase)
            Q = self.get_Qvector(phase)

            mats[key] = Matrix(T,V,G,J,Q)

        return mats

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

    def set_rcomp(self,P):
        return self.__rcomp(P)

    @property
    def shape(self):
        """Shape of the matrices of transmissibility calculations"""
        return (self.grid.numtot,self.grid.numtot)

    @property
    def size(self):
        """Shape of the vectors of transmissibility calculations"""
        return self.grid.numtot

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
        return ResRockSub(self,1,edge=True)

    @property
    def xpos(self):
        """Properties of grids that has x-positive neighbors"""
        return ResRockSub(self,2)

    @property
    def xneg(self):
        """Properties of grids that has x-negative neighbors"""
        return ResRockSub(self,1)

    @property
    def xmax(self):
        """Properties of grids on x-maximum boundary"""
        return ResRockSub(self,2,edge=True)

    @property
    def ymin(self):
        """Properties of grids on y-minimum boundary"""
        if self.flodim>1:
            return ResRockSub(self,3,edge=True)

    @property
    def ypos(self):
        """Properties of grids that has y-positive neighbors"""
        if self.flodim>1:
            return ResRockSub(self,4)

    @property
    def yneg(self):
        """Properties of grids that has y-negative neighbors"""
        if self.flodim>1:
            return ResRockSub(self,3)

    @property
    def ymax(self):
        """Properties of grids on y-maximum boundary"""
        if self.flodim>1:
            return ResRockSub(self,4,edge=True)

    @property
    def zmin(self):
        """Properties of grids on z-minimum boundary"""
        if self.flodim>2:
            return ResRockSub(self,5,edge=True)

    @property
    def zpos(self):
        """Properties of grids that has z-positive neighbors"""
        if self.flodim>2:
            return ResRockSub(self,6)

    @property
    def zneg(self):
        """Properties of grids that has z-negative neighbors"""
        if self.flodim>2:
            return ResRockSub(self,5)

    @property
    def zmax(self):
        """Properties of grids on z-maximum boundary"""
        if self.flodim>2:
            return ResRockSub(self,6,edge=True)

class ResRockSub(numpy.ndarray):

    def __new__(cls,grid,path,edge=False):
        """
        grid    : RecCuboid instance
        path    : direction (1: west, 2: east, 3: south, 4: north, 5: down, 6: up)
        edge    : boundary (True) or inner (False)
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
        """("dims","area","volume","depth","poro","perm")"""

        unit = key[1:] if key.startswith('_') else key  

        if unit in ("volume","depth","poro"):
            return getattr(self.grid,key)[self,0]

        if unit[0] not in ("x","y","z"):
            return

        if unit[1:] in ("dims","area","perm"):
            return getattr(self.grid,key)[self,0]