import numpy

# from .directory._browser import Browser

@dataclass(frozen=True)
class Well:                 # It is a well dictionary used in the simulator
    name    : str           # name of the well
    block   : tuple         # block indices containing the well 
    sort    : str           # vertical or horizontal

    radius  : float         # well radius, ft
    skin    : float = 0     # skin factor of the well, dimensionless

    conds   : list  = field(default_factory=list)

    def __post_init__(self):
        object.__setattr__(self,'_radius',self.radius*0.3048)

    def set_survey(self,survey):
        self.survey = survey

    def add_conds(self,*args,**kwargs):
        """Adds a bottom hole condition to the conds property"""
        self.conds.append(WellConds(*args,**kwargs))

@dataclass(frozen=True)
class WellConds:            # bottom hole conditions
    time    : float         # the time for implementing the condition, days
    status  : str           # status of the well, open or shut
    bhp     : float = None  # constant bottom hole pressure if that is the control
    orate   : float = None  # constant oil rate if that is the control
    wrate   : float = None  # constant water rate if that is the control
    grate   : float = None  # constant gas rate if that is the control

class WellStock():

    def __init__(self):

        self.well_grids      = np.array([],dtype=int)
        self.well_indices    = np.array([],dtype=int)
        self.well_bhpflags   = np.array([],dtype=bool)
        self.well_limits     = np.array([],dtype=float)

    def __iter__(self):

        wells = zip(self.Wells.Trajectory.tracks,
            self.Wells.consbhp,
            self.Wells.limits,
            )

        for index,(track,flag,limit) in enumerate(wells):
            pass

    def set(self):

        for index,(track,flag,limit) in enumerate(wells):

            ttrack = np.transpose(track[:,:,np.newaxis],(2,1,0))

            vector = self.PorRock.grid_centers[:,:,np.newaxis]-ttrack

            distance = np.sqrt(np.sum(vector**2,axis=1))

            well_grid = np.unique(np.argmin(distance,axis=0))

            well_index = np.full(well_grid.size,index,dtype=int)
    
            well_bhpflag = np.full(well_grid.size,flag,dtype=bool)

            well_limit = np.full(well_grid.size,limit,dtype=float)

            self.well_grids = np.append(self.well_grids,well_grid)

            self.well_indices = np.append(self.well_indices,well_index)

            self.well_bhpflags = np.append(self.well_bhpflags,well_bhpflag)

            self.well_limits = np.append(self.well_limits,well_limit)

        dx = self.PorRock.grid_sizes[self.well_grids,0]
        dy = self.PorRock.grid_sizes[self.well_grids,1]
        dz = self.PorRock.grid_sizes[self.well_grids,2]

        req = 0.14*np.sqrt(dx**2+dy**2)

        rw = self.Wells.radii[self.well_indices]

        skin = self.Wells.skinfactors[self.well_indices]

        kx = self.PorRock.permeability[:,0][self.well_grids]
        ky = self.PorRock.permeability[:,1][self.well_grids]

        dz = self.PorRock.grid_sizes[:,2][self.well_grids]

        self.JR = (2*np.pi*dz*np.sqrt(kx*ky))/(np.log(req/rw)+skin)

class Schedule(Browser): # previously it was inheritting from Browser

    def __init__(self,*args,**kwargs):

        super().__init__(*args,**kwargs)

    def write(self):

        path = os.path.join(self.workdir,self.schedule_filename)

        with open(path,"w",encoding='utf-8') as wfile:

            welspec = schedule.running[1]=="WELSPECS"
            compdat = schedule.running[1]=="COMPDATMD"
            compord = schedule.running[1]=="COMPORD"
            prodhst = schedule.running[1]=="WCONHIST"
            injdhst = schedule.running[1]=="WCONINJH"
            wefffac = schedule.running[1]=="WEFAC"
            welopen = schedule.running[1]=="WELOPEN"

            for date in numpy.unique(schedule.running[0]):

                currentdate = schedule.running[0]==date

                currentcont = schedule.running[1][currentdate]

                wfile.write("\n\n")
                wfile.write("DATES\n")
                wfile.write(self.schedule_dates.format(date.strftime("%d %b %Y").upper()))
                wfile.write("\n")
                wfile.write("/\n\n")

                if any(currentcont=="WELSPECS"):
                    indices = numpy.logical_and(currentdate,welspec)
                    wfile.write("WELSPECS\n")
                    for detail in schedule.running[2][indices]:
                        wfile.write(detail)
                        wfile.write("\n")
                    wfile.write("/\n\n")

                if any(currentcont=="COMPDATMD"):
                    indices = numpy.logical_and(currentdate,compdat)
                    wfile.write("COMPDATMD\n")
                    for detail in schedule.running[2][indices]:
                        wfile.write(detail)
                        wfile.write("\n")
                    wfile.write("/\n\n")

                if any(currentcont=="COMPORD"):
                    indices = numpy.logical_and(currentdate,compord)
                    wfile.write("COMPORD\n")
                    for detail in schedule.running[2][indices]:
                        wfile.write(detail)
                        wfile.write("\n")
                    wfile.write("/\n\n")

                if any(currentcont=="WCONHIST"):
                    indices = numpy.logical_and(currentdate,prodhst)
                    wfile.write("WCONHIST\n")
                    for detail in schedule.running[2][indices]:
                        wfile.write(detail)
                        wfile.write("\n")
                    wfile.write("/\n\n")

                if any(currentcont=="WCONINJH"):
                    indices = numpy.logical_and(currentdate,injdhst)
                    wfile.write("WCONINJH\n")
                    for detail in schedule.running[2][indices]:
                        wfile.write(detail)
                        wfile.write("\n")
                    wfile.write("/\n\n")

                if any(currentcont=="WEFAC"):
                    indices = numpy.logical_and(currentdate,wefffac)
                    wfile.write("WEFAC\n")
                    for detail in schedule.running[2][indices]:
                        wfile.write(detail)
                        wfile.write("\n")
                    wfile.write("/\n\n")

                if any(currentcont=="WELOPEN"):
                    indices = numpy.logical_and(currentdate,welopen)
                    wfile.write("WELOPEN\n")
                    for detail in schedule.running[2][indices]:
                        wfile.write(detail)
                        wfile.write("\n")
                    wfile.write("/\n\n")

def loadsched():

    pass

class SchedWorm():

    def __init__(self,filepath):

        pass

    def read(self,skiprows=0,headerline=None,comment="--",endline="/",endfile="END"):

        # While looping inside the file it does not read lines:
        # - starting with comment phrase, e.g., comment = "--"
        # - after the end of line phrase, e.g., endline = "/"
        # - after the end of file keyword e.g., endfile = "END"

        if headerline is None:
            headerline = skiprows-1
        elif headerline<skiprows:
            headerline = headerline
        else:
            headerline = skiprows-1

        _running = []

        with open(self.filepath,"r") as text:

            for line in text:

                line = line.split('\n')[0].strip()

                line = line.strip(endline)

                line = line.strip()
                line = line.strip("\t")
                line = line.strip()

                if line=="":
                    continue

                if comment is not None:
                    if line[:len(comment)] == comment:
                        continue

                if endfile is not None:
                    if line[:len(endfile)] == endfile:
                        break

                _running.append([line])

        self.title = []

        for _ in range(skiprows):
            self.title.append(_running.pop(0))

        num_cols = len(_running[0])

        if skiprows==0:
            self.set_headers(num_cols=num_cols,init=False)
        elif skiprows!=0:
            self.set_headers(headers=self.title[headerline],init=False)

        nparray = numpy.array(_running).T

    def catch(self,header_index=None,header=None,regex=None,regex_builtin="INC_HEADERS",title="SUB-HEADERS"):

        nparray = numpy.array(self._running[header_index])

        if regex is None and regex_builtin=="INC_HEADERS":
            regex = r'^[A-Z]+$'                         #for strings with only capital letters no digits
        elif regex is None and regex_builtin=="INC_DATES":
            regex = r'^\d{1,2} [A-Za-z]{3} \d{2}\d{2}?$'   #for strings with [1 or 2 digits][space][3 capital letters][space][2 or 4 digits], e.g. DATES

        vmatch = numpy.vectorize(lambda x: bool(re.compile(regex).match(x)))

        match_index = vmatch(nparray)

        firstocc = numpy.argmax(match_index)

        lower = numpy.where(match_index)[0]
        upper = numpy.append(lower[1:],nparray.size)

        repeat_count = upper-lower-1

        match_content = nparray[match_index]

        nparray[firstocc:][~match_index[firstocc:]] = numpy.repeat(match_content,repeat_count)

        self._headers.insert(header_index,title)
        self._running.insert(header_index,numpy.asarray(nparray))

        for index,datacolumn in enumerate(self._running):
            self._running[index] = numpy.array(self._running[index][firstocc:][~match_index[firstocc:]])

        self.headers = self._headers
        self.running = [numpy.asarray(datacolumn) for datacolumn in self._running]

    def program(self):

        # KEYWORDS: DATES,COMPDATMD,COMPORD,WCONHIST,WCONINJH,WEFAC,WELOPEN 

        dates      = " {} / "#.format(date)
        welspecs   = " '{}'\t1*\t2* / "
        compdatop  = " '{}'\t1*\t{}\t{}\tMD\t{}\t2*\t0.14 / "#.format(wellname,top,bottom,optype)
        compdatsh  = " '{}'\t1*\t{}\t{}\tMD\t{} / "#.format(wellname,top,bottom,optype)
        compord    = " '{}'\tINPUT\t/ "#.format(wellname)
        prodhist   = " '{}'\tOPEN\tORAT\t{}\t{}\t{} / "#.format(wellname,oilrate,waterrate,gasrate)
        injhist    = " '{}'\tWATER\tOPEN\t{}\t7*\tRATE / "#.format(wellname,waterrate)
        wefac      = " '{}'\t{} / "#.format(wellname,efficiency)
        welopen    = " '{}'\tSHUT\t3* / "#.format(wellname)