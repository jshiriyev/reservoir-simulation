from ._items import Slot

class Stock():
    """It is a binary search tree of wells"""

    def __init__(self):

        self.wells = {}

    def add(self,slot:Slot):

        self.wells[slot.name] = slot

    def remove(self,name):

        self.wells.pop(name)

    def set_names(self,*args,fstring=None,sortFlag=False):

        warnNWIF = "No well name was added or could be found."

        [self.itemnames.append(str(arg)) for arg in args]

        if len(args)==0:

            twells = np.unique(self.Trajectory.get_wellnames())
            cwells = np.unique(self.Completion.get_wellnames())
            lwells = np.unique(self.Logging.get_wellnames())
            pwells = np.unique(self.Production.get_wellnames())

            self.itemnames = np.concatenate((twells,cwells,lwells,pwells)).tolist()

        self.number = len(self.itemnames)

        if self.number==0:
            warnings.warn(warnNWIF)
            return

        if wnamefstr is not None:
            self.wnamefstr = wnamefstr      # string format to save well names

        get_digits = lambda x: re.sub("[^0-9]","",x)

        get_digitnum = lambda x: len(get_digits(x))
        arr_digitnum = np.vectorize(get_digitnum)

        max_digitnum = arr_digitnum(self.itemnames).max()

        get_wellname = lambda x: self.wnamefstr.format(get_digits(x).zfill(max_digitnum))
        arr_wellname = np.vectorize(get_wellname)

        self.itemnames = arr_wellname(self.itemnames)

        # self.itemnames = np.unique(np.array(self.itemnames)).tolist()

        if sortFlag:
            self.itemnames.sort()

    def set_schedule(self,wellname):

        flagShowSteps = False if wellname is None else True

        warnNOPROD = "{} has completion but no production data."
        warnNOCOMP = "{} has production but no completion data."

        path1 = os.path.join(self.workdir,self.filename_op+"2")
        path2 = os.path.join(self.workdir,self.filename_comp+"1")
        path3 = os.path.join(self.workdir,self.filename_comp+"uni")

        self.op_get(filending="2")
        self.comp_get(filending="1")
        self.comp_get(filending="uni")

        prodwellnames = np.unique(self.op2.running[0])
        compwellnames = np.unique(self.comp1.running[0])

        for wname in np.setdiff1d(prodwellnames,compwellnames):
            warnings.warn(warnNOCOMP.format(wname))

        for wname in np.setdiff1d(compwellnames,prodwellnames):
            warnings.warn(warnNOPROD.format(wname))

        proddata = frame(headers=self.headers_op[:7])
        schedule = frame(headers=self.schedule_headers)

        for wname in self.itemnames:

            if wellname is not None:
                if wellname!=wname:
                    continue

            self.get_conflict(wname)

            shutdates = np.array(shutdates,dtype=object)
            shutwells = np.empty(shutdates.shape,dtype=object)

            shutwells[:] = wname

            shutdays = np.zeros(shutdates.shape,dtype=int)

            shutoptype = np.empty(shutdates.shape,dtype=object)

            shutoptype[:] = "shut"

            shutroil = np.zeros(shutdates.shape,dtype=int)
            shutrwater = np.zeros(shutdates.shape,dtype=int)
            shutrgas = np.zeros(shutdates.shape,dtype=int)

            rows = np.array([shutwells,shutdates,shutdays,shutoptype,shutroil,shutrwater,shutrgas]).T.tolist()

            proddata.set_rows(rows)

            if flagShowSteps:
                print("{} check is complete.".format(wname))

        proddata.sort(header_indices=[1],inplace=True)

        toil = np.cumsum(proddata.running[4])
        twater = np.cumsum(proddata.running[5])
        tgas = np.cumsum(proddata.running[6])

        proddata.set_column(toil,header_new="TOIL")
        proddata.set_column(twater,header_new="TWATER")
        proddata.set_column(tgas,header_new="TGAS")

        proddata.astype(header=self.headers_op[2],dtype=int)

        path = os.path.join(self.workdir,self.filename_op+"3")

        fstring = "{:6s}\t{:%Y-%m-%d}\t{:2d}\t{:10s}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\n"

        proddata.write(filepath=path,fstring=fstring)

    def get_conflicts(self,wellname):

        warnCROSS = "{} production has been defined before completion."

        warnWPGPF = "{:%Y-%m-%d}: {} first perf and last plug dates do not fit production days."
        warnWPERF = "{:%Y-%m-%d}: {} first perf date does not fit production days."
        warnWPLUG = "{:%Y-%m-%d}: {} last plug date does not fit production days."
        warnWEFAC = "{:%Y-%m-%d}: {} efficiency is more than unit [{:2d} out of {:2d} days]."

        self.op2.filter(0,keywords=[wname],inplace=False)

        self.comp1.filter(0,keywords=[wname],inplace=False)
        self.compuni.filter(0,keywords=[wname],inplace=False)

        try:
            datemin = self.op2.running[1].min()
        except ValueError:
            datemin = datetime(3000,1,1)

        date = datemin+relativedelta(months=1)

        days = calendar.monthrange(date.year,date.month)[1]

        date = datetime(date.year,date.month,days)

        if self.compuni.running[1].min()>=date:
            warnings.warn(warnCROSS.format(wname))

        schedule.set_rows([[self.comp1.running[1][0],"WELSPECS",self.schedule_welspecs.format(wname)]])

        for compdate,compevent,comptop,compbottom in zip(self.comp1.running[1],self.comp1.running[2],self.comp1.running[3],self.comp1.running[4]):

            if compevent == "PERF":
                schedule.set_rows([[compdate,"COMPDATMD",self.schedule_compdatop.format(wname,comptop,compbottom,"OPEN")]])
            elif compevent == "PLUG":
                schedule.set_rows([[compdate,"COMPDATMD",self.schedule_compdatsh.format(wname,comptop,"1*","SHUT")]])

        for compunidate in self.compuni.running[1]:

            schedule.set_rows([[compunidate,"COMPORD",self.schedule_compord.format(wname)]])

        flagNoPrevProd = True

        print("{} schedule is in progress ...".format(wname))

        opdata = zip(
            self.op2.running[1],
            self.op2.running[2],
            self.op2.running[3],
            self.op2.running[4],
            self.op2.running[5],
            self.op2.running[6],
            )

        shutdates = []

        for index,(date,days,optype,oil,water,gas) in enumerate(opdata):

            prodmonthSTARTday = date+relativedelta(days=1)

            prodmonthdaycount = calendar.monthrange(prodmonthSTARTday.year,prodmonthSTARTday.month)[1]

            prodmonthENDday = datetime(prodmonthSTARTday.year,prodmonthSTARTday.month,prodmonthdaycount)

            if np.sum(self.compuni.running[1]<prodmonthSTARTday)==0:
                compSTARTindex = 0
            else:
                compSTARTindex = np.sum(self.compuni.running[1]<prodmonthSTARTday)-1

            compENDindex = np.sum(self.compuni.running[1]<=prodmonthENDday)

            compupdatedates = self.compuni.running[1][compSTARTindex:compENDindex]
            compupdatecounts = self.compuni.running[2][compSTARTindex:compENDindex]

            perfdates = compupdatedates[compupdatecounts!=0]
            plugdates = compupdatedates[compupdatecounts==0]

            try:
                flagNoPostProd = True if self.op2.running[1][index+1]-relativedelta(months=1)>prodmonthENDday else False
            except IndexError:
                flagNoPostProd = True

            if np.sum(self.compuni.running[1]<prodmonthSTARTday)==0:
                flagCompShutSTART = True
            else:
                flagCompShutSTART = compupdatecounts[0]==0

            flagCompShutEND = compupdatecounts[-1]==0

            flagPlugPerf = any([compopencount==0 for compopencount in compupdatecounts[1:-1]])

            if flagCompShutSTART and flagCompShutEND:
                compday = plugdates[-1].day-perfdates[0].day
                prodeff = days/compday
                if optype == "production":
                    schedule.set_rows([[perfdates[0],"WCONHIST",self.schedule_prodhist.format(wname,oil,water,gas)]])
                elif optype == "injection":
                    schedule.set_rows([[perfdates[0],"WCONINJH",self.schedule_injhist.format(wname,water)]])
                schedule.set_rows([[perfdates[0],"WEFAC",self.schedule_wefac.format(wname,prodeff)]])
                proddata.set_rows([[wname,perfdates[0],days,optype,oil,water,gas]])
                schedule.set_rows([[plugdates[-1],"WELOPEN",self.schedule_welopen.format(wname)]])
                shutdates.append(plugdates[-1])
                flagNoPrevProd = True
                if flagShowSteps:
                    print("{:%d %b %Y} Peforated and Plugged: OPEN ({:%d %b %Y}) and SHUT ({:%d %b %Y}) WEFAC ({:.3f})".format(prodmonthENDday,perfdates[0],plugdates[-1],prodeff))

            elif flagCompShutSTART:
                compday = prodmonthENDday.day-perfdates[0].day
                prodeff = days/compday
                if optype == "production":
                    schedule.set_rows([[perfdates[0],"WCONHIST",self.schedule_prodhist.format(wname,oil,water,gas)]])
                elif optype == "injection":
                    schedule.set_rows([[perfdates[0],"WCONINJH",self.schedule_injhist.format(wname,water)]])
                schedule.set_rows([[perfdates[0],"WEFAC",self.schedule_wefac.format(wname,prodeff)]])
                proddata.set_rows([[wname,perfdates[0],days,optype,oil,water,gas]])
                if flagNoPostProd:
                    schedule.set_rows([[prodmonthENDday,"WELOPEN",self.schedule_welopen.format(wname)]])
                    shutdates.append(prodmonthENDday)
                    flagNoPrevProd = True
                    if flagShowSteps:
                        print("{:%d %b %Y} Peforated and Open: OPEN ({:%d %b %Y}) and SHUT ({:%d %b %Y}) WEFAC ({:.3f})".format(prodmonthENDday,perfdates[0],prodmonthENDday,prodeff))
                else:                  
                    flagNoPrevProd = False
                    if flagShowSteps:
                        print("{:%d %b %Y} Peforated and Open: OPEN ({:%d %b %Y}) and CONT WEFAC ({:.3f})".format(prodmonthENDday,perfdates[0],prodeff))

            elif flagCompShutEND:
                for plugdate in plugdates:
                    if plugdate.day>=days: break
                if not plugdate.day>=days:
                    warnings.warn(warnWPLUG.format(prodmonthENDday,wname))
                compday = plugdate.day
                prodeff = days/compday
                if optype == "production":
                    schedule.set_rows([[date,"WCONHIST",self.schedule_prodhist.format(wname,oil,water,gas)]])
                elif optype == "injection":
                    schedule.set_rows([[date,"WCONINJH",self.schedule_injhist.format(wname,water)]])
                schedule.set_rows([[date,"WEFAC",self.schedule_wefac.format(wname,prodeff)]])
                proddata.set_rows([[wname,date,days,optype,oil,water,gas]])
                schedule.set_rows([[plugdate,"WELOPEN",self.schedule_welopen.format(wname)]])
                shutdates.append(plugdate)
                flagNoPrevProd = True
                if flagShowSteps:
                    print("{:%d %b %Y} Open and Plugged: CONT and SHUT ({:%d %b %Y}) WEFAC ({:.3f})".format(prodmonthENDday,plugdate,prodeff))

            elif flagPlugPerf:
                if flagNoPrevProd and flagNoPostProd:
                    # shift the start day to the first perf day
                    # shut the well at the last plug day
                    if not plugdates[-1].day-perfdates[1].day>=days:
                        warnings.warn(warnWPGPF.format(prodmonthENDday,wname))
                    compday = plugdates[-1].day-perfdates[1].day
                    prodeff = days/compday
                    if optype == "production":
                        schedule.set_rows([[perfdates[1],"WCONHIST",self.schedule_prodhist.format(wname,oil,water,gas)]])
                    elif optype == "injection":
                        schedule.set_rows([[perfdates[1],"WCONINJH",self.schedule_injhist.format(wname,water)]])
                    schedule.set_rows([[perfdates[1],"WEFAC",self.schedule_wefac.format(wname,prodeff)]])
                    proddata.set_rows([[wname,perfdates[1],days,optype,oil,water,gas]])
                    schedule.set_rows([[plugdates[-1],"WELOPEN",self.schedule_welopen.format(wname)]])
                    shutdates.append(plugdates[-1])
                    flagNoPrevProd = True
                    if flagShowSteps:
                        print("{:%d %b %Y} Plugged and Perforated: OPEN ({:%d %b %Y}) and SHUT ({:%d %b %Y}) WEFAC ({:.3f})".format(prodmonthENDday,perfdates[1],plugdates[-1],prodeff))
                elif flagNoPrevProd and not flagNoPostProd:
                    # shift the start day to the proper perf day
                    for perfdate in np.flip(perfdates[1:]):
                        if prodmonthENDday.day-perfdate.day>=days: break
                    if not prodmonthENDday.day-perfdate.day>=days:
                        warnings.warn(warnWPERF.format(prodmonthENDday,wname))
                    compday = prodmonthENDday.day-perfdate.day
                    prodeff = days/compday
                    if optype == "production":
                        schedule.set_rows([[perfdate,"WCONHIST",self.schedule_prodhist.format(wname,oil,water,gas)]])
                    elif optype == "injection":
                        schedule.set_rows([[perfdate,"WCONINJH",self.schedule_injhist.format(wname,water)]])
                    schedule.set_rows([[perfdate,"WEFAC",self.schedule_wefac.format(wname,prodeff)]])
                    proddata.set_rows([[wname,perfdate,days,optype,oil,water,gas]])
                    flagNoPrevProd = False
                    if flagShowSteps:
                        print("{:%d %b %Y} Plugged and Perforated: OPEN ({:%d %b %Y}) and CONT WEFAC ({:.3f})".format(prodmonthENDday,perfdate,prodeff))
                elif not flagNoPrevProd and flagNoPostProd:
                    # try shut the well at the proper plug day if not successful shut it at the end of month
                    for plugdate in plugdates:
                        if plugdate.day>=days: break
                    if not plugdate.day>=days:
                        plugdate = prodmonthENDday
                    compday = plugdate.day
                    prodeff = days/compday
                    if optype == "production":
                        schedule.set_rows([[date,"WCONHIST",self.schedule_prodhist.format(wname,oil,water,gas)]])
                    elif optype == "injection":
                        schedule.set_rows([[date,"WCONINJH",self.schedule_injhist.format(wname,water)]])
                    schedule.set_rows([[date,"WEFAC",self.schedule_wefac.format(wname,prodeff)]])
                    proddata.set_rows([[wname,date,days,optype,oil,water,gas]])
                    schedule.set_rows([[plugdate,"WELOPEN",self.schedule_welopen.format(wname)]])
                    shutdates.append(plugdate)
                    flagNoPrevProd = True
                    if flagShowSteps:
                        print("{:%d %b %Y} Plugged and Perforated: CONT and SHUT ({:%d %b %Y}) WEFAC ({:.3f})".format(prodmonthENDday,plugdate,prodeff))
                elif not flagNoPrevProd and not flagNoPostProd:
                    # try shut the well if not successful do nothing
                    for plugdate in plugdates:
                        if plugdate.day>=days: break
                    if not plugdate.day>=days:
                        compday = prodmonthdaycount
                        prodeff = days/compday
                        flagNoPrevProd = False
                        if flagShowSteps:
                            print("{:%d %b %Y} Plugged and Perforated: CONT and CONT WEFAC ({:.3f})".format(prodmonthENDday,prodeff))
                    else:
                        compday = plugdate.day
                        prodeff = days/compday
                        schedule.set_rows([[plugdate,"WELOPEN",self.schedule_welopen.format(wname)]])
                        shutdates.append(plugdate)
                        flagNoPrevProd = True
                        if flagShowSteps:
                            print("{:%d %b %Y} Plugged and Perforated: CONT and SHUT ({:%d %b %Y}) WEFAC ({:.3f})".format(prodmonthENDday,plugdate,prodeff))
                    if optype == "production":
                        schedule.set_rows([[date,"WCONHIST",self.schedule_prodhist.format(wname,oil,water,gas)]])
                    elif optype == "injection":
                        schedule.set_rows([[date,"WCONINJH",self.schedule_injhist.format(wname,water)]])
                    schedule.set_rows([[date,"WEFAC",self.schedule_wefac.format(wname,prodeff)]])
                    proddata.set_rows([[wname,date,days,optype,oil,water,gas]])

            else:
                compday = prodmonthdaycount
                prodeff = days/compday
                if optype == "production":
                    schedule.set_rows([[date,"WCONHIST",self.schedule_prodhist.format(wname,oil,water,gas)]])
                elif optype == "injection":
                    schedule.set_rows([[date,"WCONINJH",self.schedule_injhist.format(wname,water)]])
                schedule.set_rows([[date,"WEFAC",self.schedule_wefac.format(wname,prodeff)]])
                proddata.set_rows([[wname,date,days,optype,oil,water,gas]])
                if flagNoPostProd:
                    schedule.set_rows([[prodmonthENDday,"WELOPEN",self.schedule_welopen.format(wname)]])
                    shutdates.append(prodmonthENDday)
                    flagNoPrevProd = True
                    if flagShowSteps:
                        print("{:%d %b %Y} No completion events: CONT and SHUT ({:%d %b %Y}) WEFAC ({:.3f})".format(prodmonthENDday,prodmonthENDday,prodeff))
                else:
                    flagNoPrevProd = False
                    if flagShowSteps:
                        print("{:%d %b %Y} No completion events: CONT and CONT WEFAC ({:.3f})".format(prodmonthENDday,prodeff))

            if prodeff>1:
                warnings.warn(warnWEFAC.format(prodmonthENDday,wname,days,compday))

    def set_itemnames(self,namelist,fstring=None,zfill=3):
        
        fstring = "{}" if fstring is None else fstring

        getwname = lambda x: fstring.format(re.sub(r"[^\d]","",str(x)).zfill(zfill))

        getwname = np.vectorize(getwname)

        self.itemnames = getwname(namelist)

    def set_distance(self,depth=None):

        coords = np.zeros((len(self.files),3))

        for index,data in enumerate(self.files):

            if depth is None:
                depthIndex = 0

            coords[index,:] = data[depthIndex,:3]

        dx = coords[:,0]-coords[:,0].reshape((-1,1))
        dy = coords[:,1]-coords[:,1].reshape((-1,1))
        dz = coords[:,2]-coords[:,2].reshape((-1,1))

        self.distance = np.sqrt(dx**2+dy**2+dz**2)

    def get_kneighbors(self,k=1):

        min_indices = np.zeros((self.distance.shape[0],k),dtype=int)

        for index_self,row in enumerate(self.distance):

            indices = np.argpartition(row,range(k+1))[:k+1]

            min_indices[index_self,:] = np.delete(indices,indices==index_self)

        return min_indices

    def set_tracks(self,tracks):

        self.tracks = np.array(tracks)