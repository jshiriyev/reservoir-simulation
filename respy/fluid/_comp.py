class composition():
    """
    Gas can be defined based on specific gravity or molecular composition.
    
    Second option, molecular composition can be defined in dictionary (**kwargs) and it may contain:

        component       : abbreviation of component ('CH4','C2H6','H2S',etc.)
        mfracts         : mole fraction of each component
        mweight         : molecular weight of each component
        Tcritical       : critical temperature for each component
        Pcritical       : critical pressure for each component
        .               :
        .               :
        .               :
    """

    def __init__(self,**kwargs):
        """keys are parameters 'str', and values are fields 'float' or 'None'."""

        if len(kwargs)==0:
            raise ValueError("At least one field is required.")

        params = []
        fields = []
        fsizes = []

        for param,field in kwargs.items():

            params.append(param)

            if isinstance(field,list) or isinstance(field,tuple):
                field = list(field)
            else:
                field = [field]

            fields.append(field)
            fsizes.append(len(field))

        if len(set(fsizes))!=1:
            raise ValueError("The lengths of field are not equal!")

        super().__setattr__("params",params)
        super().__setattr__("fields",fields)

        super().__setattr__("elements",library['elements'])

    def setsum(self,param,value):

        index = self.params.index(param)

        field = self.fields[index]

        factor = value/sum(field)

        field_new = [number*factor for number in field]

        self.fields[index] = field_new

    def fillna(self):

        for rindex,row in enumerate(self):
            for cindex,param in enumerate(self.params):
                if cindex==0:
                    continue
                if row[cindex] is None:
                    self.fields[cindex][rindex] = self.elements[row[0]][param]

    def extend(self,row):

        if len(self.params)!=len(row):
            raise ValueError("The lengths of 'fields' and 'row' are not equal!")

        if isinstance(row,list) or isinstance(row,tuple):
            toextend = mixture(**dict(zip(self.params,row)))
        elif isinstance(row,dict):
            toextend = mixture(**row)
        elif isinstance(row,mixture):
            toextend = row

        for param,field in toextend.items():
            self.fields[self.params.index(param)].extend(field)

    def __setattr__(self,param,field):

        if isinstance(field,list):
            pass
        if isinstance(field,tuple):
            field = list(field)
        else:
            field = [field]

        number_of_rows = len(self)

        if len(field)!=number_of_rows:
            raise AttributeError(f"field len does not fit the len of rows in Mixture object, {len(field)}!={number_of_rows}.")
        
        if param in self.params:
            self.fields[self.params.index(param)] = field
        else:
            self.params.append(param)
            self.fields.append(field)

    def __getattr__(self,param):

        field = self.fields[self.params.index(param)]

        if len(field)==1:
            field, = field

        return field

    def __setitem__(self,key,row):

        if len(self.params)!=len(row)+1:
            raise ValueError("The lengths of 'fields' and 'row' are not equal!")

        if isinstance(row,list) or isinstance(row,tuple):
            toset = mixture(**dict(zip(self.params,row)))
        elif isinstance(row,dict):
            toset = mixture(**row)
        elif isinstance(row,mixture):
            toset = row

        for index,row in enumerate(self):
            if row[0].lower()==key.lower():
                break

        for param,field in zip(self.params,self.fields):
            field[index] = getattr(toset,param)[0]

    def __getitem__(self,key):

        if not isinstance(key,str):
            raise TypeError("key must be string!")

        for row in self:
            if row[0].lower()==key.lower():
                break
        else:
            return
        
        return mixture(**dict(zip(self.params,row)))

    def __repr__(self,comment=None):

        return self.__str__(comment)

    def __str__(self,comment=None):

        if len(self)==1:
            return repr(tuple(self.fields))

        if comment is None:
            comment = ""

        fstring = comment
        
        underline = []

        for param,field in zip(self.params,self.fields):

            field_ = list(field)
            field_.append(param)

            count_ = max([len(str(value)) for value in field_])

            fstring += f"{{:<{count_}}}   "
            
            underline.append("-"*count_)

        fstring += "\n"

        text = fstring.format(*[parm.capitalize() for parm in self.params])
        
        text += fstring.format(*underline)
        
        for row in self:
            text += fstring.format(*row)

        return text

    def __iter__(self):

        return iter([row for row in zip(*self.fields)])

    def __len__(self):

        return len(self.fields[0])

    def items(self):

        return iter([(p,f) for p,f in zip(self.params,self.fields)])

    @staticmethod
    def get_mwa(mfracts,mweight):
        """Calculates apparent molecular weight."""
        return sum([fract*weight for fract,weight in zip(mfracts,mweight)])

    @staticmethod
    def get_spgr_atsc(mwa):
        """The calculation assumes that the behavior of both the gas mixture and
        air is described by the ideal gas equation at standard conditions."""
        return mwa/28.964