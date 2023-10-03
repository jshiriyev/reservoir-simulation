import contextlib

from io import StringIO

import re

import numpy

from .directory._browser import Browser

class GridFile(Browser):

    def __init__(self,**kwargs):

        super().__init__(**kwargs)

        self.sections = {}

    def __getattr__(self,key):

        return getattr(self.sections,key)

    def __getitem__(self,key):

        return self.sections[key]

def loadgrid(*args,**kwargs):
    """
    Returns an instance of GridFile. If a filepath is specified, the instance
    represents the file.
    
    Arguments:
        filepath {str} -- path to the given grdecl file

    Keyword Arguments:
        homedir {str} -- path to the home (output) directory
        filedir {str} -- path to the file (input) directory
    
    Returns:
        GridFile -- an instance of GridFile filled with grdecl file text.

    """

    if len(args)==1:
        filepath = args[0]
    elif len(args)>1:
        raise "The function does not take more than one positional argument."

    # It creates an empty GridFile instance.
    nullfile = GridFile(filepath=filepath,**kwargs)

    # It reads grdecl file and returns GridFile instance.
    fullfile = GridWorm(nullfile).gridfile

    return fullfile

class GridWorm():

    def __init__(self,gridfile):

        self.gridfile = gridfile

        self.text(self.gridfile.filepath)

        self.process()

    def process(self):

        array = numpy.loadtxt(self.gridfile.sections['COORD'])

        self.gridfile.sections['COORD'] = array

        array = numpy.genfromtxt(StringIO(" ".join(self.gridfile.sections['ZCORN'])))

        self.gridfile.sections['ZCORN'] = array

        array = numpy.genfromtxt(StringIO(" ".join(self.gridfile.sections['ACTNUM'])),dtype=int)

        self.gridfile.sections['ACTNUM'] = array

    def headers(self,master):

        while True:

            try:
                line = next(master)
            except StopIteration:
                break

            line = line.strip()

            if len(line)<1:
                continue
            
            if line.startswith("--"):
                continue

            if line.startswith("/"):
                continue

            if line=="NOECHO":
                continue

            if line=="ECHO":
                break

            match = re.match(r'(?:[A-Z]+|[a-z]+)$',line)

            if match is not None:
                rows = self._header(master)

                if match.string=="INCLUDE":
                    [self.text(row) for row in rows]
                else:
                    self.gridfile.sections[match.string] = rows

    def _header(self,master):

        rows = []

        while True:

            line = next(master).strip()

            row = line.split("/")

            line = row[0].strip()

            if len(line)>0:
                rows.append(line)

            if len(row)>1:
                break

        return rows

    def text(self,filepath):

        with self.txtopen(filepath) as gridmaster:
            self.headers(gridmaster)

    @staticmethod
    @contextlib.contextmanager
    def txtopen(filepath):

        txtmaster = open(filepath,"r")

        try:
            yield txtmaster
        finally:
            txtmaster.close()

def starsplit(string_list,default=1.0):
    """It returns star splitted list repeating post-star pre-star times."""

    float_list = []

    for string_value in string_list:

        if "*" in string_value:

            if string_value.endswith("*"):
                mult = string_value.rstrip("*")
                mult,val = int(mult),default
            else:
                mult,val = string_value.split("*",maxsplit=1)
                mult,val = int(mult),float(val)

            for i in range(mult):
                float_list.append(val)

        else:
            float_list.append(float(string_value))

    return float_list