import contextlib

from io import StringIO

import re

import numpy

class GridRead():

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