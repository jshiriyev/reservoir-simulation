import lasio

from matplotlib import colors as mcolors
from matplotlib import gridspec
from matplotlib import pyplot
from matplotlib import transforms

from matplotlib.backends.backend_pdf import PdfPages

from matplotlib.patches import Polygon
from matplotlib.patches import Rectangle

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import LogFormatter
from matplotlib.ticker import LogFormatterExponent
from matplotlib.ticker import LogFormatterMathtext
from matplotlib.ticker import LogLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import NullLocator
from matplotlib.ticker import ScalarFormatter

import numpy

from ._items import Zone

from ._surface import Surface

class Formation():
    """It is a collection of surfaces"""

    def __init__(self):

        self.repo = {}

    def __setitem__(self,name,top):

        zone = Zone(name=name)

        zone.top = top

        self.repo[name] = zone

    def __getitem__(self,name):

        return self.repo[name]

    def tolist(self,prop):

        tops,props = [],[]

        for _,zone in self.repo.items():

            props.append(getattr(zone,prop))

            tops.append(zone.top)

        tops = numpy.array(tops)

        indices = numpy.argsort(tops)

        _props = []

        for index in indices:
            _props.append(props[index])

        return _props

    def view(self,axis=None):

        show = True if axis is None else False

        if axis is None:
            axis = pyplot.figure().add_subplot()

        tops = numpy.array(self.tolist("top"))

        depths = (tops[1:]+tops[:-1])/2

        axis.hlines(y=tops,xmin=0,xmax=1,color='k')

        axis.invert_yaxis()

        names = self.tolist("name")
        colors = self.tolist("color")
        hatches = self.tolist("hatch")

        for index,depth in enumerate(depths):

            axis.annotate(names[index],(0.5,depth),
                horizontalalignment='center',
                verticalalignment='center',
                backgroundcolor='white',)

            top,bottom = tops[index:index+2]

            axis.fill_between((0,1),top,y2=bottom,facecolor=colors[index],hatch=hatches[index])

        pyplot.setp(axis.get_xticklabels(),visible=False)
        pyplot.setp(axis.get_xticklines(),visible=False)

        if show:
            pyplot.show()

    @property
    def names(self):

        return self._names

    @property
    def tops(self):

        return self._tops

def TopsView(dictionary,axis=None,colors=None,hatches=None):

    show = True if axis is None else False

    if axis is None:
        axis = pyplot.figure().add_subplot()

    formations,tops = [],[]

    for key,value in dictionary.items():

        formations.append(key)
        
        tops.append(value)

    tops = numpy.array(tops)

    depths = (tops[1:]+tops[:-1])/2

    axis.hlines(y=tops,xmin=0,xmax=1,color='k')

    axis.invert_yaxis()

    for depth,formation in zip(depths,formations[:-1]):
        axis.annotate(formation,(0.5,depth),
                      horizontalalignment='center',
                      verticalalignment='center',
                      backgroundcolor='white',)

    for index,_ in enumerate(depths):

        x0 = (0,1)
        
        y1 = tops[index]

        y2 = tops[index+1]

        if colors is not None:
            color = colors[index]
        else:
            color = None

        if hatches is not None:
            hatch = hatches[index]
        else:
            hatch = None

        axis.fill_between(x0,y1,y2=y2,facecolor=color,hatch=hatch)

    pyplot.setp(axis.get_xticklabels(),visible=False)
    pyplot.setp(axis.get_xticklines(),visible=False)

    if show:
        pyplot.show()