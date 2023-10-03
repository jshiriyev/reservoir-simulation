from matplotlib import pyplot

import numpy

from ._items import Slot
from ._items import Zone

class Well(Slot):

	def __init__(self,*args,**kwargs):

		super().__init__(*args,**kwargs)

		self.surveys = []

		self.zones = []

		self.lasfiles = []

	def add_survey(self,survey):

		self.surveys.append(survey)

	def add_zones(self,**zones):

		for name,top in zones.items():

			zone = Zone(name=name)

			zone.top = float(top)

			self.zones.append(zone)

	def add_las(self,lasfile):

		self.lasfiles.append(lasfile)

	def viewlasmap(self,axis=None,index=0):

		if axis is None:

			self.figure = pyplot.figure()

			axis = self.figure.add_subplot()

		self.axis = axis

		zones = [zone.name for zone in self.zones]

		mnemonics,harvest = self.get_lasharvest(index)

		image = self.axis.imshow(harvest,cmap="Greens")#"YlGn")

		# Show all ticks and label them with the respective list entries
		self.axis.set_xticks(numpy.arange(len(zones)),labels=zones)
		self.axis.set_yticks(numpy.arange(len(mnemonics)),labels=mnemonics)

		self.axis.xaxis.tick_top()

		pyplot.setp(self.axis.get_xticklabels(),rotation=60,ha="left",va="center",
         	rotation_mode="anchor")

		self.figure.set_tight_layout(True)
		
		pyplot.show()

	def savelasmap(self,filepath,title=None,sizemult=0.4):

		self.figure = pyplot.figure()

		self.axis = self.figure.add_subplot()

		zones = [zone.name for zone in self.zones]

		mnemonics,harvest = self.get_lasharvest()

		image = self.axis.imshow(harvest,cmap="Greens",vmin=0,vmax=200)

		# Show all ticks and label them with the respective list entries
		self.axis.set_xticks(numpy.arange(len(zones)),labels=zones)
		self.axis.set_yticks(numpy.arange(len(mnemonics)),labels=mnemonics)

		self.axis.xaxis.tick_top()

		pyplot.setp(self.axis.get_xticklabels(),rotation=60,ha="left",va="center",
			rotation_mode="anchor")

		if title is None:
			title = self.name

		self.axis.set_title(title)

		self.figure.set_figwidth(sizemult*len(zones))

		self.figure.set_figheight(sizemult*len(mnemonics))

		self.figure.set_tight_layout(True)

		self.figure.savefig(filepath)

		pyplot.close()

	def get_lasharvest(self,index=None):

		harvest = []

		if index is None:
			lasfiles = self.lasfiles
		else:
			lasfiles = [self.lasfiles[index]]

		mnemonics = self.get_mnemonics(lasfiles)

		harvest = numpy.zeros((len(mnemonics),len(self.zones)))

		for lasfile in lasfiles:

			depth = lasfile[0]

			for curve in lasfile.curves[1:]:

				mindex = mnemonics.index(curve.mnemonic)

				row = []

				for zindex,zone in enumerate(self.zones):

					limit = depth>zone.top

					if zindex!=len(self.zones)-1:

						limit = numpy.logical_and(limit,depth<self.zones[zindex+1].top)

					data = curve.data[limit]

					if data.size == 0:
						number = 0
					else:
						number = numpy.sum(~numpy.isnan(data))/data.size*100

					row.append(number)

				harvest[mindex] = numpy.maximum(harvest[mindex],row)
					
		return mnemonics,harvest

	def get_mnemonics(self,lasfiles):

		mnemonics = []

		for lasfile in lasfiles:

			for curve in lasfile.curves[1:]:

				mnemonics.append(curve.mnemonic)

		return numpy.unique(mnemonics).tolist()