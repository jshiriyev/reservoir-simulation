from respy.fluid._phase import Phase

class OnePhase():

	def __init__(self,*args,**kwargs):

		if len(args)==0:
			self.phase = Phase(**kwargs)
		elif len(args)==1:
			self.phase = args[0]
		else:
			raise Warning("Cannot assign more than one phase!")

	def __call__(self,press,**kwargs):

		if callable(self.phase):
			return self.phase(press)

		self.phase._press = press

		return self.phase

	@property
    def isstatic(self):
        return not callable(self.phase)

    @property
    def isdynamic(self):
    	return callable(self.phase)
    