class Conds():

	def __init__(self,*args):

		self.conds = ()

		for arg in args:
			self.add(arg)

	def add(self,cond):

		self.conds + (cond,)

	def __call__(self,tcurr):

		conds = [cond for cond in self.conds if self.islive(cond,tcurr)]

		return tuple(conds)

	@staticmethod
    def islive(cond,tcurr):

        if cond._start>tcurr:
            return False

        if cond._stop is None:
            return True

        if cond._stop<=tcurr:
            return False

        return True

