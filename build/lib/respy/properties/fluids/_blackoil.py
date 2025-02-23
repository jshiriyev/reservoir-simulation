class BlackOil():

	def __init__(self,*args):

		for index,arg in enumerate(args):

			setattr(self,f"phase{index}",arg)

	def set_rperm(self,*args):

		for index in range(3):

			setattr(self,f"rperm{index}",args[index])

	def set_capil(self,*args):

		pass


