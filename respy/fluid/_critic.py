class critic:

	@staticmethod
	def natural(spgr):
		"""It calculates pseudo-critical temperature and pressure for
		natural-gas systems based on specific gravity:

		spgr    : specific gravity of the gas at
                  standard conditions.

        Output:

		ppc     : Critical Pressure in psia
		tpc     : Critical Temperature in °R
		"""
		ppc = 677+15*spgr-37.5*spgr**2
		tpc = 168+325*spgr-12.5*spgr**2

		return ppc,tpc

	@staticmethod
	def condensate(spgr):
		"""It calculates pseudo-critical temperature and pressure for
		gas-condensate systems based on specific gravity:

		spgr    : specific gravity of the gas at
                  standard conditions.

        Output:

		ppc     : Critical Pressure in psia
		tpc     : Critical Temperature in °R
		"""
		ppc = 706-51.7*spgr-11.1*spgr**2
		tpc = 187+330.0*spgr-71.5*spgr**2

		return ppc,tpc

	@staticmethod
	def highweight(spgr):
		"""It calculates pseudo-critical properties for high molecular
		weight reservoir gases based on specific gravity:

		spgr    : specific gravity of the gas at
                  standard conditions.

        Output:

		ppc     : Critical Pressure in psia
		tpc     : Critical Temperature in °R
		"""
		ppc = 756.8-131*spgr-3.6*spgr**2
		tpc = 169.2+349.5*spgr-74*spgr**2

		return ppc,tpc

if __name__ == "__main__":

	print(critic.natural(0.8))
	print(critic.condensate(0.8))
	print(critic.highweight(0.8))