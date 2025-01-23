import importlib

import numpy

def zfactor(critical_params:tuple,pressures:numpy.ndarray,temperature:float,derivative:bool=False,method:str="direct_method"):
	"""
	Dynamically calculates z-factor based on the specified method.

	Parameters:
	    critical_params (tuple): tuple of (pcrit in psi, tcrit in Rankine)

	    pressures (numpy.ndarray): Array of pressure values (psi) to calculate z-factor for.
	    temperature (float): The temperature value (Rankine) at which to calculate z-factor.
	    
	    method_name (str): The name of the method to use:
	    	(e.g., 'direct_method', 'dranchuk_abu_kassem', 'dranchuk_purvis_robinson', 'hall_yarborough').
	    
	Returns:
	    numpy.ndarray: Z-factor, and if derivative is True, Z-prime values calculated for each pressure.
	    
	"""

	method_name_parts = method.split('_')

	method_class_name = method_name_parts[0].capitalize()+''.join(word.capitalize() for word in method_name_parts[1:])
	print(f"{method_class_name=}")

	method_file_name = f"_{method}"
	print(f"{method_file_name=}")

	# Import the correct method class dynamically
	try:

		# full_module_name = f"compressibility.{method_class_name}"

		# module = importlib.import_module(full_module_name)

		module = __import__(f'respy.fluid_props.phaseg.compressibility',fromlist=[method_class_name])

		MethodClass = getattr(module, method_class_name)

	except (ImportError, AttributeError):
		raise ValueError(f"Method '{method}' not found or invalid.")

	# Create an instance of the class and calculate z-factor
	method_instance = MethodClass(critical_params,temperature)

	return method_instance(pressures,derivative)

if __name__ == "__main__":

	#phaseg.zfactor((500,200),(550,600,700),250)

	import matplotlib.pyplot as plt

	import numpy as np

	dm = DirectMethod((500,100),152)

	Pr,Tr = 2.99,1.52

	# print(Orkahub_Energy(Pr,Tr,derivative=True))
	# print(Hall_Yarborough(Pr,Tr,derivative=True))
	# print(Dranchuk_Abu_Kassem(Pr,Tr,derivative=True))
	# print(Dranchuk_Purvis_Robinson(Pr,Tr,derivative=True))

	pressure = np.linspace(1500,2000,200)

	z1,d1 = dm(pressure,derivative=True)
	# z2,d2 = HallYarborough(Pr,Tr,derivative=True)
	# z3,d3 = DranchukAbuKassem(Pr,Tr,derivative=True)
	# z4,d4 = DranchukPurvisRobinson(Pr,Tr,derivative=True)

	plt.plot(pressure,z1,label = "Direct Method")
	# plt.plot(Pr,z2,label = "Hall_Yarborough")
	# plt.plot(Pr,z3,label = "Dranchuk_Abu_Kassem")
	# plt.plot(Pr,z4,label = "Dranchuk_Purvis_Robinson")

	# cr1 = 1/(1+(0.27*Pr)/(z1**2*Tr)*d1)/Pr
	# cr2 = 1/(1+(0.27*Pr)/(z2**2*Tr)*d2)/Pr
	# cr3 = 1/(1+(0.27*Pr)/(z3**2*Tr)*d3)/Pr
	# cr4 = 1/(1+(0.27*Pr)/(z4**2*Tr)*d4)/Pr

	plt.legend()

	plt.show()
