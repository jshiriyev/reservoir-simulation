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

	method_parts = method.split('_')

	method_class = method_parts[0].capitalize()+''.join(word.capitalize() for word in method_parts[1:])

	# Import the correct method class dynamically
	try:
		module = __import__(f'respy.fluid_props.phaseg.compressibility',fromlist=[method_class])
		MethodClass = getattr(module, method_class)
	except (ImportError, AttributeError):
		raise ValueError(f"Method '{method}' not found or invalid.")

	# Create an instance of the class and calculate z-factor
	method_instance = MethodClass(critical_params,temperature)

	return method_instance(pressures,derivative)