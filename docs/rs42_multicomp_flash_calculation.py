def rachford_rice(v):

	term1 = 0.65/(1+0.55*v)
	term2 = 0.30/(1-0.36*v)
	term3 = 0.05/(1-0.72*v)

	return term1+term2+term3-1

def rachford_rice_derivative(v):

	term1 = 0.3575/(1+0.55*v)**2
	term2 = 0.1080/(1-0.36*v)**2
	term3 = 0.0360/(1-0.72*v)**2

	return -term1+term2+term3

def newton(func,der1,initial,tol=1e-5,args=None):

	if args is None:
		args = []

	d0 = func(initial,*args)
	d1 = der1(initial,*args)

	calls = 1

	while abs(d0/d1)>tol:

		print(f"{calls:2} {initial:8.6f},{d0:8.6f},{d1:8.6f},{d0/d1:8.6f}")

		initial -= d0/d1

		d0 = func(initial,*args)
		d1 = der1(initial,*args)

		calls += 1

	info = f"Converged after {calls} calls."

	minima = d0

	return initial

print(newton(rachford_rice,rachford_rice_derivative,0.6))