from sympy import *
init_printing()
x, t, sigma, b = symbols('x t sigma b')
nu = symbols('nu')
e1, e2, e3, d, r = symbols('e1 e2 e3 d r')
f = Function('f')(x, t)
# Painleve property and Backlund transform
def totalTableOfCoeffsForPoly(expr, zeta):
	generator_list = [1/zeta, zeta]
	for m in range(10):
		for n in range(10):
			if m == 0 and n == 0: pass
			elif m == 0 and n == 1:
				generator_list.append(Derivative(zeta, t))
			elif m == 1 and n == 0:
				generator_list.append(Derivative(zeta, x))
			elif m == 1 and n == 1:
				generator_list.append(Derivative(zeta, t, x))
				# generator_list.append(Derivative(zeta, t, x))
			elif m == 1 and n > 1:
				generator_list.append(Derivative(zeta, (t, n), x))
				# generator_list.append(Derivative(zeta, (t, n), x))
			elif n == 1 and m > 1:
				generator_list.append(Derivative(zeta, t, (x, m)))
				# generator_list.append(Derivative(zeta, (x, m), t))
			elif m == 0 and n > 1:
				generator_list.append(Derivative(zeta, (t, n)))
			elif n == 0 and m > 1:
				generator_list.append(Derivative(zeta, (x, m)))
			elif m > 1 and n > 1:
				generator_list.append(
					Derivative(zeta, (t, n), (x, m))
					)
				# generator_list.append(
				# 	Derivative(zeta, (t, n), (x, m))
				# 	)
	print(zeta)
	print(generator_list)
	polynomial = Poly(expr, gens=generator_list)
	degree_list = polynomial.degree_list()
	# print(degree_list)
	min_deg, max_deg = degree_list[0], degree_list[1]
	min_deg = -min_deg
	degree_coeff_tuple = [
		[k, simplify(expr.coeff(zeta, k))]
		for k in range(min_deg,max_deg+1)]
	return degree_coeff_tuple
def PainlevePDE_withBacklundTransform(
	function_PDE, # Input PDE
	x, # The input variable used for the spatial variable
	t, # The input variable used for the temporal variable
	alpha # The power balancing exponent
	):
	j = Symbol('j')
	phi = Function('phi')
	M = int(-alpha + 3) # Maximum number of terms
	U = [Function('U'+str(k)) for k in range(M+1)]
	u = 0
	for j in range(M+1):
		u = u + phi(x, t)**(j+alpha)*U[j](x, t)
	# print(u, function_PDE)
	orig_eqn = expand(function_PDE.subs(f, u).doit())
	# print(orig_eqn)
	totalTable = totalTableOfCoeffsForPoly(orig_eqn, phi(x, t))
	# print(totalTable)
	U1 = [U[k](x, t) for k in range(M+1)]
	U_final = [0 for k in range(M+1)]
	for k in range(len(totalTable)):
		# print(U_final)
		index, equation = totalTable[k]
		if equation == 0 or k >= M+1: continue
		if k == 0:
			try:
				soln = pdsolve(equation, U[k](x, t))
			except:
				soln = solve(equation, U[k](x, t))
		else:
			try:
				soln = pdsolve(
					equation.subs(
						[(U[m](x,t), U_final[m]) for m in range(k)]),
					U[k](x, t)
				)
			except:
				soln = solve(
					equation.subs(
						[(U[m](x,t), U_final[m]) for m in range(k)]),
					U[k](x, t)
				)
		if len(soln) == 0: break
		else:
			for m in range(len(soln)):
				if soln[m] == 0: continue
				else:
					U_final[k] = simplify(soln[m])
					break
	# Final expression for the Backlund transform
	backlund_transform = u.subs([
		(U[m](x,t), U_final[m]) if m<(-alpha) else \
		((U[m](x,t), U[m](x,t)) if m==(-alpha) else (U[m](x,t), 0))\
		for m in range(M+1)
		])
	# Compatibility conditions
	compatibility_Conditions = [
		[
			simplify(totalTable[k][1].subs([
			(U[m](x,t), U_final[m]) if m < cutoff_index else \
			((U[m](x,t), U[m](x,t)) if cutoff_index <= m <= (-alpha) else (U[m](x,t), 0))\
			for m in range(M+1)
			]).doit()) \
			for k in range(len(totalTable))
		] for cutoff_index in range(1, (-alpha)+1)
	]
	refined_compatibility_Conditions = [[] for cutoff_index in range(1, (-alpha)+1)]
	for cutoff_index in range(1, (-alpha)+1):
		for k in range(len(totalTable)):
			if compatibility_Conditions[cutoff_index-1][k] != 0:
				refined_compatibility_Conditions[cutoff_index-1].append(
					compatibility_Conditions[cutoff_index-1][k]
					)
	return [backlund_transform, refined_compatibility_Conditions]
def ultimatePainleve(
	function_PDE, # Input PDE
	x, # The input variable used for the spatial variable
	t # The input variable used for the temporal variable
	):
	range_of_exponents = [-1, -2, -3, -4]
	flag = False
	for alpha in range_of_exponents:
		M = int(-alpha + 3) # Maximum number of terms
		U = [Function('U'+str(k)) for k in range(M+1)]
		backlund_transform, compatibility_Conditions = PainlevePDE_withBacklundTransform(function_PDE, x, t, alpha)
		# print(alpha, compatibility_Conditions, backlund_transform)
		# print(alpha, len(compatibility_Conditions))
		if len(compatibility_Conditions[-1]) < 2: continue
		else:
			flag = expand(compatibility_Conditions[-1][-1].subs(U[-alpha](x, t),f).doit()) == expand(function_PDE)
			# print(compatibility_Conditions[-1][-1].subs(U[-alpha](x, t),f).doit(), expand(function_PDE))
			# print(alpha, flag)
			if not flag: continue
			else: break
	# print(alpha, compatibility_Conditions, backlund_transform)
	if alpha == -4 and len(compatibility_Conditions[-1]) == 1 and backlund_transform == U[4](x,t):
		return (None, backlund_transform, compatibility_Conditions)
	for k in range(len(compatibility_Conditions[-1])):
		if str(compatibility_Conditions[-1][k]).find(str(U[-alpha](x, t)))==-1:
			return (None, backlund_transform, [])				
	if not flag:
		alpha, compatibility_Conditions = None, []
	return (alpha, backlund_transform, compatibility_Conditions)
