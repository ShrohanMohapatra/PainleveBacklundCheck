from sympy import *
alpha = -1
x, t, sigma, b = symbols('x t sigma b')
f = Function('f')(x, t)
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
	# print(zeta)
	# print(generator_list)
	polynomial = Poly(expr, gens=generator_list)
	degree_list = polynomial.degree_list()
	# print(degree_list)
	min_deg, max_deg = degree_list[0], degree_list[1]
	min_deg = -min_deg
	degree_coeff_tuple = [
		[k, simplify(expr.coeff(zeta, k))]
		for k in range(min_deg,max_deg+1)]
	return degree_coeff_tuple
j = Symbol('j')
phi = Function('phi')
M = int(-alpha + 3) # Maximum number of terms
U = [Function('U'+str(k)) for k in range(M+1)]
u = 0
for j in range(M+1):
	u = u + phi(x, t)**(j+alpha)*U[j](x, t)
# print(u, function_PDE)
function_PDE = parse_expr('diff(f(x,t),t)+f(x,t)*diff(f(x,t),x)-sigma*diff(f(x,t),x,x)')
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
backlund_transform = u.subs([
		(U[m](x,t), U_final[m]) if m<(-alpha) else \
		((U[m](x,t), U[m](x,t)) if m==(-alpha) else (U[m](x,t), 0))\
		for m in range(M+1)
		])
# Compatibility conditions
compatibility_Conditions = [
	simplify(totalTable[k][1].subs([
	(U[m](x,t), U_final[m]) if m<(-alpha) else \
	((U[m](x,t), U[m](x,t)) if m==(-alpha) else (U[m](x,t), 0))\
	for m in range(M+1)
	]).doit()) \
	for k in range(len(totalTable))
]
refined_compatibility_Conditions = []
for k in range(len(totalTable)):
	print('--------------------------------------------')
	print('--------------------------------------------')
	print(k, '>>', compatibility_Conditions[k])
	print('--------------------------------------------')
	if compatibility_Conditions[k] != 0:
		Factors = factor_list(compatibility_Conditions[k])
		print('======>>', Factors)
		print('Extracted factor --->')
		print(Factors[1][-1][0])
		refined_compatibility_Conditions.append(
			compatibility_Conditions[k]
			)
	print('--------------------------------------------')
