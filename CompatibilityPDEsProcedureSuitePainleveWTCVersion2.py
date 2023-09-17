from sympy import *
from collections import Counter
import unittest

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
		simplify(totalTable[k][1].subs([
			(U[m](x,t), U_final[m]) if m<(-alpha) else \
			((U[m](x,t), U[m](x,t)) if m==(-alpha) else (U[m](x,t), 0))\
			for m in range(M+1)
		]).doit()) \
		for k in range(len(totalTable))
	]
	refined_compatibility_Conditions = []
	for k in range(len(totalTable)):
		if compatibility_Conditions[k] != 0:
			refined_compatibility_Conditions.append(
				compatibility_Conditions[k]
			)
	# refined_compatibility_Conditions.append(
		# Eq(U[-alpha](x, t), phi(x, t)))
	return [backlund_transform, refined_compatibility_Conditions]

def ultimatePainleve(
	function_PDE, # Input PDE
	x, # The input variable used for the spatial variable
	t, # The input variable used for the temporal variable
	):
	range_of_exponents = [-1, -2, -3, -4]
	flag = False
	for alpha in range_of_exponents:
		M = int(-alpha + 3) # Maximum number of terms
		U = [Function('U'+str(k)) for k in range(M+1)]
		backlund_transform, compatibility_Conditions = PainlevePDE_withBacklundTransform(function_PDE, x, t, alpha)
		print(alpha, compatibility_Conditions, backlund_transform)
		if len(compatibility_Conditions) < 2: continue
		else:
			flag = expand(compatibility_Conditions[-1].subs(U[-alpha](x, t),f).doit()) == expand(function_PDE)
			print(compatibility_Conditions[-1].subs(U[-alpha](x, t),f).doit(), expand(function_PDE))
			print(alpha, flag)
			if not flag: continue
			else: break
	print(alpha, compatibility_Conditions, backlund_transform)
	if alpha == -4 and len(compatibility_Conditions) == 1 and backlund_transform == U[4](x,t):
		return (None, backlund_transform, compatibility_Conditions)
	for k in range(len(compatibility_Conditions)):
		if str(compatibility_Conditions[k]).find(str(U[-alpha](x, t)))==-1:
			return (None, backlund_transform, [])
	if not flag:
		alpha, compatibility_Conditions = None, []
	return (alpha, backlund_transform, compatibility_Conditions)

def max_order_term(original_equation, func):
	terms_from_the_original_equation = original_equation.as_ordered_terms()
	order_of_original_equation = ode_order(original_equation, func)
	for individual_term in terms_from_the_original_equation:
		if ode_order(individual_term, func) == order_of_original_equation:
			factorList = factor_list(individual_term)
			for factor in factorList[1]:
				if ode_order(factor[0], func) == order_of_original_equation:
					actual_term = factor[0]**factor[1]
					full_term = individual_term
					break
			break
	actual_max_order_term = actual_term
	coefficient = simplify(full_term/actual_max_order_term)
	return [coefficient, actual_term, factor[0], ode_order(factor[0], func)]

def simplificationPDE(inputPDE1, inputPDE2, func):
	max_order_term1 = max_order_term(inputPDE1, func)
	max_order_term2 = max_order_term(inputPDE2, func)
	A1 = max_order_term1[0] # max_order_term1[0]/(max_order_term1[1]**max_order_term1[2])
	A2 = max_order_term2[0] # max_order_term2[0]/(max_order_term2[1]**max_order_term2[2])
	B1 = inputPDE1 - A1*max_order_term1[2]
	B2 = inputPDE2 - A2*max_order_term2[2]

	inputPDE1 = expand(-B1/A1) # (max_order_term1[1]**max_order_term1[2])
	inputPDE1, _ = fraction(together(inputPDE1.doit()))
	inputPDE1 = expand(inputPDE1)
	inputPDE2 = expand(-B2/A2) # (max_order_term2[1]**max_order_term2[2])
	inputPDE2, _ = fraction(together(inputPDE2.doit()))
	inputPDE2 = expand(inputPDE2)


	print('----------------------------------------')
	print('Here we compute the commutation of the following two PDEs:')
	print('----------------------------------------')
	print('Eqn (1) >>')
	print('----------------------------------------')
	print(inputPDE1, '=', 0)
	print('----------------------------------------')
	print('----------------------------------------')
	print('Eqn (2) >>')
	print('----------------------------------------')
	print(inputPDE2, '=', 0)
	print('----------------------------------------')
	print('The maximum order term in equation 1 is gamma_1 = ', max_order_term1[2])
	print('----------------------------------------')
	print('----------------------------------------')
	print('The maximum order term in equation 2 is gamma_2 = ', max_order_term2[2])
	print('----------------------------------------')
	print('In this scheme we represent eqn1 = A1*gamma_1 + B1 and eqn2 = A2*gamma_2 + B2 where')
	print('----------------------------------------')
	print('A1 =', A1)
	print('----------------------------------------')
	print('A2 =', A2)
	print('----------------------------------------')
	print('B1 =', B1)
	print('----------------------------------------')
	print('B2 =', B2)
	print('----------------------------------------')

	allAvailableVariables = func.free_symbols

	varsFromPDE1 = list(max_order_term1[2].variables) if max_order_term1[3] != 0 else []
	varsFromPDE2 = list(max_order_term2[2].variables) if max_order_term2[3] != 0 else []

	print('-------------------------------')
	print('The variables of differentiation in gamma_1')
	print('-------------------------------')
	print(varsFromPDE1)
	print('-------------------------------')
	print('-------------------------------')
	print('The variables of differentiation in gamma_2')
	print('-------------------------------')
	print(varsFromPDE2)
	print('-------------------------------')

	dictFromPDE1 = Counter(varsFromPDE1)
	dictFromPDE2 = Counter(varsFromPDE2)

	for variable in dictFromPDE1:
		if variable not in dictFromPDE2:
			dictFromPDE2[variable] = 0

	for variable in dictFromPDE2:
		if variable not in dictFromPDE1:
			dictFromPDE1[variable] = 0

	dictLCM_MixedPartialDeriv = {}
	dictFromPDE3 = {}
	dictFromPDE4 = {}
	for variable in dictFromPDE1:
		dictLCM_MixedPartialDeriv[variable] = max(dictFromPDE1[variable], dictFromPDE2[variable])
		dictFromPDE3[variable] = dictLCM_MixedPartialDeriv[variable] - dictFromPDE1[variable]
		dictFromPDE4[variable] = dictLCM_MixedPartialDeriv[variable] - dictFromPDE2[variable]

	varsFromPDE3 = []
	varsFromPDE4 = []

	for variable in dictFromPDE3:
		for k in range(dictFromPDE3[variable]):
			varsFromPDE3.append(variable)

	for variable in dictFromPDE4:
		for k in range(dictFromPDE4[variable]):
			varsFromPDE4.append(variable)

	print('-------------------------------')
	print('The variables of differentiation for the expression -B1/A1 from equation 1 are:')
	print('-------------------------------')
	print(varsFromPDE3)
	print('-------------------------------')
	print('The variables of differentiation for the expression -B2/A2 from equation 1 are:')
	print('-------------------------------')
	print('varsFromPDE4')
	print('-------------------------------')
	print(varsFromPDE4)
	print('-------------------------------')

	inputPDE3 = inputPDE1
	for variable in varsFromPDE3:
		inputPDE3 = diff(inputPDE3, variable)
		inputPDE3 = expand(inputPDE3.doit())

	inputPDE3 = expand(inputPDE3)
	inputPDE3, _ = fraction(together(inputPDE3.doit()))
	inputPDE3 = expand(inputPDE3)

	inputPDE4 = inputPDE2
	for variable in varsFromPDE4:
		inputPDE4 = diff(inputPDE4, variable)
		inputPDE4 = expand(inputPDE4.doit())

	inputPDE4 = expand(inputPDE4)
	inputPDE4, _ = fraction(together(inputPDE4.doit()))
	inputPDE4 = expand(inputPDE4)

	print('----------------------------')
	print('Expression after differentiation of equation 1')
	print('----------------------------')
	print(inputPDE3)
	print('----------------------------')
	print('----------------------------')
	print('----------------------------')
	print('Expression after differentiation of equation 2')
	print('----------------------------')
	print(inputPDE4)
	print('----------------------------')
	print('----------------------------')

	expression = inputPDE3 - inputPDE4
	expression, _ = fraction(together(expression.doit()))
	expression = expand(expression)

	print('----------------------------')
	print('Expression after commutation')
	print('----------------------------')
	print(expression)
	print('----------------------------')

	return expression

def CompatibilitySystemPDEs(systemsOfPDEs, func, depth_max = 15, depth = 0):

	print('------------------------------------------------------------------------------------------------')
	print('Welcome to recursion number = ',depth+1)
	print('------------------------------------------------------------------------------------------------')

	for k in range(len(systemsOfPDEs)):
		print('--------------------------------------------------------------------------------------------')
		print('>> PDE Number '+str(k+1)+' -->')
		print('--------------------------------------------------------------------------------------------')
		print(systemsOfPDEs[k])
		print('--------------------------------------------------------------------------------------------')

	orders_list = [max_order_term(equation, func)[3] for equation in systemsOfPDEs]
	if len(systemsOfPDEs) == 2:
		maximum_order = max(orders_list)
		maximum_order_index = orders_list.index(maximum_order)
		equation_2_new = expand(simplificationPDE(systemsOfPDEs[0], systemsOfPDEs[1], func))
		if equation_2_new == 0:
			return [depth, True]
		else:
			if depth >= depth_max: return [depth, False]
			try:
				return CompatibilitySystemPDEs([systemsOfPDEs[maximum_order_index], equation_2_new], func, depth_max, depth + 1)
			except Exception:
				return [depth, False]
	else:
		systemsOfPDEs_F = [None for k in range(len(systemsOfPDEs))]
		for k in range(len(systemsOfPDEs)-1):
			systemsOfPDEs_F[k] = expand(simplificationPDE(systemsOfPDEs[k], systemsOfPDEs[k+1], func))
		systemsOfPDEs_F[-1] = expand(simplificationPDE(systemsOfPDEs[-1], systemsOfPDEs[0], func))
		flag = False
		for k in range(len(systemsOfPDEs_F)):
			print('---------------------------------------------------')
			print('Number of terms in ',k, 'th PDE >>', len(Add.make_args(systemsOfPDEs_F[k])))
			print('---------------------------------------------------')
			print(systemsOfPDEs_F[k])
			print('---------------------------------------------------')
			if len(Add.make_args(systemsOfPDEs_F[k])) > 250: return [depth, False]
			flag = flag or systemsOfPDEs_F[k].equals(0)
		if flag: return [depth, flag]
		else:
			if depth >= depth_max: return [depth, False]
			try:
				return CompatibilitySystemPDEs(systemsOfPDEs_F, func, depth_max, depth + 1)
			except Exception:
				return [depth, False]

class TestCompatibilityPDE(unittest.TestCase):

	# def test_TestCaseBurgers(self):
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	print('Here I am testing for the compatibility of the system of PDEs')
	# 	print('that come from the WTC Painleve analysis of the Burgers equation ....')
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	phi, f, U1 = symbols('phi f U1', cls = Function)
	# 	BurgersPDE = parse_expr('diff(f(x,t),t)+f(x,t)*diff(f(x,t),x)-sigma*diff(f(x,t),x,x)')
	# 	alpha, transform, compatibility_equations = ultimatePainleve(
	# 			BurgersPDE, # Input PDE
	# 			x, # The input variable used for the spatial variable
	# 			t # The input variable used for the temporal variable
	# 			# alpha # The power balancing exponent
	# 		)
	# 	BurgersPDE_eqn = BurgersPDE.subs(f(x, t), transform)
	# 	BurgersPDE_eqn, _ = fraction(together(BurgersPDE_eqn.doit()))
	# 	BurgersPDE_eqn = expand(BurgersPDE_eqn)
	# 	compatibility_equations = [compatibility_equations[k] for k in range(1, len(compatibility_equations))]
	# 	compatibility_equations.insert(0, BurgersPDE_eqn)
	# 	flag = False
	# 	for k in range(len(compatibility_equations)):
	# 		flag = flag or compatibility_equations[k] == 0
	# 	if flag:
	# 		print('The PDE system from WTC Painleve analysis of the Burgers equation renders compatible')
	# 		self.assertTrue(True)
	# 	else:
	# 		for k in range(len(compatibility_equations)):
	# 			Factors = factor_list(compatibility_equations[k])
	# 			compatibility_equations[k] = Factors[1][-1][0]
	# 			print('------------------------------')
	# 			print('Compatibility equation number = ', k + 1)
	# 			print('------------------------------')
	# 			print(compatibility_equations[k])
	# 			print('------------------------------')
	# 		compatibility_checker = CompatibilitySystemPDEs(compatibility_equations, phi(x, t), 6)
	# 		depth = compatibility_checker[0]
	# 		status = compatibility_checker[1]
	# 		if status:
	# 			print('The PDE system from WTC Painleve analysis of the Burgers equation renders compatible')
	# 		else:
	# 			print('The PDE system from WTC Painleve analysis of the Burgers equation renders incompatible')
	# 		print('------------------------------------------------------------------------------------------------------------------------------------')
	# 		print('Number of cycles of reduction in Burgers equation', depth+1)
	# 		self.assertTrue(status)

	# def test_TestCaseKdV(self):
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	print('Here I am testing for the compatibility of the system of PDEs')
	# 	print('that come from the WTC Painleve analysis of the KdV equation ....')
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	phi, f, U2 = symbols('phi f U2', cls = Function)
	# 	KdVPDE = parse_expr('diff(f(x,t),t)+f(x,t)*diff(f(x,t),x)+sigma*diff(f(x,t),x,x,x)')
	# 	alpha, transform, compatibility_equations = ultimatePainleve(
	# 			KdVPDE, # Input PDE
	# 			x, # The input variable used for the spatial variable
	# 			t # The input variable used for the temporal variable
	# 			# alpha # The power balancing exponent
	# 		)
	# 	KdVPDE_eqn = KdVPDE.subs(f(x, t), transform)
	# 	KdVPDE_eqn, _ = fraction(together(KdVPDE_eqn.doit()))
	# 	KdVPDE_eqn = expand(KdVPDE_eqn)
	# 	compatibility_equations = [compatibility_equations[k] for k in range(1, len(compatibility_equations))]
	# 	compatibility_equations.insert(0, KdVPDE_eqn)
	# 	for k in range(len(compatibility_equations)):
	# 		Factors = factor_list(compatibility_equations[k])
	# 		compatibility_equations[k] = Factors[1][-1][0]
	# 		print('------------------------------')
	# 		print('Compatibility equation number = ', k + 1)
	# 		print('------------------------------')
	# 		print(compatibility_equations[k])
	# 		print('------------------------------')
	# 	compatibility_checker = CompatibilitySystemPDEs(compatibility_equations, phi(x, t), 6)
	# 	depth = compatibility_checker[0]
	# 	status = compatibility_checker[1]
	# 	if status:
	# 		print('The PDE system from WTC Painleve analysis of the KdV equation renders compatible')
	# 	else:
	# 		print('The PDE system from WTC Painleve analysis of the KdV equation renders incompatible')
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	print('Number of cycles of reduction in KdV equation', depth+1)
	# 	self.assertTrue(status)

	# def test_TestCaseMKdV(self):
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	print('Here I am testing for the compatibility of the system of PDEs')
	# 	print('that come from the WTC Painleve analysis of the modified KdV equation ....')
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	phi, f, U1 = symbols('phi f U1', cls = Function)
	# 	MKdVPDE = parse_expr('diff(f(x,t),t)+f(x,t)**2*diff(f(x,t),x)+sigma*diff(f(x,t),x,x,x)')
	# 	alpha, transform, compatibility_equations = ultimatePainleve(
	# 			MKdVPDE, # Input PDE
	# 			x, # The input variable used for the spatial variable
	# 			t # The input variable used for the temporal variable
	# 			# alpha # The power balancing exponent
	# 		)
	# 	MKdVPDE_eqn = MKdVPDE.subs(f(x, t), transform)
	# 	MKdVPDE_eqn, _ = fraction(together(MKdVPDE_eqn.doit()))
	# 	MKdVPDE_eqn = expand(MKdVPDE_eqn)
	# 	compatibility_equations = [compatibility_equations[k] for k in range(1, len(compatibility_equations))]
	# 	compatibility_equations.insert(0, MKdVPDE_eqn)
	# 	for k in range(len(compatibility_equations)):
	# 		Factors = factor_list(compatibility_equations[k])
	# 		compatibility_equations[k] = Factors[1][-1][0]
	# 		print('------------------------------')
	# 		print('Compatibility equation number = ', k + 1)
	# 		print('------------------------------')
	# 		print(compatibility_equations[k])
	# 		print('------------------------------')
	# 	compatibility_checker = CompatibilitySystemPDEs(compatibility_equations, phi(x, t), 6)
	# 	depth = compatibility_checker[0]
	# 	status = compatibility_checker[1]
	# 	if status:
	# 		print('The PDE system from WTC Painleve analysis of the modified KdV equation renders compatible')
	# 	else:
	# 		print('The PDE system from WTC Painleve analysis of the modified KdV equation renders incompatible')
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	print('Number of cycles of reduction in modified KdV equation', depth+1)
	# 	self.assertTrue(status)

	# def test_TestCaseBBM(self):
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	print('Here I am testing for the compatibility of the system of PDEs')
	# 	print('that come from the WTC Painleve analysis of the BBM (Benjamin-Bona-Mahony) equation ....')
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	phi, f, U2 = symbols('phi f U2', cls = Function)
	# 	BBM_PDE = parse_expr('diff(f(x, t), t)+diff(f(x, t), x)+f(x, t)*diff(f(x, t), x)-diff(f(x, t), x, x, t)')
	# 	alpha, transform, compatibility_equations = ultimatePainleve(
	# 			BBM_PDE, # Input PDE
	# 			x, # The input variable used for the spatial variable
	# 			t # The input variable used for the temporal variable
	# 			# alpha # The power balancing exponent
	# 		)
	# 	BBMPDE_eqn = BBM_PDE.subs(f(x, t), transform)
	# 	BBMPDE_eqn, _ = fraction(together(BBMPDE_eqn.doit()))
	# 	BBMPDE_eqn = expand(BBMPDE_eqn)
	# 	compatibility_equations = [compatibility_equations[k] for k in range(1, len(compatibility_equations))]
	# 	compatibility_equations.insert(0, BBMPDE_eqn)
	# 	for k in range(len(compatibility_equations)):
	# 		Factors = factor_list(compatibility_equations[k])
	# 		compatibility_equations[k] = Factors[1][-1][0]
	# 		print('------------------------------')
	# 		print('Compatibility equation number = ', k + 1)
	# 		print('------------------------------')
	# 		print(compatibility_equations[k])
	# 		print('------------------------------')
	# 	compatibility_checker = CompatibilitySystemPDEs(compatibility_equations, phi(x, t), 6)
	# 	depth = compatibility_checker[0]
	# 	status = compatibility_checker[1]
	# 	if status:
	# 		print('The PDE system from WTC Painleve analysis of the BBM equation renders compatible')
	# 	else:
	# 		print('The PDE system from WTC Painleve analysis of the BBM equation renders incompatible')
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	self.assertFalse(status)

	# def test_TestCaseKuramotoSivashinsky(self):
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	print('Here I am testing for the compatibility of the system of PDEs')
	# 	print('that come from the WTC Painleve analysis of the KS (Kuramoto-Sivashinsky) equation ....')
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	phi, f, U3 = symbols('phi f U3', cls = Function)
	# 	KS_PDE = parse_expr('diff(f(x, t), t)+f(x, t)*diff(f(x, t), x)+diff(f(x, t), x, x)+diff(f(x, t), x, x, x, x)')
	# 	alpha, transform, compatibility_equations = ultimatePainleve(
	# 			KS_PDE, # Input PDE
	# 			x, # The input variable used for the spatial variable
	# 			t # The input variable used for the temporal variable
	# 			# alpha # The power balancing exponent
	# 		)
	# 	KSPDE_eqn = KS_PDE.subs(f(x, t), transform)
	# 	KSPDE_eqn, _ = fraction(together(KSPDE_eqn.doit()))
	# 	KSPDE_eqn = expand(KSPDE_eqn)
	# 	compatibility_equations = [compatibility_equations[k] for k in range(1, len(compatibility_equations))]
	# 	compatibility_equations.insert(0, KSPDE_eqn)
	# 	for k in range(len(compatibility_equations)):
	# 		Factors = factor_list(compatibility_equations[k])
	# 		compatibility_equations[k] = Factors[1][-1][0]
	# 		print('------------------------------')
	# 		print('Compatibility equation number = ', k + 1)
	# 		print('------------------------------')
	# 		print(compatibility_equations[k])
	# 		print('------------------------------')
	# 	compatibility_checker = CompatibilitySystemPDEs(compatibility_equations, phi(x, t), 6)
	# 	depth = compatibility_checker[0]
	# 	status = compatibility_checker[1]
	# 	if status:
	# 		print('The PDE system from WTC Painleve analysis of the KS equation renders compatible')
	# 	else:
	# 		print('The PDE system from WTC Painleve analysis of the KS equation renders incompatible')
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	self.assertFalse(status)

	# def test_TestCaseKawahara(self):
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	print('Here I am testing for the compatibility of the system of PDEs')
	# 	print('that come from the WTC Painleve analysis of the Kawahara equation ....')
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	phi, f, U4 = symbols('phi f U4', cls = Function)
	# 	Kawahara_PDE = parse_expr('diff(f(x, t), t)+f(x, t)*diff(f(x, t), x)+diff(f(x, t), x, x, x)+diff(f(x, t), x, x, x, x, x)')
	# 	# print(integrate(Kawahara_PDE.subs(f(x, t), diff(f(x, t), x)), x))
	# 	alpha, transform, compatibility_equations = ultimatePainleve(
	# 			Kawahara_PDE, # Input PDE
	# 			x, # The input variable used for the spatial variable
	# 			t # The input variable used for the temporal variable
	# 			# alpha # The power balancing exponent
	# 		)
	# 	print('>> len(compatibility_equations)')
	# 	KawaharaPDE_eqn = Kawahara_PDE.subs(f(x, t), transform)
	# 	KawaharaPDE_eqn, _ = fraction(together(KawaharaPDE_eqn.doit()))
	# 	KawaharaPDE_eqn = expand(KawaharaPDE_eqn)
	# 	compatibility_equations = [compatibility_equations[k] for k in range(1, len(compatibility_equations))]
	# 	compatibility_equations.insert(0, KawaharaPDE_eqn)
	# 	compatibility_checker = CompatibilitySystemPDEs(compatibility_equations, phi(x, t), 6)
	# 	depth = compatibility_checker[0]
	# 	status = compatibility_checker[1]
	# 	if status:
	# 		print('The PDE system from WTC Painleve analysis of the Kawahara equation renders compatible')
	# 	else:
	# 		print('The PDE system from WTC Painleve analysis of the Kawahara equation renders incompatible')
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	self.assertFalse(status)

	# def test_TestCaseFisherKPP(self):
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	print('Here I am testing for the compatibility of the system of PDEs')
	# 	print('that come from the WTC Painleve analysis of the Fisher-KPP equation ....')
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	phi, f, U2 = symbols('phi f U2', cls = Function)
	# 	FisherKPP_PDE = parse_expr('diff(f(x, t), t)-diff(f(x, t), x, x)-f(x, t)*(1-f(x, t))')
	# 	alpha, transform, compatibility_equations = ultimatePainleve(
	# 			FisherKPP_PDE, # Input PDE
	# 			x, # The input variable used for the spatial variable
	# 			t # The input variable used for the temporal variable
	# 			# alpha # The power balancing exponent
	# 		)
	# 	FisherKPP_eqn = FisherKPP_PDE.subs(f(x, t), transform)
	# 	FisherKPP_eqn, _ = fraction(together(FisherKPP_eqn.doit()))
	# 	FisherKPP_eqn = expand(FisherKPP_eqn)
	# 	compatibility_equations = [compatibility_equations[k] for k in range(1, len(compatibility_equations))]
	# 	compatibility_equations.insert(0, FisherKPP_eqn)
	# 	compatibility_checker = CompatibilitySystemPDEs(compatibility_equations, phi(x, t), 6)
	# 	depth = compatibility_checker[0]
	# 	status = compatibility_checker[1]
	# 	if status:
	# 		print('The PDE system from WTC Painleve analysis of the Fisher-KPP equation renders compatible')
	# 	else:
	# 		print('The PDE system from WTC Painleve analysis of the Fisher-KPP equation renders incompatible')
	# 	print('------------------------------------------------------------------------------------------------------------------------------------')
	# 	self.assertFalse(status)

	def test_TestCaseFitzHughNagumo(self):
		print('------------------------------------------------------------------------------------------------------------------------------------')
		print('Here I am testing for the compatibility of the system of PDEs')
		print('that come from the WTC Painleve analysis of the FitzHugh-Nagumo equation ....')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		phi, f, U1 = symbols('phi f U1', cls = Function)
		FitzHughNagumo_PDE = parse_expr('diff(f(x, t), t)-diff(f(x, t), x, x)+f(x, t)*(1-f(x, t))*(2-f(x, t))')
		alpha, transform, compatibility_equations = ultimatePainleve(
				FitzHughNagumo_PDE, # Input PDE
				x, # The input variable used for the spatial variable
				t # The input variable used for the temporal variable
				# alpha # The power balancing exponent
			)
		FitzHughNagumoPDE_eqn = FitzHughNagumo_PDE.subs(f(x, t), transform)
		FitzHughNagumoPDE_eqn, _ = fraction(together(FitzHughNagumoPDE_eqn.doit()))
		FitzHughNagumoPDE_eqn = expand(FitzHughNagumoPDE_eqn)
		compatibility_equations = [compatibility_equations[k] for k in range(1, len(compatibility_equations))]
		compatibility_equations.insert(0, FitzHughNagumoPDE_eqn)
		compatibility_checker = CompatibilitySystemPDEs(compatibility_equations, phi(x, t), 6)
		depth = compatibility_checker[0]
		status = compatibility_checker[1]
		if status:
			print('The PDE system from WTC Painleve analysis of the FitzHugh-Nagumo equation renders compatible')
		else:
			print('The PDE system from WTC Painleve analysis of the FitzHugh-Nagumo equation renders incompatible')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		self.assertFalse(status)

unittest.main(argv=[''], verbosity=2, exit=False)
