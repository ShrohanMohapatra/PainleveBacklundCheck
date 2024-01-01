from sympy import *
from collections import Counter

def showDiffNotation(expr):
	functions = expr.atoms(Function)
	reps = {}
	
	for fun in functions:
		# Consider the case that some functions won't have the name
		# attribute e.g. Abs of an elementary function
		try:            
			reps[fun] = Symbol(fun.name) # Otherwise functions with greek symbols aren't replaced
		except AttributeError:
			continue

	dreps = [(deriv, Symbol(deriv.expr.subs(reps).name + "_{," +
							''.join(par.name for par in deriv.variables) + "}"))  \
			for deriv in expr.atoms(Derivative)]

	# Ensure that higher order derivatives are replaced first, then lower ones.
	# Otherwise you get d/dr w_r instead of w_rr

	dreps.sort(key=lambda x: len(x[0].variables), reverse=True)
	output = expr.subs(dreps).subs(reps)

	return output

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

def CompatibilitySystemPDEs(systemsOfPDEs, func, depth_max = 150, depth = 0):

	print('------------------------------------------------------------------------------------------------')
	print('Welcome to recursion number = ',depth+1)
	print('------------------------------------------------------------------------------------------------')

	for k in range(len(systemsOfPDEs)):
		print('--------------------------------------------------------------------------------------------')
		print('>> PDE Number '+str(k+1)+' -->')
		print('--------------------------------------------------------------------------------------------')
		print(systemsOfPDEs[k])
		print('--------------------------------------------------------------------------------------------')
		print('>> Latex version')
		print('--------------------------------------------------------------------------------------------')
		try:
			print_latex(showDiffNotation(systemsOfPDEs[k]))
		except:
			raise
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
			print('Latex version')
			print('---------------------------------------------------')
			print_latex(showDiffNotation(systemsOfPDEs_F[k]))
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

x, t, sigma, b, kappa = symbols('x t sigma b kappa')
phi = Function('phi')(x, t)
f = Function('f')

TransformedSineGordonEquation = parse_expr('2*f(x,t)*diff(f(x,t),x,t)-2*diff(f(x,t),x)*diff(f(x,t),t)-f(x,t)**3+f(x,t)')

u0 = Function('u0')(x, t)
u1 = Function('u1')(x, t)
u2 = Function('u2')(x, t)

u = u0/phi**2 + u1/phi + u2
TransformedSineGordonPainleveExpansion = expand(TransformedSineGordonEquation.subs(f(x, t), u).doit())

generator_list = [
	phi, 1/phi,
	diff(phi, x), diff(phi, t),
	diff(phi, x, x), diff(phi, x, t), diff(phi, t, t),
	diff(phi, x, x, x), diff(phi, x, x, t), diff(phi, x, t, t), diff(phi, t, t, t),
	diff(phi, x, x, x, x), diff(phi, x, x, x, t), diff(phi, x, x, t, t), diff(phi, x, t, t, t), diff(phi, t, t, t, t),
	diff(phi, x, x, x, x, x), diff(phi, x, x, x, x, t), diff(phi, x, x, x, t, t), diff(phi, x, x, t, t, t), diff(phi, x, t, t, t, t), diff(phi, t, t, t, t, t)
	]

TransformedSineGordonPainleveExpansion_polyInTermsOfPhi = Poly(TransformedSineGordonPainleveExpansion, gens = generator_list)
TransformedSineGordonPainleveExpansion_polyInTermsOfPhi = TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.as_expr()

Equation0 = simplify(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 6))
factors_list_of_equation0 = factor_list(Equation0)
Equation0 = factors_list_of_equation0[-1][-1][0]
u0_solution = solve(Equation0, u0)[0]
Equation0_residual = Equation0.subs(u0, u0_solution)
print('--------------------------------------------------------------------------')
print('u0 = ', u0_solution)
print('--------------------------------------------------------------------------')
print('Equation0_residual = ', Equation0_residual)
print('--------------------------------------------------------------------------')

Equation1 = expand(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 5))
Equation1 = expand(Equation1.subs(u0, u0_solution).doit())
factors_list_of_equation1 = factor_list(Equation1)
Equation1 = factors_list_of_equation1[-1][-1][0]
Equation1Storage = Equation1
u1_solution = solve(Equation1, u1)[0]
u1_solution = together(u1_solution).doit()
Equation1_residual = Equation1.subs(u1, u1_solution)
Equation1_residual, _ = fraction(together(Equation1_residual).doit())
Equation1_residual = expand(Equation1_residual)
print('--------------------------------------------------------------------------')
print('Equation for 1st power coefficient')
print('--------------------------------------------------------------------------')
print('u1 = ', u1_solution)
print('--------------------------------------------------------------------------')
print('Equation1_residual = ', Equation1_residual)
print('--------------------------------------------------------------------------')

Equation2 = expand(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 4))
Equation2 = expand(Equation2.subs(u0, u0_solution).doit())
Equation2Storage = Equation2
Equation2 = expand(Equation2.subs(u1, u1_solution).doit())
Equation2, _ = fraction(together(Equation2).doit())
Equation2 = expand(Equation2)
factors_list_of_equation2 = factor_list(Equation2)
# print('--------------------------------------------------------------------------')
# print(factors_list_of_equation2[-1][0][0])
# print('--------------------------------------------------------------------------')
# Equation2 = factors_list_of_equation2[-1][-1][0]
# u2_solution = solve(Equation2, u2)[0]
# Equation2_residual = Equation2.subs(u2, u2_solution)
# print('--------------------------------------------------------------------------')
# print('u2 = ', u2_solution)
# print('--------------------------------------------------------------------------')
print('Equation for 2nd power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation2Storage),' = 0')
print('--------------------------------------------------------------------------')

Equation3 = expand(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 3))
Equation3 = expand(Equation3.subs(u0, u0_solution).doit())
Equation3 = expand(Equation3.subs(u1, u1_solution).doit())
Equation3_compatibility_test = Equation3
# Equation3 = expand(Equation3.subs(u2, u2_solution).doit())
# Equation3, _ = fraction(together(Equation3).doit())
# Equation3 = expand(Equation3)
# factors_list_of_equation3 = factor_list(Equation3)
# Equation3 = factors_list_of_equation3[-1][-1][0]
print('--------------------------------------------------------------------------')
print('Equation for 3rd power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation3_compatibility_test), '= 0')
print('--------------------------------------------------------------------------')

Equation4 = expand(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 2))
Equation4 = expand(Equation4.subs(u0, u0_solution).doit())
Equation4 = expand(Equation4.subs(u1, u1_solution).doit())
Equation4_compatibility_test = Equation4
# Equation4 = expand(Equation4.subs(u2, u2_solution).doit())
# Equation4, _ = fraction(together(Equation4).doit())
# Equation4 = expand(Equation4)
# factors_list_of_equation4 = factor_list(Equation4)
# Equation4 = factors_list_of_equation4[-1][-1][0]
print('--------------------------------------------------------------------------')
print('Equation for 4th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation4_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation5 = expand(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 1))
Equation5 = expand(Equation5.subs(u0, u0_solution).doit())
Equation5 = expand(Equation5.subs(u1, u1_solution).doit())
Equation5_compatibility_test = Equation5
# Equation5 = expand(Equation5.subs(u2, u2_solution).doit())
# Equation5, _ = fraction(together(Equation5).doit())
# Equation5 = expand(Equation5)
# factors_list_of_equation5 = factor_list(Equation5)
# print('--------------------------------------------------------------------------')
# print(len(factors_list_of_equation5))
# print('--------------------------------------------------------------------------')
# Equation5 = factors_list_of_equation5[-1][-1][0]
print('--------------------------------------------------------------------------')
print('Equation for 5th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation5_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation6 = expand(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(phi, 0))
Equation6 = expand(Equation6.subs(u0, u0_solution).doit())
Equation6 = expand(Equation6.subs(u1, u1_solution).doit())
Equation6_compatibility_test = Equation6
# Equation5 = expand(Equation5.subs(u2, u2_solution).doit())
# Equation5, _ = fraction(together(Equation5).doit())
# Equation5 = expand(Equation5)
# factors_list_of_equation6 = factor_list(Equation6)
# print('--------------------------------------------------------------------------')
# print(len(factors_list_of_equation6))
# print('--------------------------------------------------------------------------')
# Equation6 = factors_list_of_equation6[-1][-1][0]
print('--------------------------------------------------------------------------')
print('Equation for 6th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation6_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

phi_ttxx_solution = expand(solve(Equation3_compatibility_test, diff(phi,t,t,x,x))[0])
print('--------------------------------------------------------------------------')
print('Solving equation 3 -> phi_ttxx = ', showDiffNotation(phi_ttxx_solution))
print('--------------------------------------------------------------------------')

u2_xt_solution = expand(solve(Equation6_compatibility_test, diff(u2,t,x))[0])
print('--------------------------------------------------------------------------')
print('Solving equation 6 -> u2_{ x, t} = ', showDiffNotation(u2_xt_solution))
print('--------------------------------------------------------------------------')

Equation4_compatibility_test = expand(Equation4_compatibility_test.subs(
	diff(phi,t,t,x,x), phi_ttxx_solution
	).doit())
Equation4_compatibility_test, _ = fraction(together(Equation4_compatibility_test).doit())
Equation4_compatibility_test = expand(Equation4_compatibility_test)
print('--------------------------------------------------------------------------')
print('Equation for 4th power coefficient using the solution to phi_{ x, x, t, t}')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation4_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation5_compatibility_test = expand(Equation5_compatibility_test.subs(
	diff(phi,t,t,x,x), phi_ttxx_solution
	).doit())
Equation5_compatibility_test, _ = fraction(together(Equation5_compatibility_test).doit())
Equation5_compatibility_test = expand(Equation5_compatibility_test)
print('--------------------------------------------------------------------------')
print('Equation for 5th power coefficient using the solution to phi_{ x, x, t, t}')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation5_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation4_compatibility_test = expand(Equation4_compatibility_test.subs(
	diff(u2,x,t), u2_xt_solution
	).doit())
Equation4_compatibility_test, _ = fraction(together(Equation4_compatibility_test).doit())
Equation4_compatibility_test = expand(Equation4_compatibility_test)
print('--------------------------------------------------------------------------')
print('Equation for 4th power coefficient using the solution to u2_{ x, t}')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation4_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation5_compatibility_test = expand(Equation5_compatibility_test.subs(
	diff(u2,x,t), u2_xt_solution
	).doit())
Equation5_compatibility_test, _ = fraction(together(Equation5_compatibility_test).doit())
Equation5_compatibility_test = expand(Equation5_compatibility_test)
print('--------------------------------------------------------------------------')
print('Equation for 5th power coefficient using the solution to phi_{ x, x, t, t}')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation5_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

u2_xt_solution = expand(solve(Equation6_compatibility_test, diff(u2,t,x))[0])
print('--------------------------------------------------------------------------')
print('Solving equation 6 -> u2_{ x, t} = ', showDiffNotation(u2_xt_solution))
print('--------------------------------------------------------------------------')

Equation7_compatibility_test = expand(resultant(
	Equation4_compatibility_test,
	Equation5_compatibility_test,
	diff(phi, t, t, x)))
print('--------------------------------------------------------------------------')
print('Equation7_compatibility_test (polynomial in terms of diff(phi, t, x, x))')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation7_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation8_compatibility_test = expand(resultant(
	Equation4_compatibility_test,
	Equation5_compatibility_test,
	diff(phi, t, x, x)))
print('--------------------------------------------------------------------------')
print('Equation8_compatibility_test (polynomial in terms of diff(phi, t, t, x))')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation8_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

phi_txx_solution = solve(Equation7_compatibility_test, diff(phi, t, x, x))[0]
print('--------------------------------------------------------------------------')
print('phi_txx_solution = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_txx_solution))
print('--------------------------------------------------------------------------')
print_latex(phi_txx_solution)
print('--------------------------------------------------------------------------')

phi_ttx_solution = solve(expand(Equation4_compatibility_test.subs(diff(phi,t,x,x),phi_txx_solution).doit()),diff(phi,t,t,x))[-1]
print('--------------------------------------------------------------------------')
print('phi_ttx_solution = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_ttx_solution))
print('--------------------------------------------------------------------------')
print_latex(phi_ttx_solution)
print('--------------------------------------------------------------------------')

phi_ttxx_solution = expand(
	phi_ttxx_solution.subs(
		diff(phi,t,x,x),phi_txx_solution
		).subs(
		diff(phi,t,t,x),phi_ttx_solution
		).doit())

print('--------------------------------------------------------------------------')
print('phi_ttxx_solution = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_ttxx_solution))
print('--------------------------------------------------------------------------')
print_latex(phi_ttxx_solution)
print('--------------------------------------------------------------------------')

dictionary_of_solutions = {}
dictionary_of_solutions[diff(phi,x,x,t)] = phi_txx_solution
dictionary_of_solutions[diff(phi,x,t,t)] = phi_ttx_solution
dictionary_of_solutions[diff(u2,x,t)] = u2_xt_solution

phi_ttxx_solution_1 = expand(diff(dictionary_of_solutions[diff(phi,x,x,t)],t).doit())
for term in dictionary_of_solutions:
	phi_ttxx_solution_1 = phi_ttxx_solution_1.subs(term, dictionary_of_solutions[term]).doit()
	phi_ttxx_solution_1 = expand(phi_ttxx_solution_1).doit()
print('--------------------------------------------------------------------------')
print('phi_ttxx_solution_1 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_ttxx_solution_1))
print('--------------------------------------------------------------------------')
print_latex(phi_ttxx_solution_1)
print('--------------------------------------------------------------------------')

phi_ttxx_solution_2 = expand(diff(dictionary_of_solutions[diff(phi,x,t,t)],x).doit())
for term in dictionary_of_solutions:
	phi_ttxx_solution_2 = phi_ttxx_solution_2.subs(term, dictionary_of_solutions[term]).doit()
	phi_ttxx_solution_2 = expand(phi_ttxx_solution_2).doit()
print('--------------------------------------------------------------------------')
print('phi_ttxx_solution_2 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_ttxx_solution_2))
print('--------------------------------------------------------------------------')
print_latex(phi_ttxx_solution_2)
print('--------------------------------------------------------------------------')

compatibility_check_1 = phi_ttxx_solution_1 - phi_ttxx_solution_2
compatibility_check_1, _ = fraction(together(compatibility_check_1).doit())
compatibility_check_1 = expand(compatibility_check_1).doit()
print('--------------------------------------------------------------------------')
print('phi_ttxx_solution_1 - phi_ttxx_solution_2 = ')
print('--------------------------------------------------------------------------')
print(compatibility_check_1)
print('--------------------------------------------------------------------------')
print_latex(compatibility_check_1)
print('--------------------------------------------------------------------------')

compatibility_check_2 = phi_ttxx_solution_1 - phi_ttxx_solution
compatibility_check_2, _ = fraction(together(compatibility_check_2).doit())
compatibility_check_2 = expand(compatibility_check_2).doit()
print('--------------------------------------------------------------------------')
print('phi_ttxx_solution_1 - phi_ttxx_solution = ')
print('--------------------------------------------------------------------------')
print(compatibility_check_2)
print('--------------------------------------------------------------------------')
print_latex(compatibility_check_2)
print('--------------------------------------------------------------------------')

compatibility_check_3 = phi_ttxx_solution_2 - phi_ttxx_solution
compatibility_check_3, _ = fraction(together(compatibility_check_3).doit())
compatibility_check_3 = expand(compatibility_check_3).doit()
print('--------------------------------------------------------------------------')
print('phi_ttxx_solution_2 - phi_ttxx_solution = ')
print('--------------------------------------------------------------------------')
print(compatibility_check_3)
print('--------------------------------------------------------------------------')
print_latex(compatibility_check_3)
print('--------------------------------------------------------------------------')

compatibility_check_1_new = diff(phi,x,t)**2-u2*diff(phi,x)*diff(phi,t) - (solve(compatibility_check_1,sqrt(diff(phi,x,t)**2-u2*diff(phi,x)*diff(phi,t)))[0])**2
compatibility_check_1_new, _ = fraction(together(compatibility_check_1_new).doit())
compatibility_check_1_new = expand(compatibility_check_1_new)

compatibility_check_2_new = diff(phi,x,t)**2-u2*diff(phi,x)*diff(phi,t) - (solve(compatibility_check_2,sqrt(diff(phi,x,t)**2-u2*diff(phi,x)*diff(phi,t)))[0])**2
compatibility_check_2_new, _ = fraction(together(compatibility_check_2_new).doit())
compatibility_check_2_new = expand(compatibility_check_2_new)

compatibility_check_3_new = diff(phi,x,t)**2-u2*diff(phi,x)*diff(phi,t) - (solve(compatibility_check_3,sqrt(diff(phi,x,t)**2-u2*diff(phi,x)*diff(phi,t)))[0])**2
compatibility_check_3_new, _ = fraction(together(compatibility_check_3_new).doit())
compatibility_check_3_new = expand(compatibility_check_3_new)

compatibility_check_wrap_1 = expand(compatibility_check_1_new + compatibility_check_2_new**2)
print('--------------------------------------------------------------------------')
print('compatibility_check_wrap_1 :')
print('--------------------------------------------------------------------------')
print(showDiffNotation(compatibility_check_wrap_1))
print('--------------------------------------------------------------------------')

compatibility_check_wrap_2 = expand(compatibility_check_2_new**5 - compatibility_check_3_new**3)
print('--------------------------------------------------------------------------')
print('compatibility_check_wrap_2 :')
print('--------------------------------------------------------------------------')
print(showDiffNotation(compatibility_check_wrap_2))
print('--------------------------------------------------------------------------')

compatibility_check_wrap_3 = expand(compatibility_check_1_new**5 + compatibility_check_3_new**6)
print('--------------------------------------------------------------------------')
print('compatibility_check_wrap_3 :')
print('--------------------------------------------------------------------------')
print(showDiffNotation(compatibility_check_wrap_3))
print('--------------------------------------------------------------------------')
