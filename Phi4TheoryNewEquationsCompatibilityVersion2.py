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

Phi4TheoryEquation = parse_expr('diff(f(x,t),t,t)-diff(f(x,t),x,x)+2*f(x,t)-2*f(x,t)**3')

u0 = Function('u0')(x, t)
u1 = Function('u1')(x, t)

u = u0/phi + u1
Phi4TheoryPainleveExpansion = expand(Phi4TheoryEquation.subs(f(x, t), u).doit())

generator_list = [
	phi, 1/phi,
	diff(phi, x), diff(phi, t),
	diff(phi, x, x), diff(phi, x, t), diff(phi, t, t),
	diff(phi, x, x, x), diff(phi, x, x, t), diff(phi, x, t, t), diff(phi, t, t, t),
	diff(phi, x, x, x, x), diff(phi, x, x, x, t), diff(phi, x, x, t, t), diff(phi, x, t, t, t), diff(phi, t, t, t, t),
	diff(phi, x, x, x, x, x), diff(phi, x, x, x, x, t), diff(phi, x, x, x, t, t), diff(phi, x, x, t, t, t), diff(phi, x, t, t, t, t), diff(phi, t, t, t, t, t)
	]

Phi4TheoryPainleveExpansion_polyInTermsOfPhi = Poly(Phi4TheoryPainleveExpansion, gens = generator_list)
Phi4TheoryPainleveExpansion_polyInTermsOfPhi = Phi4TheoryPainleveExpansion_polyInTermsOfPhi.as_expr()

Equation0 = simplify(Phi4TheoryPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 3))
factors_list_of_equation0 = factor_list(Equation0)
Equation0 = factors_list_of_equation0[-1][-1][0]
u0_solution = solve(Equation0, u0)[0]
Equation0_residual = Equation0.subs(u0, u0_solution)
# print('--------------------------------------------------------------------------')
# print('u0 = ', u0_solution)
# print('--------------------------------------------------------------------------')
# print('Equation0_residual = ', Equation0_residual)
# print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')
print('Equation for 0th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation0),' = 0')
print('--------------------------------------------------------------------------')

Equation1 = expand(Phi4TheoryPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 2))
# Equation1 = expand(Equation1.subs(u0, u0_solution).doit())
factors_list_of_equation1 = factor_list(Equation1)
Equation1 = factors_list_of_equation1[-1][-1][0]
Equation1Storage = Equation1
# u1_solution = solve(Equation1, u1)[0]
# u1_solution = together(u1_solution).doit()
# Equation1_residual = Equation1.subs(u1, u1_solution)
# Equation1_residual, _ = fraction(together(Equation1_residual).doit())
# Equation1_residual = expand(Equation1_residual)
# print('--------------------------------------------------------------------------')
# print('Equation for 1st power coefficient')
# print('--------------------------------------------------------------------------')
# print('u1 = ', u1_solution)
print('--------------------------------------------------------------------------')
print('Equation for 1st power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation1Storage),' = 0')
print('--------------------------------------------------------------------------')

Equation2 = expand(Phi4TheoryPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 1))
# Equation2 = expand(Equation2.subs(u0, u0_solution).doit())
Equation2Storage = Equation2
# Equation2 = expand(Equation2.subs(u1, u1_solution).doit())
# Equation2, _ = fraction(together(Equation2).doit())
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
print('--------------------------------------------------------------------------')
print('Equation for 2nd power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation2Storage),' = 0')
print('--------------------------------------------------------------------------')

Equation3 = expand(Phi4TheoryPainleveExpansion_polyInTermsOfPhi.coeff(phi, 0))
# Equation3 = expand(Equation3.subs(u0, u0_solution).doit())
# Equation3 = expand(Equation3.subs(u1, u1_solution).doit())
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

print(
	CompatibilitySystemPDEs([
		Equation0,
		Equation1Storage,
		Equation2Storage,
		Equation3_compatibility_test
		], phi)
)

phi_x_solution = solve(Equation0, diff(phi,x))[-1]
print('--------------------------------------------------------------------------')
print('phi_x_solution = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_x_solution))
print('--------------------------------------------------------------------------')

phi_xx_solution = diff(phi_x_solution, x)
phi_xx_solution = expand(phi_xx_solution.doit())

phi_xt_solution = diff(phi_x_solution, t)
phi_xt_solution = expand(phi_xt_solution.doit())

phi_tt_solution = solve(Equation1Storage, diff(phi,t,t))[-1]
phi_tt_solution = expand(phi_tt_solution.subs(diff(phi,x),phi_x_solution).subs(diff(phi,x,x),phi_xx_solution).doit())

solutions = solve([expand(phi_xt_solution-diff(phi,x,t)),expand(phi_tt_solution-diff(phi,t,t))],[diff(phi,x,t),diff(phi,t,t)])
phi_xt_solution = solutions[diff(phi,x,t)]
phi_tt_solution = solutions[diff(phi,t,t)]

print('--------------------------------------------------------------------------')
print('phi_xt_solution = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_xt_solution))
print('--------------------------------------------------------------------------')

print('--------------------------------------------------------------------------')
print('phi_tt_solution = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_tt_solution))
print('--------------------------------------------------------------------------')

phi_xx_solution = expand(phi_xx_solution.subs(diff(phi, x,t),phi_xt_solution).doit())
print('--------------------------------------------------------------------------')
print('phi_xx_solution = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_xx_solution))
print('--------------------------------------------------------------------------')


u0_tt_solution = solve(Equation2Storage, diff(u0,t,t))[-1]
print('--------------------------------------------------------------------------')
print('u0_tt_solution = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(u0_tt_solution))
print('--------------------------------------------------------------------------')

u1_tt_solution = solve(Equation3_compatibility_test, diff(u1,t,t))[-1]
print('--------------------------------------------------------------------------')
print('u1_tt_solution = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(u1_tt_solution))
print('--------------------------------------------------------------------------')

dictionary_of_solutions = {}
dictionary_of_solutions[diff(phi,x)] = phi_x_solution
dictionary_of_solutions[diff(phi,x,t)] = phi_xt_solution
dictionary_of_solutions[diff(phi,x,x)] = phi_xx_solution
dictionary_of_solutions[diff(phi,t,t)] = phi_tt_solution
dictionary_of_solutions[diff(u0,t,t)] = u0_tt_solution
dictionary_of_solutions[diff(u1,t,t)] = u1_tt_solution

phi_xtt_solution_1 = expand(diff(dictionary_of_solutions[diff(phi,x,t)],t).doit())
for term in dictionary_of_solutions:
	phi_xtt_solution_1 = phi_xtt_solution_1.subs(term, dictionary_of_solutions[term]).doit()
	phi_xtt_solution_1 = expand(phi_xtt_solution_1)
print('--------------------------------------------------------------------------')
print('phi_xtt_solution_1 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_xtt_solution_1))
print('--------------------------------------------------------------------------')

phi_xtt_solution_2 = expand(diff(dictionary_of_solutions[diff(phi,t,t)],x).doit())
for term in dictionary_of_solutions:
	phi_xtt_solution_2 = phi_xtt_solution_2.subs(term, dictionary_of_solutions[term]).doit()
	phi_xtt_solution_2 = expand(phi_xtt_solution_2)
print('--------------------------------------------------------------------------')
print('phi_xtt_solution_2 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_xtt_solution_2))
print('--------------------------------------------------------------------------')

compatibility_check_1 = phi_xtt_solution_1 - phi_xtt_solution_2
compatibility_check_1, _ = fraction(together(compatibility_check_1).doit())
compatibility_check_1 = expand(compatibility_check_1.doit())
print('--------------------------------------------------------------------------')
print('phi_xtt_solution_1 - phi_xtt_solution_2 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(compatibility_check_1))
print('--------------------------------------------------------------------------')

phi_xxt_solution_1 = expand(diff(dictionary_of_solutions[diff(phi,x,x)],t).doit())
for term in dictionary_of_solutions:
	phi_xxt_solution_1 = phi_xtt_solution_1.subs(term, dictionary_of_solutions[term]).doit()
	phi_xxt_solution_1 = expand(phi_xtt_solution_1)
print('--------------------------------------------------------------------------')
print('phi_xxt_solution_1 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_xxt_solution_1))
print('--------------------------------------------------------------------------')

phi_xxt_solution_2 = expand(diff(dictionary_of_solutions[diff(phi,x,t)],x).doit())
for term in dictionary_of_solutions:
	phi_xxt_solution_2 = phi_xxt_solution_2.subs(term, dictionary_of_solutions[term]).doit()
	phi_xxt_solution_2 = expand(phi_xxt_solution_2)
print('--------------------------------------------------------------------------')
print('phi_xxt_solution_2 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_xxt_solution_2))
print('--------------------------------------------------------------------------')

compatibility_check_2 = phi_xxt_solution_1 - phi_xxt_solution_2
compatibility_check_2, _ = fraction(together(compatibility_check_2).doit())
compatibility_check_2 = expand(compatibility_check_2.doit())
print('--------------------------------------------------------------------------')
print('phi_xxt_solution_1 - phi_xxt_solution_2 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(compatibility_check_2))
print('--------------------------------------------------------------------------')

# A1 = compatibility_check_1.coeff(sqrt(diff(phi, t)**2-u0**2))
# B1 = compatibility_check_1 - sqrt(diff(phi, t)**2-u0**2)*compatibility_check_1.coeff(sqrt(diff(phi, t)**2-u0**2))
# A2 = compatibility_check_2.coeff(sqrt(diff(phi, t)**2-u0**2))
# B2 = compatibility_check_2 - sqrt(diff(phi, t)**2-u0**2)*compatibility_check_2.coeff(sqrt(diff(phi, t)**2-u0**2))

# print('--------------------------------------------------------------------------')
# print('A1 = ')
# print('--------------------------------------------------------------------------')
# print(showDiffNotation(A1))
# print('--------------------------------------------------------------------------')

# print('--------------------------------------------------------------------------')
# print('B1 = ')
# print('--------------------------------------------------------------------------')
# print(showDiffNotation(B1))
# print('--------------------------------------------------------------------------')

# print('--------------------------------------------------------------------------')
# print('A2 = ')
# print('--------------------------------------------------------------------------')
# print(showDiffNotation(A2))
# print('--------------------------------------------------------------------------')

# print('--------------------------------------------------------------------------')
# print('B2 = ')
# print('--------------------------------------------------------------------------')
# print(showDiffNotation(B2))
# print('--------------------------------------------------------------------------')

# d1 = diff(phi, t)**2-u0**2-B1**2/A1**2
# d1, _ = fraction(together(d1).doit())
# d1 = expand(d1)
# print('--------------------------------------------------------------------------')
# print('d1 = ')
# print('--------------------------------------------------------------------------')
# print(showDiffNotation(d1))
# print('--------------------------------------------------------------------------')

# d2 = diff(phi, t)**2-u0**2-B2**2/A2**2
# d2, _ = fraction(together(d2).doit())
# d2 = expand(d2)
# print('--------------------------------------------------------------------------')
# print('d2 = ')
# print('--------------------------------------------------------------------------')
# print(showDiffNotation(d2))
# print('--------------------------------------------------------------------------')

# generator_list = [
# 	phi, 1/phi, u0, u1,
# 	diff(phi, x), diff(phi, t),
# 	diff(phi, x, x), diff(phi, x, t), diff(phi, t, t),
# 	diff(phi, x, x, x), diff(phi, x, x, t), diff(phi, x, t, t), diff(phi, t, t, t),
# 	diff(phi, x, x, x, x), diff(phi, x, x, x, t), diff(phi, x, x, t, t), diff(phi, x, t, t, t), diff(phi, t, t, t, t),
# 	diff(phi, x, x, x, x, x), diff(phi, x, x, x, x, t), diff(phi, x, x, x, t, t), diff(phi, x, x, t, t, t), diff(phi, x, t, t, t, t), diff(phi, t, t, t, t, t),
# 	diff(u0, x), diff(u0, t), diff(u0, x, x), diff(u0, x, t), diff(u0, t, t),
# 	diff(u1, x), diff(u1, t), diff(u1, x, x), diff(u1, x, t), diff(u1, t, t)
# 	]

# d1 = Poly(d1, gens = generator_list)
# d2 = Poly(d2, gens = generator_list)

# print('--------------------------------------------------------------------------')
# print('resultant(d1, d2) = ')
# print('--------------------------------------------------------------------------')
# print(resultant(d1, d2))
# print('--------------------------------------------------------------------------')
