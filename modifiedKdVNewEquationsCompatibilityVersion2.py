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


	# print('----------------------------------------')
	# print('Here we compute the commutation of the following two PDEs:')
	# print('----------------------------------------')
	# print('Eqn (1) >>')
	# print('----------------------------------------')
	# print(inputPDE1, '=', 0)
	# print('----------------------------------------')
	# print('----------------------------------------')
	# print('Eqn (2) >>')
	# print('----------------------------------------')
	# print(inputPDE2, '=', 0)
	# print('----------------------------------------')
	# print('The maximum order term in equation 1 is gamma_1 = ', max_order_term1[2])
	# print('----------------------------------------')
	# print('----------------------------------------')
	# print('The maximum order term in equation 2 is gamma_2 = ', max_order_term2[2])
	print('----------------------------------------')
	print('In this scheme we represent eqn1 = A1*gamma_1 + B1 and eqn2 = A2*gamma_2 + B2 where')
	print('----------------------------------------')
	print('A1 =', A1)
	print('----------------------------------------')
	print('A2 =', A2)
	print('----------------------------------------')
	# print('B1 =', B1)
	# print('----------------------------------------')
	# print('B2 =', B2)
	# print('----------------------------------------')

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

	# print('-------------------------------')
	# print('The variables of differentiation for the expression -B1/A1 from equation 1 are:')
	# print('-------------------------------')
	# print(varsFromPDE3)
	# print('-------------------------------')
	# print('The variables of differentiation for the expression -B2/A2 from equation 1 are:')
	# print('-------------------------------')
	# print('varsFromPDE4')
	# print('-------------------------------')
	# print(varsFromPDE4)
	# print('-------------------------------')

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

	# print('----------------------------')
	# print('Expression after differentiation of equation 1')
	# print('----------------------------')
	# print(inputPDE3)
	# print('----------------------------')
	# print('----------------------------')
	# print('----------------------------')
	# print('Expression after differentiation of equation 2')
	# print('----------------------------')
	# print(inputPDE4)
	# print('----------------------------')
	# print('----------------------------')

	expression = inputPDE3 - inputPDE4
	expression, _ = fraction(together(expression.doit()))
	expression = expand(expression)

	# print('----------------------------')
	# print('Expression after commutation')
	# print('----------------------------')
	# print(expression)
	# print('----------------------------')

	return expression

def CompatibilitySystemPDEs(systemsOfPDEs, func, depth_max = 150, depth = 0):

	print('------------------------------------------------------------------------------------------------')
	print('Welcome to recursion number = ',depth+1)
	print('------------------------------------------------------------------------------------------------')

	# for k in range(len(systemsOfPDEs)):
	# 	print('--------------------------------------------------------------------------------------------')
	# 	print('>> PDE Number '+str(k+1)+' -->')
	# 	print('--------------------------------------------------------------------------------------------')
	# 	print(systemsOfPDEs[k])
	# 	print('--------------------------------------------------------------------------------------------')
		
		# print('>> Latex version')
		# print('--------------------------------------------------------------------------------------------')
		# try:
		# 	print_latex(showDiffNotation(systemsOfPDEs[k]))
		# except:
		# 	raise
		# print('--------------------------------------------------------------------------------------------')

	orders_list = [max_order_term(equation, func)[3] for equation in systemsOfPDEs]
	if len(systemsOfPDEs) == 2:
		maximum_order = max(orders_list)
		maximum_order_index = orders_list.index(maximum_order)
		equation_2_new = expand(simplificationPDE(systemsOfPDEs[0], systemsOfPDEs[1], func))
		systemsOfPDEs_F = [equation_2_new, systemsOfPDEs[maximum_order_index]]
		for k in range(len(systemsOfPDEs_F)):
			print('---------------------------------------------------')
			print('Number of terms in ',k, 'th PDE >>', len(Add.make_args(systemsOfPDEs_F[k])))
			print('---------------------------------------------------')
		if equation_2_new == 0:
			return [depth, True]
		else:
			if depth >= depth_max: return [depth, False]
			try:
				return CompatibilitySystemPDEs([equation_2_new, systemsOfPDEs[maximum_order_index]], func, depth_max, depth + 1)
			except Exception:
				return [depth+1, False]
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
			# print(systemsOfPDEs_F[k])
			print('---------------------------------------------------')
			# print('Latex version')
			# print('---------------------------------------------------')
			# print_latex(showDiffNotation(systemsOfPDEs_F[k]))
			# print('---------------------------------------------------')
			flag = flag or systemsOfPDEs_F[k].equals(0)
		if flag: return [depth, flag]
		else:
			if depth >= depth_max: return [depth+1, False]
			try:
				return CompatibilitySystemPDEs(systemsOfPDEs_F, func, depth_max, depth + 1)
			except Exception:
				return [depth+1, False]

x, t, sigma, b, kappa = symbols('x t sigma b kappa')
phi = Function('phi')(x, t)
f = Function('f')

modifiedKdVEquation = parse_expr('diff(f(x,t),t)+2*sigma**2*diff(f(x,t),x,x,x)-3*f(x,t)**2*diff(f(x,t),x)')

u0 = Function('u0')(x, t)
u1 = Function('u1')(x, t)

u = u0/phi + u1
modifiedKdVPainleveExpansion = expand(modifiedKdVEquation.subs(f(x, t), u).doit())

generator_list = [
	phi, 1/phi,
	diff(phi, x), diff(phi, t),
	diff(phi, x, x), diff(phi, x, t), diff(phi, t, t),
	diff(phi, x, x, x), diff(phi, x, x, t), diff(phi, x, t, t), diff(phi, t, t, t),
	diff(phi, x, x, x, x), diff(phi, x, x, x, t), diff(phi, x, x, t, t), diff(phi, x, t, t, t), diff(phi, t, t, t, t),
	diff(phi, x, x, x, x, x), diff(phi, x, x, x, x, t), diff(phi, x, x, x, t, t), diff(phi, x, x, t, t, t), diff(phi, x, t, t, t, t), diff(phi, t, t, t, t, t)
	]

modifiedKdVPainleveExpansion_polyInTermsOfPhi = Poly(modifiedKdVPainleveExpansion, gens = generator_list)
modifiedKdVPainleveExpansion_polyInTermsOfPhi = modifiedKdVPainleveExpansion_polyInTermsOfPhi.as_expr()
Equation0 = simplify(modifiedKdVPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 4))
print('--------------------------------------------------------------------------')
print('P0 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation0), ' = 0')
print('--------------------------------------------------------------------------')
factors_list_of_equation0 = factor_list(Equation0)
# print('--------------------------------------------------------------------------')
# print('factors_list_of_equation0 = ')
# print('--------------------------------------------------------------------------')
# print(factors_list_of_equation0, ' = 0')
# print('--------------------------------------------------------------------------')
Equation0 = factors_list_of_equation0[-1][3][0]
u0_solution = solve(Equation0, u0)[0]
Equation0_residual = Equation0.subs(u0, u0_solution)
print('--------------------------------------------------------------------------')
print('u0 = ', showDiffNotation(u0_solution))
print('--------------------------------------------------------------------------')
print('Equation0_residual = ', Equation0_residual)
print('--------------------------------------------------------------------------')

Equation1 = expand(modifiedKdVPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 3))
Equation1 = expand(Equation1.subs(u0, u0_solution).doit())
print('--------------------------------------------------------------------------')
print('P1 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation1), ' = 0')
print('--------------------------------------------------------------------------')
Equation1_compatibility_test = expand(Equation1/(24*sigma**2*diff(phi,x)**2))
print('--------------------------------------------------------------------------')
print('Equation for 1st power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation1_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation2 = expand(modifiedKdVPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 2))
print('--------------------------------------------------------------------------')
print('P2 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation2), ' = 0')
print('--------------------------------------------------------------------------')
Equation2 = expand(Equation2.subs(u0, u0_solution).doit())
Equation2_compatibility_test = expand(Equation2/(2*sigma))
print('--------------------------------------------------------------------------')
print('Equation for 2nd power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation2_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation3 = expand(modifiedKdVPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 1))
print('--------------------------------------------------------------------------')
print('P3 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation3), ' = 0')
print('--------------------------------------------------------------------------')
Equation3 = expand(Equation3.subs(u0, u0_solution).doit())
Equation3_compatibility_test = expand(Equation3/(2*sigma))
print('--------------------------------------------------------------------------')
print('Equation for 3rd power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation3_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation4_compatibility_test = expand(modifiedKdVPainleveExpansion_polyInTermsOfPhi.coeff(phi, 0))
print('--------------------------------------------------------------------------')
print('P4 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation4_compatibility_test), ' = 0')
print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')
print('Equation for 4th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation4_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

# print(CompatibilitySystemPDEs([
# 	Equation1_compatibility_test,
# 	Equation2_compatibility_test,
# 	Equation3_compatibility_test,
# 	Equation4_compatibility_test
# 	], phi))

# 0 \in g_9

phi_xx_solution = solve(Equation1_compatibility_test, diff(phi,x,x))[0]

print('--------------------------------------------------------------------------')
print('phi_xx_solution =')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_xx_solution))
print('--------------------------------------------------------------------------')

phi_t_solution = solve(Equation2_compatibility_test, diff(phi,t))[0]

print('--------------------------------------------------------------------------')
print('phi_t_solution =')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_t_solution))
print('--------------------------------------------------------------------------')

u1_t_solution = solve(Equation4_compatibility_test, diff(u1,t))[0]

print('--------------------------------------------------------------------------')
print('u1_t_solution =')
print('--------------------------------------------------------------------------')
print(showDiffNotation(u1_t_solution))
print('--------------------------------------------------------------------------')

dictionary_of_solutions = {}
dictionary_of_solutions[diff(phi,x,x)] = phi_xx_solution
dictionary_of_solutions[diff(phi,t)] = phi_t_solution
dictionary_of_solutions[diff(u1,t)] = u1_t_solution

phi_xxx_solution = diff(dictionary_of_solutions[diff(phi,x,x)],x)
phi_xxx_solution = expand(phi_xxx_solution.subs(diff(phi,x,x),dictionary_of_solutions[diff(phi,x,x)]))
print('--------------------------------------------------------------------------')
print('phi_xxx_solution =')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_xxx_solution))
print('--------------------------------------------------------------------------')

dictionary_of_solutions[diff(phi,x,x,x)] = phi_xxx_solution

dictionary_of_solutions[diff(phi,t)] = expand(dictionary_of_solutions[diff(phi,t)].subs(diff(phi,x,x,x),dictionary_of_solutions[diff(phi,x,x,x)]).doit())
dictionary_of_solutions[diff(phi,t)] = expand(dictionary_of_solutions[diff(phi,t)].subs(diff(phi,x,x),dictionary_of_solutions[diff(phi,x,x)]).doit())
print('--------------------------------------------------------------------------')
print('phi_t_solution =')
print('--------------------------------------------------------------------------')
print(showDiffNotation(dictionary_of_solutions[diff(phi,t)]))
print('--------------------------------------------------------------------------')

phi_xxt_solution_1 = diff(dictionary_of_solutions[diff(phi,x,x)],t)
for already_var in dictionary_of_solutions:
	phi_xxt_solution_1 = expand(phi_xxt_solution_1.subs(already_var,dictionary_of_solutions[already_var]).doit())
for already_var in dictionary_of_solutions:
	phi_xxt_solution_1 = expand(phi_xxt_solution_1.subs(already_var,dictionary_of_solutions[already_var]).doit())
print('--------------------------------------------------------------------------')
print('phi_xxt_solution_1 = diff(phi_xx, t) =')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_xxt_solution_1))
print('--------------------------------------------------------------------------')

phi_xt_solution = diff(dictionary_of_solutions[diff(phi,t)],x)
for already_var in dictionary_of_solutions:
	phi_xt_solution = expand(phi_xt_solution.subs(already_var,dictionary_of_solutions[already_var]).doit())
for already_var in dictionary_of_solutions:
	phi_xt_solution = expand(phi_xt_solution.subs(already_var,dictionary_of_solutions[already_var]).doit())
print('--------------------------------------------------------------------------')
print('phi_xt_solution = diff(phi_t, x) =')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_xt_solution))
print('--------------------------------------------------------------------------')

dictionary_of_solutions[diff(phi,x,t)] = phi_xt_solution

phi_xxt_solution_2 = diff(dictionary_of_solutions[diff(phi,x,t)],x)
for already_var in dictionary_of_solutions:
	phi_xxt_solution_2 = expand(phi_xxt_solution_2.subs(already_var,dictionary_of_solutions[already_var]).doit())
print('--------------------------------------------------------------------------')
print('phi_xxt_solution_2 = diff(phi_xt, x) =')
print('--------------------------------------------------------------------------')
print(showDiffNotation(phi_xxt_solution_2))
print('--------------------------------------------------------------------------')

compatibility_test = expand(phi_xxt_solution_1 - phi_xxt_solution_2)
print('--------------------------------------------------------------------------')
print('phi_xxt_solution_1 - phi_xxt_solution_2 =')
print('--------------------------------------------------------------------------')
print(showDiffNotation(compatibility_test))
print('--------------------------------------------------------------------------')

# print('--------------------------------------------------------------------------')
# print('Reconciliation equation work below ')
# print('--------------------------------------------------------------------------')
# exact_partial_differential_from_equation_3 = diff(phi,t) + 2*sigma**2*diff(phi,x,x,x) - 3*diff(phi,x)*u1**2
# print(showDiffNotation(Equation3_compatibility_test-diff(exact_partial_differential_from_equation_3,x)))
# print('--------------------------------------------------------------------------')

# print('--------------------------------------------------------------------------')
# print('\\frac\\{ \\partial \\}\\{ \\partial x\\}(exact_partial_differential_from_equation_3) = Equation3_compatibility_test'+
# 	', exact_partial_differential_from_equation_3 = ')
# print('--------------------------------------------------------------------------')
# print(showDiffNotation(exact_partial_differential_from_equation_3))
# print('--------------------------------------------------------------------------')

# reconciliation_test = expand(
# 	diff(phi,x)*exact_partial_differential_from_equation_3+\
# 	Equation2_compatibility_test+\
# 	diff(6*sigma*diff(phi,x)*Equation1_compatibility_test,x)
# 	)
# print('--------------------------------------------------------------------------')
# print('reconciliation_test = ')
# print('--------------------------------------------------------------------------')
# print(showDiffNotation(reconciliation_test))
# print('--------------------------------------------------------------------------')

