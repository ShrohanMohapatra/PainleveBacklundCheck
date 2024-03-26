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
u2 = Function('u2')(x, t)

u = u0/phi + u1 + u2*phi
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

# print('--------------------------------------------------------------------------')
# print('modifiedKdVPainleveExpansion_polyInTermsOfPhi =')
# print('--------------------------------------------------------------------------')
# print(showDiffNotation(modifiedKdVPainleveExpansion_polyInTermsOfPhi))
# print('--------------------------------------------------------------------------')

Equation0 = simplify(modifiedKdVPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 4))
print('--------------------------------------------------------------------------')
print('P0 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation0), ' = 0')
print('--------------------------------------------------------------------------')
factors_list_of_equation0 = factor_list(Equation0)
Equation0 = factors_list_of_equation0[-1][2][0]
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
factors_list_of_equation1 = factor_list(Equation1)
Equation1 = factors_list_of_equation1[-1][-1][0]
u1_solution = solve(Equation1, u1)[0]
Equation1_residual = Equation1.subs(u1, u1_solution)
print('--------------------------------------------------------------------------')
print('u1 = ', showDiffNotation(u1_solution))
print('--------------------------------------------------------------------------')
print('Equation1_residual = ', Equation1_residual)
print('--------------------------------------------------------------------------')

Equation2 = expand(modifiedKdVPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 2))
print('--------------------------------------------------------------------------')
print('P2 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation2), ' = 0')
print('--------------------------------------------------------------------------')
Equation2 = expand(Equation2.subs(u0, u0_solution).doit())
Equation2 = expand(Equation2.subs(u1, u1_solution).doit())
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
Equation3 = expand(Equation3.subs(u1, u1_solution).doit())
Equation3_compatibility_test = expand(Equation3/(2*sigma))
Equation3_compatibility_test, _ = fraction(together(Equation3_compatibility_test).doit())
Equation3_compatibility_test = expand(Equation3_compatibility_test)
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
Equation4_compatibility_test = expand(Equation4_compatibility_test.subs(u0, u0_solution).doit())
Equation4_compatibility_test = expand(Equation4_compatibility_test.subs(u1, u1_solution).doit())
Equation4_compatibility_test = expand(Equation4_compatibility_test/(2*sigma))
Equation4_compatibility_test, _ = fraction(together(Equation4_compatibility_test).doit())
Equation4_compatibility_test = expand(Equation4_compatibility_test)
print('--------------------------------------------------------------------------')
print('Equation for 4th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation4_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation5_compatibility_test = expand(modifiedKdVPainleveExpansion_polyInTermsOfPhi.coeff(phi, 1))
print('--------------------------------------------------------------------------')
print('P5 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation5_compatibility_test), ' = 0')
print('--------------------------------------------------------------------------')
Equation5_compatibility_test = expand(Equation5_compatibility_test.subs(u0, u0_solution).doit())
Equation5_compatibility_test = expand(Equation5_compatibility_test.subs(u1, u1_solution).doit())
Equation5_compatibility_test = expand(Equation5_compatibility_test/(2*sigma))
Equation5_compatibility_test, _ = fraction(together(Equation5_compatibility_test).doit())
Equation5_compatibility_test = expand(Equation5_compatibility_test)
print('--------------------------------------------------------------------------')
print('Equation for 5th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation5_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation6_compatibility_test = expand(modifiedKdVPainleveExpansion_polyInTermsOfPhi.coeff(phi, 2))
print('--------------------------------------------------------------------------')
print('P6 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation6_compatibility_test), ' = 0')
print('--------------------------------------------------------------------------')
Equation6_compatibility_test = expand(Equation6_compatibility_test.subs(u0, u0_solution).doit())
Equation6_compatibility_test = expand(Equation6_compatibility_test.subs(u1, u1_solution).doit())
Equation6_compatibility_test = expand(Equation6_compatibility_test/(2*sigma))
Equation6_compatibility_test, _ = fraction(together(Equation6_compatibility_test).doit())
Equation6_compatibility_test = expand(Equation6_compatibility_test)
print('--------------------------------------------------------------------------')
print('Equation for 6th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation6_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation7_compatibility_test = expand(modifiedKdVPainleveExpansion_polyInTermsOfPhi.coeff(phi, 3))
print('--------------------------------------------------------------------------')
print('P7 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation7_compatibility_test), ' = 0')
print('--------------------------------------------------------------------------')
Equation7_compatibility_test = expand(Equation7_compatibility_test.subs(u0, u0_solution).doit())
Equation7_compatibility_test = expand(Equation7_compatibility_test.subs(u1, u1_solution).doit())
Equation7_compatibility_test = expand(Equation7_compatibility_test/(2*sigma))
Equation7_compatibility_test, _ = fraction(together(Equation7_compatibility_test).doit())
Equation7_compatibility_test = expand(Equation7_compatibility_test)
print('--------------------------------------------------------------------------')
print('Equation for 7th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation7_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation3_compatibility_test = expand(Equation3_compatibility_test.subs(diff(u2, x), 0))
print('--------------------------------------------------------------------------')
print('Equation3_compatibility_test = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation3_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation4_compatibility_test = expand(Equation4_compatibility_test.subs(diff(u2, x), 0))
print('--------------------------------------------------------------------------')
print('Equation4_compatibility_test = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation4_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation5_compatibility_test = expand(Equation5_compatibility_test.subs(diff(u2, x), 0).doit())
print('--------------------------------------------------------------------------')
print('Equation5_compatibility_test = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation5_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

Equation6_compatibility_test = expand(Equation6_compatibility_test.subs(diff(u2, x), 0).doit())
print('--------------------------------------------------------------------------')
print('Equation6_compatibility_test = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation6_compatibility_test),' = 0')
print('--------------------------------------------------------------------------')

# print(CompatibilitySystemPDEs([
# 	# Equation1,
# 	# Equation2_compatibility_test,
# 	Equation3_compatibility_test,
# 	Equation4_compatibility_test,
# 	Equation5_compatibility_test,
# 	Equation6_compatibility_test#,
# 	# Equation7_compatibility_test
# 	], phi))

P3_P4 = expand(diff(Equation3_compatibility_test*diff(phi,x),x)+Equation4_compatibility_test)
print('--------------------------------------------------------------------------')
print('P3_P4 = diff(Equation3_compatibility_test*diff(phi,x),x)+Equation4_compatibility_test = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(P3_P4),' = 0')
print('--------------------------------------------------------------------------')

P3_P5 = expand(diff(Equation3_compatibility_test,t)+12*sigma*diff(phi,x,x)*Equation5_compatibility_test)
print('--------------------------------------------------------------------------')
print('P3_P5 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(P3_P5),' = 0')
print('--------------------------------------------------------------------------')

P6 = expand(Equation6_compatibility_test/(3*u2**2))
print('--------------------------------------------------------------------------')
print('P6 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(P6),' = 0')
print('--------------------------------------------------------------------------')

print('--------------------------------------------------------------------------')
print('P3_P4 = diff(Equation3_compatibility_test*diff(phi,x),x)+Equation4_compatibility_test = ')
print('--------------------------------------------------------------------------')
print(P3_P4)
print('--------------------------------------------------------------------------')
print(showDiffNotation(P3_P4),' = 0')
print('--------------------------------------------------------------------------')

P3_P5_P6 = expand(P3_P5 + 6*sigma*diff(phi,t,x,x)*P6)
print('--------------------------------------------------------------------------')
print('P3_P5_P6 = ')
print('--------------------------------------------------------------------------')
print(P3_P5_P6)
print('--------------------------------------------------------------------------')
print(showDiffNotation(P3_P5_P6),' = 0')
print('--------------------------------------------------------------------------')

reconciliation1 = expand(
	P3_P4*(-72*sigma**3*u2*diff(phi, x)*diff(phi, (x, 2))**2*diff(phi, (x, 3)) + 72*sigma**3*u2*diff(phi, (x, 2))**4 + 2*sigma**2*diff(phi, x)**2*diff(phi, t, (x, 4)) - 6*sigma**2*diff(phi, x)*diff(phi, (x, 2))*diff(phi, t, (x, 3)) + 4*sigma**2*diff(phi, x)*diff(phi, (x, 4))*diff(phi, t, x) + 3*sigma**2*diff(phi, (x, 2))**2*diff(phi, t, (x, 2)) - 6*sigma**2*diff(phi, (x, 2))*diff(phi, (x, 3))*diff(phi, t, x) - 18*sigma*u2*diff(phi, x)**3*diff(phi, t, (x, 2)) - 36*sigma*u2*diff(phi, x)**2*diff(phi, (x, 2))*diff(phi, t, x) + diff(phi, x)**2*diff(phi, (t, 2), x) + 2*diff(phi, x)*diff(phi, t, x)**2) - \
	P3_P5_P6*(8*sigma**2*diff(phi, x)**2*diff(phi, (x, 2))*diff(phi, (x, 4)) - 24*sigma**2*diff(phi, x)*diff(phi, (x, 2))**2*diff(phi, (x, 3)) + 12*sigma**2*diff(phi, (x, 2))**4 - 48*sigma*u2*diff(phi, x)**3*diff(phi, (x, 2))**2 + 4*diff(phi, x)**2*diff(phi, (x, 2))*diff(phi, t, x))
	)
print('--------------------------------------------------------------------------')
print('reconciliation1 =')
print('--------------------------------------------------------------------------')
print(showDiffNotation(reconciliation1),' = 0')
print('--------------------------------------------------------------------------')

P4_P5 = expand(diff(Equation4_compatibility_test,t)-12*sigma*diff(phi,x,x,x)*diff(phi,x)*Equation5_compatibility_test)
print('--------------------------------------------------------------------------')
print('P4_P5 = ')
print('--------------------------------------------------------------------------')
print(P4_P5)
print('--------------------------------------------------------------------------')
print(showDiffNotation(P4_P5),' = 0')
print('--------------------------------------------------------------------------')

P4_P5_P6 = expand(
	P4_P5*(diff(phi,x,x,x)*diff(phi,x)*sigma - diff(phi,x,x)**2*sigma) - \
	P6*( - 2*sigma**2*diff(phi, x)**3*diff(phi, t, (x, 5)) + 8*sigma**2*diff(phi, x)**2*diff(phi, (x, 2))*diff(phi, t, (x, 4)) + 12*sigma**2*diff(phi, x)**2*diff(phi, (x, 3))*diff(phi, t, (x, 3)) + 8*sigma**2*diff(phi, x)**2*diff(phi, (x, 4))*diff(phi, t, (x, 2)) - 6*sigma**2*diff(phi, x)**2*diff(phi, (x, 5))*diff(phi, t, x) - 21*sigma**2*diff(phi, x)*diff(phi, (x, 2))**2*diff(phi, t, (x, 3)) - 42*sigma**2*diff(phi, x)*diff(phi, (x, 2))*diff(phi, (x, 3))*diff(phi, t, (x, 2)) + 16*sigma**2*diff(phi, x)*diff(phi, (x, 2))*diff(phi, (x, 4))*diff(phi, t, x) + 12*sigma**2*diff(phi, x)*diff(phi, (x, 3))**2*diff(phi, t, x) + 36*sigma**2*diff(phi, (x, 2))**3*diff(phi, t, (x, 2)) - 21*sigma**2*diff(phi, (x, 2))**2*diff(phi, (x, 3))*diff(phi, t, x) - diff(phi, x)**3*diff(phi, (t, 2), (x, 2)) + diff(phi, x)**2*diff(phi, (x, 2))*diff(phi, (t, 2), x) - 2*diff(phi, x)**2*diff(phi, t, x)*diff(phi, t, (x, 2)) + 2*diff(phi, x)*diff(phi, (x, 2))*diff(phi, t, x)**2)
	)
print('--------------------------------------------------------------------------')
print('P4_P5_P6 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(P4_P5_P6),' = 0')
print('--------------------------------------------------------------------------')

divisor_to_P4_P5_P6 = expand(P4_P5_P6/u2)
print('--------------------------------------------------------------------------')
print('divisor_to_P4_P5_P6 = expand(P4_P5_P6/u2) = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(divisor_to_P4_P5_P6))
print('--------------------------------------------------------------------------')

print('--------------------------------------------------------------------------')
print('reconciliation2 is simply diff(P4_P5_P6/divisor_to_P4_P5_P6,x) = 0')
print('--------------------------------------------------------------------------')

f3mKdV_higher = Function('f3mKdV_higher')(x, t)
f4mKdV_higher = Function('f4mKdV_higher')(x, t)
f5mKdV_higher = Function('f5mKdV_higher')(x, t)
f6mKdV_higher = Function('f6mKdV_higher')(x, t)

P3_P4 = expand(diff(f3mKdV_higher*diff(phi,x),x)+f4mKdV_higher)
P3_P5 = expand(diff(f3mKdV_higher,t)+12*sigma*diff(phi,x,x)*f5mKdV_higher)
P3_P5_P6 = expand(P3_P5 + 6*sigma*diff(phi,t,x,x)*f6mKdV_higher/(3*u2**2))

reconciliation1_sub = expand(
	P3_P4*(-72*sigma**3*u2*diff(phi, x)*diff(phi, (x, 2))**2*diff(phi, (x, 3)) + 72*sigma**3*u2*diff(phi, (x, 2))**4 + 2*sigma**2*diff(phi, x)**2*diff(phi, t, (x, 4)) - 6*sigma**2*diff(phi, x)*diff(phi, (x, 2))*diff(phi, t, (x, 3)) + 4*sigma**2*diff(phi, x)*diff(phi, (x, 4))*diff(phi, t, x) + 3*sigma**2*diff(phi, (x, 2))**2*diff(phi, t, (x, 2)) - 6*sigma**2*diff(phi, (x, 2))*diff(phi, (x, 3))*diff(phi, t, x) - 18*sigma*u2*diff(phi, x)**3*diff(phi, t, (x, 2)) - 36*sigma*u2*diff(phi, x)**2*diff(phi, (x, 2))*diff(phi, t, x) + diff(phi, x)**2*diff(phi, (t, 2), x) + 2*diff(phi, x)*diff(phi, t, x)**2) - \
	P3_P5_P6*(8*sigma**2*diff(phi, x)**2*diff(phi, (x, 2))*diff(phi, (x, 4)) - 24*sigma**2*diff(phi, x)*diff(phi, (x, 2))**2*diff(phi, (x, 3)) + 12*sigma**2*diff(phi, (x, 2))**4 - 48*sigma*u2*diff(phi, x)**3*diff(phi, (x, 2))**2 + 4*diff(phi, x)**2*diff(phi, (x, 2))*diff(phi, t, x))
	)
print('--------------------------------------------------------------------------')
print('reconciliation1_sub =')
print('--------------------------------------------------------------------------')
print(showDiffNotation(reconciliation1_sub),' = ',showDiffNotation(reconciliation1))
print('--------------------------------------------------------------------------')

P4_P5 = expand(diff(f4mKdV_higher,t)-12*sigma*diff(phi,x,x,x)*diff(phi,x)*f5mKdV_higher)
P4_P5_P6 = expand(
	P4_P5*(diff(phi,x,x,x)*diff(phi,x)*sigma - diff(phi,x,x)**2*sigma) - \
	f6mKdV_higher/(3*u2**2)*( - 2*sigma**2*diff(phi, x)**3*diff(phi, t, (x, 5)) + 8*sigma**2*diff(phi, x)**2*diff(phi, (x, 2))*diff(phi, t, (x, 4)) + 12*sigma**2*diff(phi, x)**2*diff(phi, (x, 3))*diff(phi, t, (x, 3)) + 8*sigma**2*diff(phi, x)**2*diff(phi, (x, 4))*diff(phi, t, (x, 2)) - 6*sigma**2*diff(phi, x)**2*diff(phi, (x, 5))*diff(phi, t, x) - 21*sigma**2*diff(phi, x)*diff(phi, (x, 2))**2*diff(phi, t, (x, 3)) - 42*sigma**2*diff(phi, x)*diff(phi, (x, 2))*diff(phi, (x, 3))*diff(phi, t, (x, 2)) + 16*sigma**2*diff(phi, x)*diff(phi, (x, 2))*diff(phi, (x, 4))*diff(phi, t, x) + 12*sigma**2*diff(phi, x)*diff(phi, (x, 3))**2*diff(phi, t, x) + 36*sigma**2*diff(phi, (x, 2))**3*diff(phi, t, (x, 2)) - 21*sigma**2*diff(phi, (x, 2))**2*diff(phi, (x, 3))*diff(phi, t, x) - diff(phi, x)**3*diff(phi, (t, 2), (x, 2)) + diff(phi, x)**2*diff(phi, (x, 2))*diff(phi, (t, 2), x) - 2*diff(phi, x)**2*diff(phi, t, x)*diff(phi, t, (x, 2)) + 2*diff(phi, x)*diff(phi, (x, 2))*diff(phi, t, x)**2)
	)

print('--------------------------------------------------------------------------')
print('u2**2*diff(P4_P5_P6/divisor_to_P4_P5_P6,x) = f7mKdV_Higher, where P4_P5_P6 ->')
print('--------------------------------------------------------------------------')
print(showDiffNotation(P4_P5_P6))
print('--------------------------------------------------------------------------')
print('and divisor_to_P4_P5_P6 =')
print('--------------------------------------------------------------------------')
print(showDiffNotation(divisor_to_P4_P5_P6))
print('--------------------------------------------------------------------------')


