# KdV one order higher

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

def CompatibilitySystemPDEs(systemsOfPDEs, func, depth_max = 25, depth = 0):

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

KdVTheoryEquation = parse_expr('diff(f(x,t),t)+f(x,t)*diff(f(x,t),x)+sigma*diff(f(x,t),x,x,x)')

u0 = Function('u0')(x, t)
u1 = Function('u1')(x, t)
u2 = Function('u2')(x, t)
u3 = Function('u3')(x, t)

u = u0/phi**2 + u1/phi + u2 + u3*phi
KdVPainleveExpansion = expand(KdVTheoryEquation.subs(f(x, t), u).doit())

generator_list = [
	phi, 1/phi,
	diff(phi, x), diff(phi, t),
	diff(phi, x, x), diff(phi, x, t), diff(phi, t, t),
	diff(phi, x, x, x), diff(phi, x, x, t), diff(phi, x, t, t), diff(phi, t, t, t),
	diff(phi, x, x, x, x), diff(phi, x, x, x, t), diff(phi, x, x, t, t), diff(phi, x, t, t, t), diff(phi, t, t, t, t),
	diff(phi, x, x, x, x, x), diff(phi, x, x, x, x, t), diff(phi, x, x, x, t, t), diff(phi, x, x, t, t, t), diff(phi, x, t, t, t, t), diff(phi, t, t, t, t, t)
	]

KdVPainleveExpansion_polyInTermsOfPhi = Poly(KdVPainleveExpansion, gens = generator_list)
KdVPainleveExpansion_polyInTermsOfPhi = KdVPainleveExpansion_polyInTermsOfPhi.as_expr()

Equation0 = simplify(KdVPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 5))
factors_list_of_equation0 = factor_list(Equation0)
Equation0 = factors_list_of_equation0[-1][-1][0]
u0_solution = solve(Equation0, u0)[0]
Equation0_residual = Equation0.subs(u0, u0_solution)
print('--------------------------------------------------------------------------')
print('Equation for 0th power coefficient')
print('--------------------------------------------------------------------------')
print('u0 = ', u0_solution)
print('--------------------------------------------------------------------------')
print('Equation0_residual = ', Equation0_residual)
print('--------------------------------------------------------------------------')
# print('--------------------------------------------------------------------------')
# print('Equation for 0th power coefficient')
# print('--------------------------------------------------------------------------')
# print(showDiffNotation(Equation0),' = 0')
# print('--------------------------------------------------------------------------')

Equation1 = expand(KdVPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 4))
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
# print('Equation for 1st power coefficient')
# print('--------------------------------------------------------------------------')
# print(showDiffNotation(Equation1Storage),' = 0')
# print('--------------------------------------------------------------------------')

Equation2 = expand(KdVPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 3))
Equation2 = expand(Equation2.subs(u0, u0_solution).doit())
Equation2 = expand(Equation2.subs(u1, u1_solution).doit())
Equation2, _ = fraction(together(Equation2).doit())
Equation2 = expand(Equation2)
Equation2Storage = Equation2
factors_list_of_equation2 = factor_list(Equation2)
# print('--------------------------------------------------------------------------')
# print(factors_list_of_equation2[-1][0][0])
# print('--------------------------------------------------------------------------')
Equation2 = factors_list_of_equation2[-1][-1][0]
u2_solution = solve(Equation2, u2)[0]
# Equation2_residual = Equation2.subs(u2, u2_solution)
print('--------------------------------------------------------------------------')
print('u2 = ', u2_solution)
print('--------------------------------------------------------------------------')
print('Equation for 2nd power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(expand(Equation2Storage/(24*sigma*diff(phi, x)))),' = 0')
print('--------------------------------------------------------------------------')

Equation3 = expand(KdVPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 2))
Equation3 = expand(Equation3.subs(u0, u0_solution).doit())
Equation3 = expand(Equation3.subs(u1, u1_solution).doit())
Equation3 = expand(Equation3.subs(u2, u2_solution).doit())
Equation3, _ = fraction(together(Equation3).doit())
Equation3 = expand(Equation3)
Equation3_compatibility_test = Equation3
# factors_list_of_equation3 = factor_list(Equation3)
# Equation3 = factors_list_of_equation3[-1][-1][0]
print('--------------------------------------------------------------------------')
print('Equation for 3rd power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation3_compatibility_test), '= 0')
print('--------------------------------------------------------------------------')

Equation4 = expand(KdVPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 1))
Equation4 = expand(Equation4.subs(u0, u0_solution).doit())
Equation4 = expand(Equation4.subs(u1, u1_solution).doit())
# Equation4 = expand(Equation3.subs(u2, u2_solution).doit())
Equation4, _ = fraction(together(Equation4).doit())
Equation4 = expand(Equation4)
Equation4_compatibility_test = Equation4
# factors_list_of_equation3 = factor_list(Equation3)
# Equation3 = factors_list_of_equation3[-1][-1][0]
print('--------------------------------------------------------------------------')
print('Equation for 4th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation4_compatibility_test)/12, '= 0')
print('--------------------------------------------------------------------------')

Equation5 = expand(KdVPainleveExpansion_polyInTermsOfPhi.coeff(phi, 0))
Equation5 = expand(Equation5.subs(u0, u0_solution).doit())
Equation5 = expand(Equation5.subs(u1, u1_solution).doit())
# Equation5 = expand(Equation5.subs(u2, u2_solution).doit())
Equation5, _ = fraction(together(Equation5).doit())
Equation5 = expand(Equation5)
Equation5_compatibility_test = Equation5
# factors_list_of_equation3 = factor_list(Equation3)
# Equation3 = factors_list_of_equation3[-1][-1][0]
print('--------------------------------------------------------------------------')
print('Equation for 5th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation5_compatibility_test), '= 0')
print('--------------------------------------------------------------------------')

Equation6 = expand(KdVPainleveExpansion_polyInTermsOfPhi.coeff(phi, 1))
Equation6 = expand(Equation6.subs(u0, u0_solution).doit())
Equation6 = expand(Equation6.subs(u1, u1_solution).doit())
# Equation6 = expand(Equation6.subs(u2, u2_solution).doit())
Equation6, _ = fraction(together(Equation6).doit())
Equation6 = expand(Equation6)
Equation6_compatibility_test = Equation6
# factors_list_of_equation3 = factor_list(Equation3)
# Equation3 = factors_list_of_equation3[-1][-1][0]
print('--------------------------------------------------------------------------')
print('Equation for 6th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation6_compatibility_test), '= 0')
print('--------------------------------------------------------------------------')

Equation7 = expand(KdVPainleveExpansion_polyInTermsOfPhi.coeff(phi, 2))
Equation7 = expand(Equation7.subs(u0, u0_solution).doit())
Equation7 = expand(Equation7.subs(u1, u1_solution).doit())
# Equation7 = expand(Equation7.subs(u2, u2_solution).doit())
Equation7, _ = fraction(together(Equation7).doit())
Equation7 = expand(Equation7)
Equation7_compatibility_test = Equation7
# factors_list_of_equation3 = factor_list(Equation3)
# Equation3 = factors_list_of_equation3[-1][-1][0]
print('--------------------------------------------------------------------------')
print('Equation for 7th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation7_compatibility_test), '= 0')
print('--------------------------------------------------------------------------')


# print(
# 	CompatibilitySystemPDEs([
# 		Equation5_compatibility_test,
# 		Equation6_compatibility_test
# 		], phi)
# )

solution_to_phi_t = solve(Equation2Storage, diff(phi, t))[0]
print('--------------------------------------------------------------------------')
print('phi_t = ', solution_to_phi_t)
print('--------------------------------------------------------------------------')
solution_to_phi_xt_first = expand(diff(solution_to_phi_t, x))
print('--------------------------------------------------------------------------')
print('phi_xt first solution = ', solution_to_phi_xt_first)
print('--------------------------------------------------------------------------')
Equation3_compatibility_test_transformed = expand(Equation3_compatibility_test.subs(diff(phi, t), solution_to_phi_t))
solution_dictionary = solve(
		 [diff(phi, x, t) - solution_to_phi_xt_first, diff(phi, (x, 4)) - Equation3_compatibility_test_transformed],
		 [diff(phi, x, t), diff(phi, (x, 4))]
		 )
solution_to_phi_xt = expand(solution_dictionary[diff(phi, x, t)])
solution_to_phi_x4 = expand(solution_dictionary[diff(phi, (x, 4))])
solution_to_phi_x5 = expand(diff(solution_to_phi_x4, x).doit())
solution_to_phi_x6 = expand(diff(solution_to_phi_x5, x).doit())
solution_dictionary = solve(
		 [
		 	diff(phi, x, t) - solution_to_phi_xt_first,
		 	diff(phi, (x, 4)) - solution_to_phi_x4,
		 	diff(phi, (x, 5)) - solution_to_phi_x5,
		 	diff(phi, (x, 6)) - solution_to_phi_x6
		 ],
		 [diff(phi, x, t), diff(phi, (x, 4)), diff(phi, (x, 5)), diff(phi, (x, 6))]
		 )
solution_to_phi_xt = expand(solution_dictionary[diff(phi, x, t)])
solution_to_phi_x4 = expand(solution_dictionary[diff(phi, (x, 4))])
solution_to_phi_x5 = expand(solution_dictionary[diff(phi, (x, 5))])
solution_to_phi_x6 = expand(solution_dictionary[diff(phi, (x, 6))])
print('--------------------------------------------------------------------------')
print('phi_xt = ', solution_to_phi_xt)
print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')
print('phi_\\{xxxx\\} = ', solution_to_phi_x4)
print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')
print('phi_\\{xxxxx\\} = ', solution_to_phi_x5)
print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')
print('phi_\\{xxxxxx\\} = ', solution_to_phi_x6)
print('--------------------------------------------------------------------------')
solution_to_u2t = expand(solve(Equation5_compatibility_test, diff(u2, t))[0])
print('--------------------------------------------------------------------------')
print('u2_t = ', solution_to_u2t)
print('--------------------------------------------------------------------------')
solution_to_u3t = expand(solve(Equation6_compatibility_test, diff(u3, t))[0])
print('--------------------------------------------------------------------------')
print('u3_t = ', solution_to_u3t)
print('--------------------------------------------------------------------------')
solution_to_phi_xxt = expand(diff(solution_to_phi_xt,x).doit())
solution_to_phi_xxt = expand(solution_to_phi_xxt.subs(diff(phi,(x,6)), solution_to_phi_x6).doit())
solution_to_phi_xxt = expand(solution_to_phi_xxt.subs(diff(phi,(x,5)), solution_to_phi_x5).doit())
solution_to_phi_xxt = expand(solution_to_phi_xxt.subs(diff(phi,(x,4)), solution_to_phi_x4).doit())
print('--------------------------------------------------------------------------')
print('phi_xxt = ', solution_to_phi_xxt)
print('--------------------------------------------------------------------------')
solution_to_phi_xxxt = expand(diff(solution_to_phi_xxt,x).doit())
solution_to_phi_xxxt = expand(solution_to_phi_xxxt.subs(diff(phi,(x,6)), solution_to_phi_x6).doit())
solution_to_phi_xxxt = expand(solution_to_phi_xxxt.subs(diff(phi,(x,5)), solution_to_phi_x5).doit())
solution_to_phi_xxxt = expand(solution_to_phi_xxxt.subs(diff(phi,(x,4)), solution_to_phi_x4).doit())
print('--------------------------------------------------------------------------')
print('phi_xxxt = ', solution_to_phi_xxxt)
print('--------------------------------------------------------------------------')
solution_to_phi_xxxxt_1 = expand(diff(solution_to_phi_xxxt,x).doit())
solution_to_phi_xxxxt_1 = expand(solution_to_phi_xxxxt_1.subs(diff(phi,(x,6)), solution_to_phi_x6).doit())
solution_to_phi_xxxxt_1 = expand(solution_to_phi_xxxxt_1.subs(diff(phi,(x,5)), solution_to_phi_x5).doit())
solution_to_phi_xxxxt_1 = expand(solution_to_phi_xxxxt_1.subs(diff(phi,(x,4)), solution_to_phi_x4).doit())
print('--------------------------------------------------------------------------')
print('phi_xxxt from route 1 = ', solution_to_phi_xxxxt_1)
print('--------------------------------------------------------------------------')
solution_to_phi_xxxxt_2 = expand(diff(solution_to_phi_x4,t).doit())
solution_to_phi_xxxxt_2 = expand(solution_to_phi_xxxxt_2.subs(diff(u2, t), solution_to_u2t).doit())
solution_to_phi_xxxxt_2 = expand(solution_to_phi_xxxxt_2.subs(diff(u3, t), solution_to_u3t).doit())
solution_to_phi_xxxxt_2 = expand(solution_to_phi_xxxxt_2.subs(diff(phi,(x,6)), solution_to_phi_x6).doit())
solution_to_phi_xxxxt_2 = expand(solution_to_phi_xxxxt_2.subs(diff(phi,(x,5)), solution_to_phi_x5).doit())
solution_to_phi_xxxxt_2 = expand(solution_to_phi_xxxxt_2.subs(diff(phi,(x,4)), solution_to_phi_x4).doit())
print('--------------------------------------------------------------------------')
print('phi_xxxt from route 2 = ', solution_to_phi_xxxxt_2)
print('--------------------------------------------------------------------------')
compatibility_Test = expand(solution_to_phi_xxxxt_2-solution_to_phi_xxxxt_1)
# compatibility_Test, _ = fraction(together(compatibility_Test).doit())
# compatibility_Test = expand(compatibility_Test)
print('--------------------------------------------------------------------------')
print('compatibility_Test = ', compatibility_Test)
print('--------------------------------------------------------------------------')
