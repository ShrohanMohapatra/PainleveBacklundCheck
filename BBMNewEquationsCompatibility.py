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

x, t, sigma, b, alpha, beta = symbols('x t sigma b alpha beta')
phi = Function('phi')(x, t)
f = Function('f')

BBMEquation = parse_expr('diff(f(x,t),t)+diff(f(x,t),x)+f(x,t)*diff(f(x,t),x)-diff(f(x,t),x,x,t)')

u0 = Function('u0')(x, t)
u1 = Function('u1')(x, t)
u2 = Function('u2')(x, t)

u = u0/phi**2 + u1/phi + u2
BBMPainleveExpansion = expand(BBMEquation.subs(f(x, t), u).doit())

generator_list = [
	phi, 1/phi,
	diff(phi, x), diff(phi, t),
	diff(phi, x, x), diff(phi, x, t), diff(phi, t, t),
	diff(phi, x, x, x), diff(phi, x, x, t), diff(phi, x, t, t), diff(phi, t, t, t),
	diff(phi, x, x, x, x), diff(phi, x, x, x, t), diff(phi, x, x, t, t), diff(phi, x, t, t, t), diff(phi, t, t, t, t),
	diff(phi, x, x, x, x, x), diff(phi, x, x, x, x, t), diff(phi, x, x, x, t, t), diff(phi, x, x, t, t, t), diff(phi, x, t, t, t, t), diff(phi, t, t, t, t, t)
	]

BBMPainleveExpansion_polyInTermsOfPhi = Poly(BBMPainleveExpansion, gens = generator_list)
BBMPainleveExpansion_polyInTermsOfPhi = BBMPainleveExpansion_polyInTermsOfPhi.as_expr()

Equation0 = simplify(BBMPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 5))
factors_list_of_equation0 = factor_list(Equation0)
Equation0 = factors_list_of_equation0[-1][-1][0]
u0_solution = solve(Equation0, u0)[0]
Equation0_residual = Equation0.subs(u0, u0_solution)
print('--------------------------------------------------------------------------')
print('u0 = ', u0_solution)
print('--------------------------------------------------------------------------')
print('Equation0_residual = ', Equation0_residual)
print('--------------------------------------------------------------------------')

Equation1 = expand(BBMPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 4))
Equation1 = expand(Equation1.subs(u0, u0_solution).doit())
factors_list_of_equation1 = factor_list(Equation1)
Equation1 = factors_list_of_equation1[-1][-1][0]
# u1_solution = solve(Equation1, u1)[0]
# u1_solution = together(u1_solution).doit()
# Equation1_residual = Equation1.subs(u1, u1_solution)
# Equation1_residual, _ = fraction(together(Equation1_residual).doit())
# Equation1_residual = expand(Equation1_residual)
print('--------------------------------------------------------------------------')
print('Equation for 1st power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation1), '= 0')
print('--------------------------------------------------------------------------')
# print('Equation1_residual = ', Equation1_residual)
# print('--------------------------------------------------------------------------')

Equation2 = expand(BBMPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 3))
Equation2 = expand(Equation2.subs(u0, u0_solution).doit())
# Equation2 = expand(Equation2.subs(u1, u1_solution).doit())
Equation2, _ = fraction(together(Equation2).doit())
Equation2 = expand(Equation2)
factors_list_of_equation2 = factor_list(Equation2)
# print('--------------------------------------------------------------------------')
# print(factors_list_of_equation2[-1][0][0])
# print('--------------------------------------------------------------------------')
Equation2 = factors_list_of_equation2[-1][-1][0]
# u2_solution = solve(Equation2, u2)[0]
# Equation2_residual = Equation2.subs(u2, u2_solution)
# print('--------------------------------------------------------------------------')
# print('u2 = ', u2_solution)
print('--------------------------------------------------------------------------')
print('Equation for 2nd power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation2),'= 0')
print('--------------------------------------------------------------------------')

Equation3 = expand(BBMPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 2))
Equation3 = expand(Equation3.subs(u0, u0_solution).doit())
# Equation3 = expand(Equation3.subs(u1, u1_solution).doit())
# Equation3_compatibility_test = Equation3
# Equation3 = expand(Equation3.subs(u2, u2_solution).doit())
Equation3, _ = fraction(together(Equation3).doit())
Equation3 = expand(Equation3)
factors_list_of_equation3 = factor_list(Equation3)
Equation3 = factors_list_of_equation3[-1][-1][0]
print('--------------------------------------------------------------------------')
print('Equation for 3rd power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation3), '= 0')
print('--------------------------------------------------------------------------')

Equation4 = expand(BBMPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 1))
# Equation4 = expand(Equation4.subs(u0, u0_solution).doit())
# Equation4 = expand(Equation4.subs(u1, u1_solution).doit())
# Equation4_compatibility_test = Equation4
# Equation4 = expand(Equation4.subs(u2, u2_solution).doit())
Equation4, _ = fraction(together(Equation4).doit())
Equation4 = expand(Equation4)
factors_list_of_equation4 = factor_list(Equation4)
Equation4 = factors_list_of_equation4[-1][-1][0]
print('--------------------------------------------------------------------------')
print('Equation for 4th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation4),' = 0')
print('--------------------------------------------------------------------------')

Equation5 = expand(BBMPainleveExpansion_polyInTermsOfPhi.coeff(phi, 0))
# Equation5 = expand(Equation5.subs(u0, u0_solution).doit())
# Equation5 = expand(Equation5.subs(u1, u1_solution).doit())
# Equation5_compatibility_test = Equation5
# Equation5 = expand(Equation5.subs(u2, u2_solution).doit())
Equation5, _ = fraction(together(Equation5).doit())
Equation5 = expand(Equation5)
# factors_list_of_equation5 = factor_list(Equation5)
# print('--------------------------------------------------------------------------')
# print(len(factors_list_of_equation5))
# print('--------------------------------------------------------------------------')
# Equation5 = factors_list_of_equation5[-1][-1][0]
print('--------------------------------------------------------------------------')
print('Equation for 5th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation5),' = 0')
print('--------------------------------------------------------------------------')

eqn51 = Equation1
eqn54 = Equation3
print('--------------------------------------------------------------------------')
print('Equation 51 from the paper')
print('--------------------------------------------------------------------------')
print(showDiffNotation(eqn51),'= 0')
print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')
print('Equation 54 from the paper')
print('--------------------------------------------------------------------------')
print(showDiffNotation(eqn54),'= 0')
print('--------------------------------------------------------------------------')
solution = solve([eqn51,eqn54],[diff(phi,t,t),diff(phi,t,x,x,x)])

print('--------------------------------------------------------------------------')
print('phi_{t,t} = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(solution[diff(phi,t,t)]))
print('--------------------------------------------------------------------------')

print('--------------------------------------------------------------------------')
print('phi_{t,x,x,x} = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(solution[diff(phi,t,x,x,x)]))
print('--------------------------------------------------------------------------')

eqn55 = Equation4
eqn56 = Equation5
print('--------------------------------------------------------------------------')
print('Equation 55 from the paper')
print('--------------------------------------------------------------------------')
print(showDiffNotation(eqn55),' = 0')
print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')
print('Equation 56 from the paper')
print('--------------------------------------------------------------------------')
print(showDiffNotation(eqn56),' = 0')
print('--------------------------------------------------------------------------')

u1_txx_solution = solve(eqn55,diff(u1,t,x,x))[0] # [diff(u1,t,x,x)]
u2_txx_solution = solve(eqn56,diff(u2,t,x,x))[0] # [diff(u2,t,x,x)]

print('--------------------------------------------------------------------------')
print('u1_{t,x,x} = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(u1_txx_solution))
print('--------------------------------------------------------------------------')

print('--------------------------------------------------------------------------')
print('u2_{t,x,x} = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(u2_txx_solution))
print('--------------------------------------------------------------------------')

dictionary_of_solutions = {} # Initializing the dictionary of solutions
dictionary_of_solutions[diff(phi,t,x,x,x)] = solution[diff(phi,t,x,x,x)]
dictionary_of_solutions[diff(phi,t,t)] = solution[diff(phi,t,t)]
dictionary_of_solutions[diff(u1,t,x,x)] = u1_txx_solution
dictionary_of_solutions[diff(u2,t,x,x)] = u2_txx_solution
# print('--------------------------------------------------------------------------')
# print(dictionary_of_solutions)
# print('--------------------------------------------------------------------------')

max_number_of_derivatives = 3

for m in range(max_number_of_derivatives+1):
	for n in range(max_number_of_derivatives+1):
		print('m = ', m, ' n = ', n,
			diff(u2,(x,(m+2)),(t,(n+1))) )
		if m == 0 and n == 0:
			print('if m == 0 or n == 0 statement entry')
			pass
		elif n == 0:
			print('elif n == 0 statement entry')
			dictionary_of_solutions[diff(u2,(x,(m+2)),t)] = diff(dictionary_of_solutions[diff(u2,(x,(m+1)),t)],x)
			dictionary_of_solutions[diff(u2,(x,(m+2)),t)] = expand(dictionary_of_solutions[diff(u2,(x,(m+2)),t)]).doit()
		elif n >= 1:
			print('elif n >= 1 statement entry')
			dictionary_of_solutions[diff(u2,(x,(m+2)),(t,n+1))] = diff(dictionary_of_solutions[diff(u2,(x,(m+2)),(t,n))],t)
			dictionary_of_solutions[diff(u2,(x,(m+2)),(t,n+1))] = expand(dictionary_of_solutions[diff(u2,(x,(m+2)),(t,n+1))]).doit()
		for old_variable in dictionary_of_solutions:
			if old_variable != diff(u2,(x,m+2),(t,n+1)):
				print('m = ',m,'n = ',n)
				print('new_variable = ',diff(u2,(x,(m+2)),(t,n+1)))
				print('old_variable = ',old_variable)
				dictionary_of_solutions[diff(u2,(x,(m+2)),(t,n+1))] = dictionary_of_solutions[diff(u2,(x,(m+2)),(t,n+1))].subs(
																		old_variable,
																		dictionary_of_solutions[old_variable]
																		)
				dictionary_of_solutions[diff(u2,(x,(m+2)),(t,n+1))] = expand(dictionary_of_solutions[diff(u2,(x,(m+2)),(t,n+1))]).doit()

for m in range(max_number_of_derivatives+1):
	for n in range(max_number_of_derivatives+1):
		print('m = ', m, ' n = ', n,
			diff(u1,(x,(m+2)),(t,(n+1))) )
		if m == 0 and n == 0:
			print('if m == 0 or n == 0 statement entry')
			pass
		elif n == 0:
			print('elif n == 0 statement entry')
			dictionary_of_solutions[diff(u1,(x,(m+2)),t)] = diff(dictionary_of_solutions[diff(u1,(x,(m+1)),t)],x)
			dictionary_of_solutions[diff(u1,(x,(m+2)),t)] = expand(dictionary_of_solutions[diff(u1,(x,(m+2)),t)]).doit()
		elif n >= 1:
			print('elif n >= 1 statement entry')
			dictionary_of_solutions[diff(u1,(x,(m+2)),(t,n+1))] = diff(dictionary_of_solutions[diff(u1,(x,(m+2)),(t,n))],t)
			dictionary_of_solutions[diff(u1,(x,(m+2)),(t,n+1))] = expand(dictionary_of_solutions[diff(u1,(x,(m+2)),(t,n+1))]).doit()
		for old_variable in dictionary_of_solutions:
			if old_variable != diff(u1,(x,(m+2)),(t,n+1)):
				print('m = ',m,'n = ',n)
				print('new_variable = ',diff(u1,(x,(m+2)),(t,n+1)))
				print('old_variable = ',old_variable)
				dictionary_of_solutions[diff(u1,(x,(m+2)),(t,n+1))] = dictionary_of_solutions[diff(u1,(x,(m+2)),(t,n+1))].subs(
																		old_variable,
																		dictionary_of_solutions[old_variable]
																		)
				dictionary_of_solutions[diff(u1,(x,(m+2)),(t,n+1))] = expand(dictionary_of_solutions[diff(u1,(x,(m+2)),(t,n+1))]).doit()

for m in range(1,max_number_of_derivatives+1):
	dictionary_of_solutions[diff(phi,t,(x,m+3))] = diff(dictionary_of_solutions[diff(phi,t,(x,m+2))],x)
	dictionary_of_solutions[diff(phi,t,(x,m+3))] = expand(dictionary_of_solutions[diff(phi,t,(x,m+3))])
	for term in dictionary_of_solutions:
		if term != diff(phi,t,(x,m+3)):
			print('Currently processing the equation for', diff(phi,t,(x,m+3)), 'by the substitution of the term = ',term)
			dictionary_of_solutions[diff(phi,t,(x,m+3))] = dictionary_of_solutions[diff(phi,t,(x,m+3))].subs(
																	term,
																	dictionary_of_solutions[term]
																	)
			dictionary_of_solutions[diff(phi,t,(x,m+3))] = expand(dictionary_of_solutions[diff(phi,t,(x,m+3))]).doit()

equation_list_1 = [expand(diff(phi,t,(x,m+3)) - dictionary_of_solutions[diff(phi,t,(x,m+3))]) for m in range(1,max_number_of_derivatives+1)]
variable_list_1 = [diff(phi,t,(x,m+3)) for m in range(1,max_number_of_derivatives+1)]

solution_list = solve(equation_list_1, variable_list_1)
for m in range(1,max_number_of_derivatives+1):
	dictionary_of_solutions[diff(phi,t,(x,m+3))] = solution_list[diff(phi,t,(x,m+3))]

dictionary_of_solutions[diff(phi,t,t,x)] = diff(solution[diff(phi,t,t)],x)
dictionary_of_solutions[diff(phi,t,t,x)] = expand(dictionary_of_solutions[diff(phi,t,t,x)]).doit()


for term in dictionary_of_solutions:
	if term != diff(phi,t,t,x):
		print('Currently processing the equation for', diff(phi,t,t,x), 'by the substitution of the term = ',term)
		dictionary_of_solutions[diff(phi,t,t,x)] = dictionary_of_solutions[diff(phi,t,t,x)].subs(
																term,
																dictionary_of_solutions[term]
																)
		dictionary_of_solutions[diff(phi,t,t,x)] = expand(dictionary_of_solutions[diff(phi,t,t,x)]).doit()

dictionary_of_solutions[diff(phi,t,t,x,x)] = diff(dictionary_of_solutions[diff(phi,t,t,x)],x)
dictionary_of_solutions[diff(phi,t,t,x,x)] = expand(dictionary_of_solutions[diff(phi,t,t,x,x)]).doit()

for term in dictionary_of_solutions:
	if term != diff(phi,t,t,x,x):
		print('Currently processing the equation for', diff(phi,t,t,x,x), 'by the substitution of the term = ',term)
		dictionary_of_solutions[diff(phi,t,t,x,x)] = dictionary_of_solutions[diff(phi,t,t,x,x)].subs(
																term,
																dictionary_of_solutions[term]
																)
		dictionary_of_solutions[diff(phi,t,t,x,x)] = expand(dictionary_of_solutions[diff(phi,t,t,x,x)]).doit()

equation_list_2 = [expand(diff(phi,t,t,x,x) - dictionary_of_solutions[diff(phi,t,t,x,x)]), expand(diff(phi,t,x,x,x) - dictionary_of_solutions[diff(phi,t,x,x,x)]), ]
variable_list_2 = [diff(phi,t,t,x,x), diff(phi,t,x,x,x)]
solution_list = solve(equation_list_2, variable_list_2)
dictionary_of_solutions[diff(phi,t,t,x,x)] = solution_list[diff(phi,t,t,x,x)]
dictionary_of_solutions[diff(phi,t,x,x,x)] = solution_list[diff(phi,t,x,x,x)]

print('--------------------------------------------------------------------------')
print(diff(phi,t,t,x,x))
print('--------------------------------------------------------------------------')
print(showDiffNotation(solution_list[diff(phi,t,t,x,x)]))
print('--------------------------------------------------------------------------')
print(diff(phi,t,x,x,x))
print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')
print(showDiffNotation(solution_list[diff(phi,t,x,x,x)]))
print('--------------------------------------------------------------------------')

print('--------------------------------------------------------------------------')
for term in dictionary_of_solutions:
	print('--------------------------------------------------------------------------')
	print(showDiffNotation(term), ' = ')
	print('--------------------------------------------------------------------------')
	print(showDiffNotation(dictionary_of_solutions[term]))
	print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')


#########################################################################################################################

equation1 = diff(dictionary_of_solutions[diff(phi,t,t,x,x)],x)
equation1 = expand(equation1).doit()

equation_list_3= [
	expand(diff(phi,t,t,x,x,x) - equation1),
	expand(diff(phi,t,t) - dictionary_of_solutions[diff(phi,t,t)]),
	expand(diff(phi,t,t,x) - dictionary_of_solutions[diff(phi,t,t,x)]),
	expand(diff(phi,t,t,x,x) - dictionary_of_solutions[diff(phi,t,t,x,x)]),
	expand(diff(phi,t,x,x,x) - dictionary_of_solutions[diff(phi,t,x,x,x)])
	]
variable_list_3 = [
	diff(phi,t,t,x,x,x),
	diff(phi,t,t),
	diff(phi,t,t,x),
	diff(phi,t,t,x,x),
	diff(phi,t,x,x,x)
	]

solution_list = solve(equation_list_3, variable_list_3)
equation1 = solution_list[diff(phi,t,t,x,x,x)]
dictionary_of_solutions[diff(phi,t,t)] = solution_list[diff(phi,t,t)]
dictionary_of_solutions[diff(phi,t,t,x)] = solution_list[diff(phi,t,t,x)]
dictionary_of_solutions[diff(phi,t,t,x,x)] = solution_list[diff(phi,t,t,x,x)]
dictionary_of_solutions[diff(phi,t,x,x,x)] = solution_list[diff(phi,t,x,x,x)]

print('--------------------------------------------------------------------------')
for term in dictionary_of_solutions:
	print('--------------------------------------------------------------------------')
	print(showDiffNotation(term), ' = ')
	print('--------------------------------------------------------------------------')
	print(showDiffNotation(dictionary_of_solutions[term]))
	print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')


print('--------------------------------------------------------------------------')
print('Equation 1 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(equation1))
print('--------------------------------------------------------------------------')


equation2 = diff(dictionary_of_solutions[diff(phi,t,x,x,x)],t)
equation2 = expand(equation2).doit()

for term in dictionary_of_solutions:
	print('Equation 2 processing the substitution of the term ',term)
	equation2 = equation2.subs(term,dictionary_of_solutions[term])
	equation2 = expand(equation2).doit()

equation_list_4 = [
	expand(diff(phi,t,t,x,x,x) - equation2),
	expand(diff(phi,t,t) - dictionary_of_solutions[diff(phi,t,t)]),
	expand(diff(phi,t,t,x) - dictionary_of_solutions[diff(phi,t,t,x)]),
	expand(diff(phi,t,t,x,x) - dictionary_of_solutions[diff(phi,t,t,x,x)]),
	expand(diff(phi,t,x,x,x) - dictionary_of_solutions[diff(phi,t,x,x,x)])
	]
variable_list_4 = [
	diff(phi,t,t,x,x,x),
	diff(phi,t,t),
	diff(phi,t,t,x),
	diff(phi,t,t,x,x),
	diff(phi,t,x,x,x)
	]

solution_list = solve(equation_list_4, variable_list_4)
equation2 = solution_list[diff(phi,t,t,x,x,x)]
dictionary_of_solutions[diff(phi,t,t)] = solution_list[diff(phi,t,t)]
dictionary_of_solutions[diff(phi,t,t,x)] = solution_list[diff(phi,t,t,x)]
dictionary_of_solutions[diff(phi,t,t,x,x)] = solution_list[diff(phi,t,t,x,x)]
dictionary_of_solutions[diff(phi,t,x,x,x)] = solution_list[diff(phi,t,x,x,x)]

print('--------------------------------------------------------------------------')
for term in dictionary_of_solutions:
	print('--------------------------------------------------------------------------')
	print(showDiffNotation(term), ' = ')
	print('--------------------------------------------------------------------------')
	print(showDiffNotation(dictionary_of_solutions[term]))
	print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')


print('--------------------------------------------------------------------------')
print('Equation 2 = ')
print('--------------------------------------------------------------------------')
print(showDiffNotation(equation2))
print('--------------------------------------------------------------------------')

compatibility_test = equation2 - equation1
compatibility_test = expand(compatibility_test)
print('--------------------------------')
print('Compatibility test equation2 - equation1 =')
print('--------------------------------')
print('--------------------------------')
print(showDiffNotation(compatibility_test))
print('--------------------------------')


