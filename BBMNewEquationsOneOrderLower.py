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

u = u0/phi**2 + u1/phi
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
u1_solution = solve(Equation1, u1)[0]
u1_solution = together(u1_solution).doit()
Equation1_residual = Equation1.subs(u1, u1_solution)
Equation1_residual, _ = fraction(together(Equation1_residual).doit())
Equation1_residual = expand(Equation1_residual)
print('--------------------------------------------------------------------------')
print('u1 = ', u1_solution)
print('--------------------------------------------------------------------------')
print('Equation1_residual = ', Equation1_residual)
print('--------------------------------------------------------------------------')

Equation2 = expand(BBMPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 3))
Equation2 = expand(Equation2.subs(u0, u0_solution).doit())
Equation2 = expand(Equation2.subs(u1, u1_solution).doit())
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
Equation3 = expand(Equation3.subs(u1, u1_solution).doit())
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
Equation4 = expand(Equation4.subs(u0, u0_solution).doit())
Equation4 = expand(Equation4.subs(u1, u1_solution).doit())
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
Equation5 = expand(Equation5.subs(u0, u0_solution).doit())
Equation5 = expand(Equation5.subs(u1, u1_solution).doit())
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

Equation6 = expand(BBMPainleveExpansion_polyInTermsOfPhi.coeff(phi, 1))
Equation6 = expand(Equation6.subs(u0, u0_solution).doit())
Equation6 = expand(Equation6.subs(u1, u1_solution).doit())
# Equation5_compatibility_test = Equation5
# Equation5 = expand(Equation5.subs(u2, u2_solution).doit())
Equation6, _ = fraction(together(Equation6).doit())
Equation6 = expand(Equation6)
# factors_list_of_equation5 = factor_list(Equation5)
# print('--------------------------------------------------------------------------')
# print(len(factors_list_of_equation5))
# print('--------------------------------------------------------------------------')
# Equation5 = factors_list_of_equation5[-1][-1][0]
print('--------------------------------------------------------------------------')
print('Equation for 6th power coefficient')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation6),' = 0')
print('--------------------------------------------------------------------------')