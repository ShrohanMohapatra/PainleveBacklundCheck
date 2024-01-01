from sympy import *

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

x, t = symbols('x t')
TransformedSineGordonEquation = parse_expr('2*f(x,t)*diff(f(x,t),x,t)-2*diff(f(x,t),x)*diff(f(x,t),t)-f(x,t)**3+f(x,t)')
u0 = Function('u0')(x, t)
u1 = Function('u1')(x, t)
phi = Function('phi')(x, t)
f = Function('f')
u = u0/phi**2 + u1/phi
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

Equation2 = simplify(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 4))
Equation2 = expand(Equation2.subs(u0, u0_solution).doit())
# Equation2 = expand(Equation2.subs(u1, u1_solution).doit())
print('--------------------------------------------------------------------------')
print('Equation for 2nd power coefficient')
print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation2),' = 0')
print('--------------------------------------------------------------------------')

Equation3 = simplify(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 3))
Equation3 = expand(Equation3.subs(u0, u0_solution).doit())
# Equation3 = expand(Equation3.subs(u1, u1_solution).doit())
Equation3, _ = fraction(together(Equation3.doit()))
Equation3 = expand(Equation3)
print('--------------------------------------------------------------------------')
print('Equation for 3rd power coefficient')
print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation3),' = 0')
print('--------------------------------------------------------------------------')

Equation4 = simplify(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 2))
Equation4 = expand(Equation4.subs(u0, u0_solution).doit())
# Equation4 = expand(Equation4.subs(u1, u1_solution).doit())
Equation4, _ = fraction(together(Equation4.doit()))
Equation4 = expand(Equation4)
print('--------------------------------------------------------------------------')
print('Equation for 4th power coefficient')
print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation4),' = 0')
print('--------------------------------------------------------------------------')

Equation5 = simplify(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 1))
Equation5 = expand(Equation5.subs(u0, u0_solution).doit())
# Equation5 = expand(Equation5.subs(u1, u1_solution).doit())
Equation5, _ = fraction(together(Equation5.doit()))
Equation5 = expand(Equation5)
print('--------------------------------------------------------------------------')
print('Equation for 5th power coefficient')
print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation5),' = 0')
print('--------------------------------------------------------------------------')

Equation6 = simplify(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 0))
Equation6 = expand(Equation6.subs(u0, u0_solution).doit())
# Equation6 = expand(Equation6.subs(u1, u1_solution).doit())
Equation6, _ = fraction(together(Equation6.doit()))
Equation6 = expand(Equation6)
print('--------------------------------------------------------------------------')
print('Equation for 6th power coefficient')
print('--------------------------------------------------------------------------')
print('--------------------------------------------------------------------------')
print(showDiffNotation(Equation6),' = 0')
print('--------------------------------------------------------------------------')
