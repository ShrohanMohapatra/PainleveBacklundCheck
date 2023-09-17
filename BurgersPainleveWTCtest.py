from sympy import *
x, t, sigma = symbols('x t sigma')
phi, u0 = symbols('phi u0', cls = Function)
u = (phi(x,t)*(sigma*diff(phi(x, t), (x, 2)) - diff(phi(x, t), t)) - 2*sigma*diff(phi(x, t), x)**2)/phi(x, t)/diff(phi(x, t), x)
u2 = (sigma*diff(phi(x, t), (x, 2)) - diff(phi(x, t), t))/diff(phi(x, t), x) 
original_equation = diff(u, t) + u*diff(u, x) - sigma*diff(u, (x, 2))
original_equation, _ = fraction(together(original_equation.doit()))
original_equation = expand(original_equation)
print('------------------------------------------------------------')
print('original_equation')
print('------------------------------------------------------------')
print(original_equation)
print('------------------------------------------------------------')
backlund_equation = diff(u2, t) + u2*diff(u2, x) - sigma*diff(u2, (x, 2))
backlund_equation, _ = fraction(together(backlund_equation.doit()))
backlund_equation = expand(backlund_equation)
print('------------------------------------------------------------')
print('backlund_equation')
print('------------------------------------------------------------')
print(backlund_equation)
print('------------------------------------------------------------')
################## Let's preprocess the original equation a tad bit
terms_from_the_original_equation = original_equation.as_ordered_terms()
print('------------------------------------------------------------')
print('terms_from_the_original_equation')
print('------------------------------------------------------------')
print(terms_from_the_original_equation)
print('------------------------------------------------------------')
########## We will extract the maximum order term of the original equation
order_of_original_equation = ode_order(original_equation, phi)
for individual_term in terms_from_the_original_equation:
	if ode_order(individual_term, phi) == order_of_original_equation:
		factorList = factor_list(individual_term)
		for factor in factorList[1]:
			if ode_order(factor[0], phi) == order_of_original_equation:
				actual_term = factor[0]**factor[1]
				full_term = individual_term
				break
		break
actual_max_order_term = actual_term
coefficient = simplify(full_term/actual_max_order_term)
print('------------------------------------------------------------')
print('Order of the original equation')
print('------------------------------------------------------------')
print(order_of_original_equation)
print('------------------------------------------------------------')
print('------------------------------------------------------------')
print('Maximum order term occuring in the original equation')
print('------------------------------------------------------------')
print(actual_max_order_term)
print('------------------------------------------------------------')
print('------------------------------------------------------------')
print('The full term of maximum order in the original equation')
print('------------------------------------------------------------')
print('(', actual_max_order_term, ')* (', coefficient, ') =', full_term)
print('------------------------------------------------------------')
################## Let's preprocess the Backlund equation a little bit more
terms_from_the_backlund_equation = backlund_equation.as_ordered_terms()
print('------------------------------------------------------------')
print('terms_from_the_backlund_equation')
print('------------------------------------------------------------')
print(terms_from_the_backlund_equation)
print('------------------------------------------------------------')
order_of_backlund_equation = ode_order(backlund_equation, phi)
for individual_term in terms_from_the_backlund_equation:
	if ode_order(individual_term, phi) == order_of_backlund_equation:
		factorList = factor_list(individual_term)
		for factor in factorList[1]:
			if ode_order(factor[0], phi) == order_of_backlund_equation:
				actual_term = factor[0]**factor[1]
				full_term = individual_term
				break
		break
actual_max_order_term = actual_term
coefficient = simplify(full_term/actual_max_order_term)
print('------------------------------------------------------------')
print('Order of the Backlund equation')
print('------------------------------------------------------------')
print(order_of_backlund_equation)
print('------------------------------------------------------------')
print('------------------------------------------------------------')
print('Maximum order term occuring in the Backlund equation')
print('------------------------------------------------------------')
print(actual_max_order_term)
print('------------------------------------------------------------')
print('------------------------------------------------------------')
print('The full term of maximum order in the Backlund equation')
print('------------------------------------------------------------')
print('(', actual_max_order_term, ')* (', coefficient, ') =', full_term)
print('------------------------------------------------------------')
############# Let's execute a while