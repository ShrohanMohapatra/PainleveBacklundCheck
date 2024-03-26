from sympy import *
import pprint

x, t, sigma  = symbols('x t sigma')
u2 = Function('u2')(x, t)
phi = Function('phi')(x, t)
eqn1 = diff(phi,x)*diff(phi,t)+u2*diff(phi,x)**2+4*sigma*diff(phi,x)*diff(phi,(x,3))-3*sigma*diff(phi,(x,2))**2
eqn2 = diff(phi,x,t)+u2*diff(phi,x,x)+sigma*diff(phi,x,x,x,x)
solution = solve([eqn1,eqn2],[diff(phi,t),diff(phi,x,x,x,x)])

phit_soln = solution[diff(phi,t)]
phix4_soln = solution[diff(phi,x,x,x,x)]
u2_t_soln = -u2*diff(u2,x)-sigma*diff(u2,(x,3))

dictionary_of_solutions = {} # Initializing the dictionary of solutions
dictionary_of_solutions[diff(phi,t)] = phit_soln
dictionary_of_solutions[diff(phi,(x, 4))] = phix4_soln
dictionary_of_solutions[diff(u2,t)] = u2_t_soln

equation_new = dictionary_of_solutions[diff(phi,(x,4))].subs(diff(phi,x,t),diff(dictionary_of_solutions[diff(phi,t)],x))
equation_new = expand(equation_new)
dictionary_of_solutions[diff(phi,(x,4))] = equation_new

# Here I am solving for diff(phi,(x,4)) .....
equation_new = together(solve(diff(phi,(x,4))-dictionary_of_solutions[diff(phi,(x,4))],diff(phi,(x,4)))[0]).doit()
dictionary_of_solutions[diff(phi,(x,4))] = equation_new

max_number_of_derivatives = 11

for k in range(5,max_number_of_derivatives):
	equation_new = diff(dictionary_of_solutions[diff(phi,(x,k-1))],x)
	numerator, denominator = fraction(together(equation_new).doit())
	numerator = expand(numerator)
	equation_new = together(numerator/denominator).doit()
	equation_new = equation_new.subs(diff(phi,(x,4)),dictionary_of_solutions[diff(phi,(x,4))])
	numerator, denominator = fraction(together(equation_new).doit())
	numerator = expand(numerator)
	assert ode_order(numerator, phi) == 3
	equation_new = together(numerator/denominator).doit()
	# print('--------------------------------')
	# print(equation_new)
	# print('--------------------------------')
	dictionary_of_solutions[diff(phi,(x,k))] = equation_new

for k in range(1, max_number_of_derivatives):
	if k == 1:
		equation_new = diff(dictionary_of_solutions[diff(u2,t)],x)
	else:
		equation_new = diff(dictionary_of_solutions[diff(u2,(x,k-1),t)],x)
	equation_new = expand(equation_new)
	# print('--------------------------------')
	# print(equation_new)
	# print('--------------------------------')
	dictionary_of_solutions[diff(u2,(x,k),t)] = equation_new

equation_new = diff(dictionary_of_solutions[diff(phi,t)],x)
equation_new = expand(equation_new)
equation_new = equation_new.subs(diff(phi,(x,4)),dictionary_of_solutions[diff(phi,(x,4))])
equation_new = expand(equation_new)
numerator, denominator = fraction(together(equation_new).doit())
numerator = expand(numerator)
assert ode_order(numerator, phi) == 3
dictionary_of_solutions[diff(phi,x,t)] = equation_new

equation_new = diff(dictionary_of_solutions[diff(phi,x,t)],x)
equation_new = expand(equation_new)
equation_new = equation_new.subs(diff(phi,(x,4)),dictionary_of_solutions[diff(phi,(x,4))])
equation_new = expand(equation_new)
numerator, denominator = fraction(together(equation_new).doit())
numerator = expand(numerator)
assert ode_order(numerator, phi) == 3
dictionary_of_solutions[diff(phi,x,x,t)] = equation_new

equation_new = diff(dictionary_of_solutions[diff(phi,x,x,t)],x)
equation_new = expand(equation_new)
equation_new = equation_new.subs(diff(phi,(x,4)),dictionary_of_solutions[diff(phi,(x,4))])
equation_new = expand(equation_new)
numerator, denominator = fraction(together(equation_new).doit())
numerator = expand(numerator)
assert ode_order(numerator, phi) == 3
dictionary_of_solutions[diff(phi,x,x,x,t)] = equation_new

# for term in dictionary_of_solutions:
# 	print('--------------------------------')
# 	print(term,':', dictionary_of_solutions[term])
# 	print('--------------------------------')

equation1 = diff(dictionary_of_solutions[diff(phi,x,x,x,t)],x)
equation1 = expand(equation1)
equation1 = equation1.subs(diff(phi,(x,4)),dictionary_of_solutions[diff(phi,(x,4))])
equation1 = expand(equation1)
# print('--------------------------------')
# print(equation1)
# print('--------------------------------')

equation2 = diff(dictionary_of_solutions[diff(phi,x,x,x,x)],t)
recursion_track = 0
while True:
	recursion_track = recursion_track + 1
	# print('--------------------------------')
	# print('Recursion number for substitution for time derivative of phi_{x,4}', recursion_track)
	# print('--------------------------------')
	for term in dictionary_of_solutions:
		equation2 = equation2.subs(term,dictionary_of_solutions[term])
		equation2 = together(equation2).doit()
		equation2 = expand(equation2)
	numerator, denominator = fraction(together(equation2).doit())
	numerator = expand(numerator)
	order_param = ode_order(numerator, phi)
	# print('--------------------------------')
	# print('ode_order(equation2, phi) = ', order_param)
	# print('--------------------------------')
	if order_param < 4: break

compatibility_test = equation2 - equation1
compatibility_test = expand(compatibility_test)
print('--------------------------------')
print(compatibility_test)
print('--------------------------------')
