from sympy import *
x, t = symbols('x t')
u = Function('u')(x,t)
pde_equation = u*diff(u,x)-diff(u,t,x,x)+diff(u,x)**2+diff(u,t)/diff(u,(x,3))
expression_terms = pde_equation.as_coefficients_dict()
list_of_differential_terms = []
for term, coefficient in expression_terms.items():
	print('term ->', term, ' type(term) ->', type(term))
	if isinstance(term, Derivative) and term not in list_of_differential_terms:
		list_of_differential_terms.append(term)
	elif isinstance(term, Mul):
		print('term ->', term,' term.args ->', term.args)
		for sub_term in term.args:
			print('sub_term ->', sub_term, ' type(sub_term) ->', type(sub_term), ' sub_term.args ->', sub_term.args)
			if isinstance(sub_term, Derivative) and sub_term not in list_of_differential_terms:
				list_of_differential_terms.append(sub_term)
			elif isinstance(sub_term, Pow):
				for sub_sub_term in sub_term.args:
					if isinstance(sub_sub_term, Derivative) and sub_sub_term not in list_of_differential_terms:
						list_of_differential_terms.append(sub_sub_term)
print(list_of_differential_terms)
max_order_derivative_in_t = 0
max_total_order = 0
term_seek = list_of_differential_terms[0]
for derivative_term in list_of_differential_terms:
	list_of_partial_derivatives = derivative_term.variables
	if t in list_of_partial_derivatives:
		count_of_t = 0
		for variable in list_of_partial_derivatives:
			if variable == t: count_of_t = count_of_t + 1
		if max_order_derivative_in_t <= count_of_t and max_total_order <= len(list_of_partial_derivatives):
			max_order_derivative_in_t = count_of_t
			max_total_order = len(list_of_partial_derivatives)
			term_seek = derivative_term
print(term_seek)

