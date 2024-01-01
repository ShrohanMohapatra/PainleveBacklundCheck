from sympy import *
from random import choice
x, t, sigma, kappa = symbols('x t sigma kappa')
def ContextSensitiveDifferentialCommutationCompatibilitySystemPDEs(systemsOfPDEs, contextualPDEs, phi, list_of_u_functions):
	assert len(contextualPDEs) == len(list_of_u_functions)
	num_kappa = len(systemsOfPDEs)
	n = len(contextualPDEs)
	solution_dictionary_for_all_higher_derivs = {}
	list_of_all_highest_derivative_terms_in_uk = []
	list_of_all_highest_derivative_terms_in_phi = []
	for k in range(n):
		equation_for_u_k = contextualPDEs[k]
		u_function_k = list_of_u_functions[k]
		list_of_terms = equation_for_u_k.args
		list_of_differential_terms = []
		for m in range(len(list_of_terms)):
			list_of_sub_terms_in_a_term = list_of_terms[m].args
			for differential_term in list_of_sub_terms_in_a_term:
				if differential_term not in list_of_differential_terms:
					list_of_differential_terms.append(differential_term)
		expression_terms = equation_for_u_k.as_coefficients_dict()
		list_of_differential_terms = []
		for term, coefficient in expression_terms.items():
			if isinstance(term, Derivative) and term not in list_of_differential_terms:
				list_of_differential_terms.append(term)
			elif isinstance(term, Mul):
				for sub_term in term.args:
					if isinstance(sub_term, Derivative) and sub_term not in list_of_differential_terms:
						list_of_differential_terms.append(sub_term)
		max_order_derivative_in_t = 0
		max_total_order = 0
		highest_order_term_of_u_function_k = list_of_differential_terms[0]
		for derivative_term in list_of_differential_terms:
			list_of_partial_derivatives = derivative_term.variables
			if t in list_of_partial_derivatives:
				count_of_t = 0
				for variable in list_of_partial_derivatives:
					if variable == t: count_of_t = count_of_t + 1
				if max_order_derivative_in_t <= count_of_t and max_total_order <= len(list_of_partial_derivatives):
					max_order_derivative_in_t = count_of_t
					max_total_order = len(list_of_partial_derivatives)
					highest_order_term_of_u_function_k = derivative_term
		solution_of_u_function_k = solve(equation_for_u_k, highest_order_term_of_u_function_k)[0]
		solution_of_u_function_k = expand(solution_of_u_function_k)
		solution_dictionary_for_all_higher_derivs[highest_order_term_of_u_function_k] = solution_of_u_function_k
		print('---------------------------------------------------')
		print(solution_dictionary_for_all_higher_derivs)
		print('---------------------------------------------------')
		if highest_order_term_of_u_function_k not in list_of_all_highest_derivative_terms_in_uk:
			list_of_all_highest_derivative_terms_in_uk.append(highest_order_term_of_u_function_k)
	for k in range(num_kappa):
		equation_k_from_the_phi_system = systemsOfPDEs[k]
		list_of_terms = equation_k_from_the_phi_system.args
		list_of_differential_terms_from_phi = []
		for m in range(len(list_of_terms)):
			list_of_sub_terms_in_a_term = list_of_terms[m].args
			for differential_term in list_of_sub_terms_in_a_term:
				if differential_term not in list_of_differential_terms_from_phi:
					list_of_differential_terms_from_phi.append(differential_term)
		if k != num_kappa - 1:
			expression_terms = equation_k_from_the_phi_system.as_coefficients_dict()
			list_of_differential_terms = []
			for term, coefficient in expression_terms.items():
				if isinstance(term, Derivative) and term not in list_of_differential_terms:
					list_of_differential_terms.append(term)
				elif isinstance(term, Mul):
					for sub_term in term.args:
						if isinstance(sub_term, Derivative) and sub_term not in list_of_differential_terms:
							list_of_differential_terms.append(sub_term)
			max_order_derivative_in_t = 0
			max_total_order = 0
			highest_order_term_of_phi_from_equation_k = list_of_differential_terms[0]
			for derivative_term in list_of_differential_terms:
				list_of_partial_derivatives = derivative_term.variables
				if t in list_of_partial_derivatives:
					count_of_t = 0
					for variable in list_of_partial_derivatives:
						if variable == t: count_of_t = count_of_t + 1
					if max_order_derivative_in_t <= count_of_t and max_total_order <= len(list_of_partial_derivatives):
						max_order_derivative_in_t = count_of_t
						max_total_order = len(list_of_partial_derivatives)
						highest_order_term_of_phi_from_equation_k = derivative_term
			if highest_order_term_of_phi_from_equation_k in solution_dictionary_for_all_higher_derivs:
				while True:
					print('Detection stage')
					new_equation = expand(equation_k_from_the_phi_system.subs(highest_order_term_of_phi_from_equation_k,
						solution_dictionary_for_all_higher_derivs[highest_order_term_of_phi_from_equation_k]).doit())
					order_of_phi_from_new_equation = ode_order(new_equation, phi)
					list_of_terms_from_the_new_equation = new_equation.args
					for term in list_of_terms_from_the_new_equation:
						if isinstance(term, Derivative) and ode_order(term, phi) == order_of_phi_from_new_equation:
							highest_order_term_of_phi_from_equation_k = term
							break
						else:
							list_of_sub_terms = term.args
							for sub_term in list_of_sub_terms:
								if isinstance(sub_term, Derivative) and ode_order(sub_term, phi) == order_of_phi_from_new_equation:
									highest_order_term_of_phi_from_equation_k = sub_term
									break
								else:
									list_of_sub_sub_terms = sub_term.args
									for sub_sub_term in list_of_sub_sub_terms:
										if isinstance(sub_sub_term, Derivative) and ode_order(sub_sub_term, phi) == order_of_phi_from_new_equation:
											highest_order_term_of_phi_from_equation_k = sub_sub_term
											break
					if highest_order_term_of_phi_from_equation_k in solution_dictionary_for_all_higher_derivs: pass
					else:
						new_equation, _ = fraction(together(new_equation).doit())
						new_equation = expand(new_equation)
						solution_of_highest_order_term_of_phi_from_equation_k = expand(solve(new_equation, highest_order_term_of_phi_from_equation_k)[0])
						solution_of_highest_order_term_of_phi_from_equation_k = expand(solution_of_highest_order_term_of_phi_from_equation_k)
						solution_dictionary_for_all_higher_derivs[highest_order_term_of_phi_from_equation_k] = solution_of_highest_order_term_of_phi_from_equation_k
						print('---------------------------------------------------')
						print(solution_dictionary_for_all_higher_derivs)
						print('---------------------------------------------------')
						if highest_order_term_of_phi_from_equation_k not in list_of_all_highest_derivative_terms_in_phi:
							list_of_all_highest_derivative_terms_in_phi.append(highest_order_term_of_phi_from_equation_k)
						break
			else:
				solution_of_highest_order_term_of_phi_from_equation_k = solve(equation_k_from_the_phi_system, highest_order_term_of_phi_from_equation_k)[0]
				solution_of_highest_order_term_of_phi_from_equation_k = expand(solution_of_highest_order_term_of_phi_from_equation_k)
				solution_dictionary_for_all_higher_derivs[highest_order_term_of_phi_from_equation_k] = solution_of_highest_order_term_of_phi_from_equation_k
				print('---------------------------------------------------')
				print(solution_dictionary_for_all_higher_derivs)
				print('---------------------------------------------------')
				if highest_order_term_of_phi_from_equation_k not in list_of_all_highest_derivative_terms_in_phi:
					list_of_all_highest_derivative_terms_in_phi.append(highest_order_term_of_phi_from_equation_k)
		else:
			order_of_phi_from_equation_k = ode_order(equation_k_from_the_phi_system, phi)
			for term in list_of_differential_terms_from_phi:
				if ode_order(term, phi) == order_of_phi_from_equation_k:
					highest_order_term_of_phi_from_equation_k = term
					break
			if highest_order_term_of_phi_from_equation_k in solution_dictionary_for_all_higher_derivs:
				while True:
					print('Detection stage')
					new_equation = expand(equation_k_from_the_phi_system.subs(highest_order_term_of_phi_from_equation_k,
						solution_dictionary_for_all_higher_derivs[highest_order_term_of_phi_from_equation_k]).doit())
					order_of_phi_from_new_equation = ode_order(new_equation, phi)
					list_of_terms_from_the_new_equation = new_equation.args
					for term in list_of_terms_from_the_new_equation:
						if isinstance(term, Derivative) and ode_order(term, phi) == order_of_phi_from_new_equation:
							highest_order_term_of_phi_from_equation_k = term
							break
						else:
							list_of_sub_terms = term.args
							for sub_term in list_of_sub_terms:
								if isinstance(sub_term, Derivative) and ode_order(sub_term, phi) == order_of_phi_from_new_equation:
									highest_order_term_of_phi_from_equation_k = sub_term
									break
								else:
									list_of_sub_sub_terms = sub_term.args
									for sub_sub_term in list_of_sub_sub_terms:
										if isinstance(sub_sub_term, Derivative) and ode_order(sub_sub_term, phi) == order_of_phi_from_new_equation:
											highest_order_term_of_phi_from_equation_k = sub_sub_term
											break
					if highest_order_term_of_phi_from_equation_k in solution_dictionary_for_all_higher_derivs: pass
					else:
						new_equation, _ = fraction(together(new_equation).doit())
						new_equation = expand(new_equation)
						print('highest_order_term_of_phi_from_equation_k -> ', highest_order_term_of_phi_from_equation_k)
						print('new_equation -> ', new_equation)
						solution_of_highest_order_term_of_phi_from_equation_k = solve(new_equation, highest_order_term_of_phi_from_equation_k)[0]
						solution_of_highest_order_term_of_phi_from_equation_k = expand(solution_of_highest_order_term_of_phi_from_equation_k)
						solution_dictionary_for_all_higher_derivs[highest_order_term_of_phi_from_equation_k] = solution_of_highest_order_term_of_phi_from_equation_k
						print('---------------------------------------------------')
						print(solution_dictionary_for_all_higher_derivs)
						print('---------------------------------------------------')
						if highest_order_term_of_phi_from_equation_k not in list_of_all_highest_derivative_terms_in_phi:
							list_of_all_highest_derivative_terms_in_phi.append(highest_order_term_of_phi_from_equation_k)
						break
			else:
				solution_of_highest_order_term_of_phi_from_equation_k = solve(equation_k_from_the_phi_system, highest_order_term_of_phi_from_equation_k)[0]
				solution_of_highest_order_term_of_phi_from_equation_k = expand(solution_of_highest_order_term_of_phi_from_equation_k)
				solution_dictionary_for_all_higher_derivs[highest_order_term_of_phi_from_equation_k] = solution_of_highest_order_term_of_phi_from_equation_k
				print('---------------------------------------------------')
				print(solution_dictionary_for_all_higher_derivs)
				print('---------------------------------------------------')
				if highest_order_term_of_phi_from_equation_k not in list_of_all_highest_derivative_terms_in_phi:
					list_of_all_highest_derivative_terms_in_phi.append(highest_order_term_of_phi_from_equation_k)
	# From here on I have gotten the basic preparatory lists maintained, now I am following the rest of it from the pseudocode of the paper
	S_equations, S_variables = [], []
	for variable, substitution in solution_dictionary_for_all_higher_derivs.items():
		equation = expand(variable - substitution)
		# equation, _ = fraction(together(equation).doit())
		# equation = expand(equation)
		if variable not in S_variables:
			S_variables.append(variable)
			S_equations.append(equation)
	print('---------------------------------------------------')
	print('list_of_all_highest_derivative_terms_in_phi = ')
	print('---------------------------------------------------')
	print(list_of_all_highest_derivative_terms_in_phi)
	print('---------------------------------------------------')
	for k in range(len(S_equations)):
		print('---------------------------------------------------')
		print('S_equations[k] = ')
		print('---------------------------------------------------')
		print(S_equations[k])
		print('---------------------------------------------------')
	print('---------------------------------------------------')
	print('List of variables ->')
	print('---------------------------------------------------')
	print(S_variables)
	print('---------------------------------------------------')
	
	for k in range(len(S_variables)):
		variable = S_variables[k]
		equation = S_equations[k]
		solution = expand(variable-equation)
		# solution, _ = fraction(together(solution).doit())
		# solution = expand(solution)
		solution_dictionary_for_all_higher_derivs[variable] = solution
		print('---------------------------------------------------')
		print(S_equations)
		print('---------------------------------------------------')
		print(S_variables)
		print('---------------------------------------------------')

	list_of_A_orders, list_of_B_orders = [], []
	for term in list_of_all_highest_derivative_terms_in_phi:
		list_of_partial_derivatives = term.variables
		if t in list_of_partial_derivatives:
			count_of_t, count_of_x = 0, 0
			for variable in list_of_partial_derivatives:
				if variable == t: count_of_t = count_of_t + 1
				elif variable == x: count_of_x = count_of_x + 1
		list_of_A_orders.append(count_of_x)
		list_of_B_orders.append(count_of_t)
	A_max, A_min, B_max, B_min = max(list_of_A_orders), min(list_of_A_orders), max(list_of_B_orders), min(list_of_B_orders)
	list_of_C_orders, list_of_D_orders = [], []
	for term in list_of_all_highest_derivative_terms_in_uk:
		list_of_partial_derivatives = term.variables
		if t in list_of_partial_derivatives:
			count_of_t, count_of_x = 0, 0
			for variable in list_of_partial_derivatives:
				if variable == t: count_of_t = count_of_t + 1
				elif variable == x: count_of_x = count_of_x + 1
		list_of_C_orders.append(count_of_x)
		list_of_D_orders.append(count_of_t)
	# print('len(list_of_all_highest_derivative_terms_in_uk) = ', len(list_of_all_highest_derivative_terms_in_uk))
	# print('len(list_of_C_orders) = ',len(list_of_C_orders))
	# print('len(list_of_D_orders) = ',len(list_of_D_orders))
	# print('num_kappa = ', num_kappa)
	for v in range(len(list_of_all_highest_derivative_terms_in_uk)):
		Cv = list_of_C_orders[v]
		Dv = list_of_D_orders[v]
		uv = list_of_all_highest_derivative_terms_in_uk[v]
		gv = solution_dictionary_for_all_higher_derivs[uv]
		for i in range(1, A_max-A_min+1):
			for j in range(1, B_max-B_min+1):
				equation = expand(diff(uv,(x,i+Cv),(t,j+Dv))-diff(gv,(x,i),(t,j)))
				# equation, _ = fraction(together(equation).doit())
				# equation = expand(equation)
				variable = diff(uv,(x,i+Cv),(t,j+Dv))
				if equation not in S_equations: S_equations.append(equation)
				if variable not in S_variables: S_variables.append(variable)
				print('---------------------------------------------------')
				print(S_equations)
				print('---------------------------------------------------')
				print(S_variables)
				print('---------------------------------------------------')
	for w in range(len(list_of_A_orders)):
		Aw = list_of_A_orders[w]
		Bw = list_of_B_orders[w]
		phi_term_w = list_of_all_highest_derivative_terms_in_phi[w]
		fw = solution_dictionary_for_all_higher_derivs[phi_term_w]
		for derivsInX in range(1, A_max-Aw+1):
			for derivsInT in range(1, B_max-Bw+1):
				equation = expand(diff(phi,(x,Aw+derivsInX),(t,Bw+derivsInT))-diff(fw,(x,derivsInX),(t,derivsInT)))
				# equation, _ = fraction(together(equation).doit())
				# equation = expand(equation)
				variable = diff(phi,(x,Aw+derivsInX),(t,Bw+derivsInT))
				if equation not in S_equations: S_equations.append(equation)
				if variable not in S_variables: S_variables.append(variable)
				print('---------------------------------------------------')
				print(S_equations)
				print('---------------------------------------------------')
				print(S_variables)
				print('---------------------------------------------------')
	# Now having set all the equations and the variables I will have to be a little careful
	# setting up the equation solver both in terms of the data structure and also
	# the discontinuation of the old solutions.
	
	while True:
		S_new_equations = []
		solution_to_data_structure = list(nonlinsolve(S_equations, S_variables))
		try:
			for variable, solution in solution_to_data_structure.items():
				if variable in S_variables:
					equation = expand(variable-solution)
					# equation, _ = fraction(together(equation).doit())
					# equation = expand(equation)
					if equation not in S_new_equations:
						S_new_equations.append(equation)
		except:
			new_solution_set = solution_to_data_structure[0]
			for k in range(len(S_variables)):
				variable = S_variables[k]
				solution = new_solution_set[k]
				equation = expand(variable-solution)
				# equation, _ = fraction(together(equation).doit())
				# equation = expand(equation)
				if equation not in S_new_equations:
					S_new_equations.append(equation)
		flag = len(list(set(S_new_equations) - set(S_equations))) == 0
		S_equations = S_new_equations
		print('---------------------------------------------------')
		print(S_equations)
		print('---------------------------------------------------')
		print(S_variables)
		print('---------------------------------------------------')
		if flag: break

	W_set, R_set = [None for k in range(n)], [None for k in range(n)]
	flag1 = True
	for k in range(n-1):
		Ak, Ak_plus_1 = list_of_A_orders[k], list_of_A_orders[k+1]
		Bk, Bk_plus_1 = list_of_B_orders[k], list_of_B_orders[k+1]
		phi_term_k = list_of_all_highest_derivative_terms_in_phi[k]
		fk = solution_dictionary_for_all_higher_derivs[phi_term_k]
		phi_term_k_plus_1 = list_of_all_highest_derivative_terms_in_phi[k+1]
		fk_plus_1 = solution_dictionary_for_all_higher_derivs[phi_term_k_plus_1]
		W_set[k] = expand(diff(fk, (x,A_max-Ak), (t, B_max-Bk))-diff(fk_plus_1, (x,A_max-Ak_plus_1), (t, B_max-Bk_plus_1)))
		W_set[k], _ = fraction(together(W_set[k]).doit())
		W_set[k] = expand(W_set[k])
		for variable, solution in solution_dictionary_for_all_higher_derivs.items():
			W_set[k] = expand(W_set[k].subs(variable, solution))
			W_set[k], _ = fraction(together(W_set[k]).doit())
			W_set[k] = expand(W_set[k])
		# Let us see if we can extract the nth roots (for e.g. square roots etc.) in W_set[k] and make it easier ....
		target_root_order = 2
		while target_root_order < 10:
			print('target_root_order = ',target_root_order)
			nth_root_instances = [term for term in W_set[k].find(lambda x: isinstance(x, Pow) and x.exp == 1/target_root_order)]
			has_nth_root = len(nth_root_instances) > 0
			if has_nth_root:
				reference_expression = expand(W_set[k]**target_root_order)
				radical = nth_root_instances[0]
				W_set[k] = resultant(W_set[k], reference_expression, radical)
				break
			target_root_order = target_root_order + 1
		flag1 = flag1 and W_set[k] == 0
		print('---------------------------------------------------')
		print(W_set[k])
		print('---------------------------------------------------')
	A_n_minus_1, A_0 = list_of_A_orders[n-1], list_of_A_orders[0]
	B_n_minus_1, B_0 = list_of_B_orders[n-1], list_of_B_orders[0]
	phi_term_n_minus_1 = list_of_all_highest_derivative_terms_in_phi[n-1]
	f_n_minus_1 = solution_dictionary_for_all_higher_derivs[phi_term_n_minus_1]
	phi_term_0 = list_of_all_highest_derivative_terms_in_phi[0]
	f_0 = solution_dictionary_for_all_higher_derivs[phi_term_0]
	W_set[n-1] = expand(diff(f_n_minus_1, (x,A_max-A_n_minus_1), (t, B_max-B_n_minus_1))-diff(f_0, (x,A_max-A_0), (t, B_max-B_0)))
	W_set[n-1], _ = fraction(together(W_set[n-1]).doit())
	W_set[n-1] = expand(W_set[n-1])
	for variable, solution in solution_dictionary_for_all_higher_derivs.items():
		W_set[n-1] = expand(W_set[n-1].subs(variable, solution))
		W_set[n-1], _ = fraction(together(W_set[n-1]).doit())
		W_set[n-1] = expand(W_set[n-1])
	target_root_order = 2
	while target_root_order < 10:
		nth_root_instances = [term for term in W_set[n-1].find(lambda x: isinstance(x, Pow) and x.exp == 1/target_root_order)]
		has_nth_root = len(nth_root_instances) > 0
		if has_nth_root:
			reference_expression = expand(W_set[n-1]**target_root_order)
			radical = nth_root_instances[0]
			W_set[n-1] = resultant(W_set[n-1], reference_expression, radical)
			break
		target_root_order = target_root_order + 1
	print('---------------------------------------------------')
	print(W_set[n-1])
	print('---------------------------------------------------')
	flag1 = flag1 and W_set[n-1] == 0
	flag2 = True
	phi_term_random = choice(list_of_all_highest_derivative_terms_in_phi)
	for k in range(n-1):
		R_set[k] = expand(resultant(
			W_set[k],W_set[k+1],\
			phi_term_random
			))
		flag2 = flag2 and R_set[k] == 0
		print('---------------------------------------------------')
		print(R_set[k])
		print('---------------------------------------------------')
	R_set[n-1] = expand(resultant(
		W_set[n-1],W_set[0],\
		phi_term_random
		))
	flag2 = flag2 and R_set[n-1] == 0
	print('---------------------------------------------------')
	print(R_set[n-1])
	print('---------------------------------------------------')
	print('---------------------------------------------------')
	print(flag1)
	print('---------------------------------------------------')
	print('---------------------------------------------------')
	print(flag2)
	print('---------------------------------------------------')
	return (flag1 or flag2)


