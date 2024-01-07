from kivy.lang import Builder
from kivymd.app import MDApp
from kivy.uix.boxlayout import BoxLayout
from kivymd.uix.dialog import MDDialog
from kivy.clock import Clock
from sympy import *
from time import sleep
import threading
from collections import Counter
from random import choice

x, t, sigma, b, kappa = symbols('x t sigma b kappa')
f = Function('f')

def showDiffNotation(expr):
	functions = expr.atoms(Function)
	reps = {}
	for fun in functions:
		try:            
			reps[fun] = Symbol(fun.name)
		except AttributeError:
			continue
	dreps = [(deriv, Symbol(deriv.expr.subs(reps).name + "_{," +
							''.join(par.name for par in deriv.variables) + "}"))  \
			for deriv in expr.atoms(Derivative)]
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

KV = '''
<Content>
	orientation: "vertical"
	spacing: "12dp"
	size_hint_y: None
	height: "120dp"

	MDTextField:
		id: process_displayer
		hint_text: "Process will be shown here"
		font_name: "Georgia"
		font_size: 30
		text: "Process will be shown here"

	MDRectangleFlatButton:
		id: close_button
		font_name: "Georgia"
		font_size: 30
		text: "Close"
		on_release: app.close_handler()

MDScreen:

	MDLabel:
		id: main_handler
		multiline: True
		text: "Welcome to the \\nPainleve-Backlund \\n check app"
		halign: "center"
		font_name: "Georgia"
		font_size: 30
		pos_hint: {"center_x": 0.2, "center_y": 0.9}

	MDTextField:
		id: input_pde
		text: "Enter the PDE"
		font_name: "Georgia"
		font_size: 30
		halign: "center"
		pos_hint: {"center_x": 0.65, "center_y": 0.9}
		size_hint: 0.5, 0.125

	MDTextField:
		multiline: True
		id: integrability_handler
		font_name: "Georgia"
		font_size: 30
		text: "The result of the integrability test will be shown here."
		halign: "center"
		pos_hint: {"center_x": 0.5, "center_y": 0.4}
		size_hint: 0.3, 0.3

	MDRectangleFlatButton:
		id: integrability_window
		text: "Integrability test"
		font_name: "Georgia"
		font_size: 30
		theme_text_color: "Custom"
		pos_hint: {"center_x": 0.4, "center_y": 0.1}
		on_release: app.show_dialog()

'''

class Content(BoxLayout):
    pass

class PainleveBacklundCheckApp(MDApp):
	dialog = None
	string_storer = ''
	context_free_compatibility_test = True
	context_sensitive_compatibility_test = True
	compatibility_test = True
	backlund_transform = None
	associated_conditions = None

	def build(self):
		self.root = Builder.load_string(KV)
		self.theme_cls.theme_style = 'Dark'
		self.theme_cls.primary_palette = 'Pink'
		return self.root

	def on_start(self):
		self.fps_monitor_start()

	def scheduleOnceTester(self, dt):
		self.dialog.content_cls.ids.process_displayer.text = self.string_storer+'\n'

	def scheduleFinalizer(self, dt):
		self.dialog.content_cls.ids.process_displayer.text = self.string_storer+'\n'

	def totalTableOfCoeffsForPoly(self, expr, zeta):
		generator_list = [1/zeta, zeta]
		for m in range(10):
			for n in range(10):
				if m == 0 and n == 0: pass
				elif m == 0 and n == 1:
					generator_list.append(Derivative(zeta, t))
				elif m == 1 and n == 0:
					generator_list.append(Derivative(zeta, x))
				elif m == 1 and n == 1:
					generator_list.append(Derivative(zeta, t, x))
					# generator_list.append(Derivative(zeta, t, x))
				elif m == 1 and n > 1:
					generator_list.append(Derivative(zeta, (t, n), x))
					# generator_list.append(Derivative(zeta, (t, n), x))
				elif n == 1 and m > 1:
					generator_list.append(Derivative(zeta, t, (x, m)))
					# generator_list.append(Derivative(zeta, (x, m), t))
				elif m == 0 and n > 1:
					generator_list.append(Derivative(zeta, (t, n)))
				elif n == 0 and m > 1:
					generator_list.append(Derivative(zeta, (x, m)))
				elif m > 1 and n > 1:
					generator_list.append(
						Derivative(zeta, (t, n), (x, m))
						)
					# generator_list.append(
					# 	Derivative(zeta, (t, n), (x, m))
					# 	)
		print(zeta)
		print(generator_list)
		polynomial = Poly(expr, gens=generator_list)
		degree_list = polynomial.degree_list()
		# print(degree_list)
		min_deg, max_deg = degree_list[0], degree_list[1]
		min_deg = -min_deg
		degree_coeff_tuple = [
			[k, simplify(expr.coeff(zeta, k))]
			for k in range(min_deg,max_deg+1)]
		return degree_coeff_tuple
	def PainlevePDE_withBacklundTransform(
		self,
		function_PDE, # Input PDE
		x, # The input variable used for the spatial variable
		t, # The input variable used for the temporal variable
		alpha # The power balancing exponent
		):
		f = Function('f')(x, t)
		j = Symbol('j')
		phi = Function('phi')
		M = int(-alpha + 3) # Maximum number of terms
		U = [Function('U'+str(k)) for k in range(M+1)]
		u = 0
		for j in range(M+1):
			u = u + phi(x, t)**(j+alpha)*U[j](x, t)
		# print(u, function_PDE)
		orig_eqn = expand(function_PDE.subs(f, u).doit())
		# print(orig_eqn)
		totalTable = self.totalTableOfCoeffsForPoly(orig_eqn, phi(x, t))
		# print(totalTable)
		U1 = [U[k](x, t) for k in range(M+1)]
		U_final = [0 for k in range(M+1)]
		for k in range(len(totalTable)):
			# print(U_final)
			index, equation = totalTable[k]
			if equation == 0 or k >= M+1: continue
			if k == 0:
				try:
					soln = pdsolve(equation, U[k](x, t))
				except:
					soln = solve(equation, U[k](x, t))
			else:
				try:
					soln = pdsolve(
						equation.subs(
							[(U[m](x,t), U_final[m]) for m in range(k)]),
						U[k](x, t)
					)
				except:
					soln = solve(
						equation.subs(
							[(U[m](x,t), U_final[m]) for m in range(k)]),
						U[k](x, t)
					)
			if len(soln) == 0: break
			else:
				for m in range(len(soln)):
					if soln[m] == 0: continue
					else:
						U_final[k] = simplify(soln[m])
						break
		# Final expression for the Backlund transform
		backlund_transform = u.subs([
			(U[m](x,t), U_final[m]) if m<(-alpha) else \
			((U[m](x,t), U[m](x,t)) if m==(-alpha) else (U[m](x,t), 0))\
			for m in range(M+1)
			])
		# Compatibility conditions
		compatibility_Conditions = [
			[
				simplify(totalTable[k][1].subs([
				(U[m](x,t), U_final[m]) if m < cutoff_index else \
				((U[m](x,t), U[m](x,t)) if cutoff_index <= m <= (-alpha) else (U[m](x,t), 0))\
				for m in range(M+1)
				]).doit()) \
				for k in range(len(totalTable))
			] for cutoff_index in range(1, (-alpha)+1)
		]
		refined_compatibility_Conditions = [[] for cutoff_index in range(1, (-alpha)+1)]
		for cutoff_index in range(1, (-alpha)+1):
			for k in range(len(totalTable)):
				if compatibility_Conditions[cutoff_index-1][k] != 0:
					refined_compatibility_Conditions[cutoff_index-1].append(
						compatibility_Conditions[cutoff_index-1][k]
						)
		return [backlund_transform, refined_compatibility_Conditions]
	def ultimatePainleve(
		self,
		function_PDE, # Input PDE
		x, # The input variable used for the spatial variable
		t # The input variable used for the temporal variable
		):
		f = Function('f')(x, t)
		range_of_exponents = [-1, -2, -3, -4, -5]
		flag = False
		for alpha in range_of_exponents:
			M = int(-alpha + 3) # Maximum number of terms
			U = [Function('U'+str(k)) for k in range(M+1)]
			backlund_transform, compatibility_Conditions = self.PainlevePDE_withBacklundTransform(function_PDE, x, t, alpha)
			# print(alpha, compatibility_Conditions, backlund_transform)
			# print(alpha, len(compatibility_Conditions))
			if len(compatibility_Conditions[-1]) < 2: continue
			else:
				flag = expand(compatibility_Conditions[-1][-1].subs(U[-alpha](x, t),f).doit()) == expand(function_PDE)
				# print(compatibility_Conditions[-1][-1].subs(U[-alpha](x, t),f).doit(), expand(function_PDE))
				# print(alpha, flag)
				if not flag: continue
				else: break
		# print(alpha, compatibility_Conditions, backlund_transform)
		if alpha == -4 and len(compatibility_Conditions[-1]) == 1 and backlund_transform == U[4](x,t):
			return (None, backlund_transform, compatibility_Conditions)
		for k in range(len(compatibility_Conditions[-1])):
			if str(compatibility_Conditions[-1][k]).find(str(U[-alpha](x, t)))==-1:
				return (None, backlund_transform, [])				
		if not flag:
			alpha, compatibility_Conditions = None, []
		return (alpha, backlund_transform, compatibility_Conditions)

	def simplificationPDE(self, inputPDE1, inputPDE2, func):
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

		
		self.string_storer = str(inputPDE1)+' = 0'
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = str(inputPDE2)+' = 0'
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = 'The maximum order term in equation 1 is gamma_1 = '+str(max_order_term1[2])
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = 'The maximum order term in equation 2 is gamma_2 = '+str(max_order_term2[2])
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = 'A1 = '+str(A1)
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = 'A2 = '+str(A2)
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = 'B1 = '+str(B1)
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = 'B2 = '+str(B2)
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		
		allAvailableVariables = func.free_symbols

		varsFromPDE1 = list(max_order_term1[2].variables) if max_order_term1[3] != 0 else []
		varsFromPDE2 = list(max_order_term2[2].variables) if max_order_term2[3] != 0 else []

		self.string_storer = str(varsFromPDE1)
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = str(varsFromPDE2)
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		
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

		self.string_storer = str(varsFromPDE3)
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = str(varsFromPDE4)
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		
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

		self.string_storer = str(inputPDE3)
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = str(inputPDE4)
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		
		expression = inputPDE3 - inputPDE4
		expression, _ = fraction(together(expression.doit()))
		expression = expand(expression)

		self.string_storer = str(expression)
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		
		return expression

	def ContextFreeDifferentialCommutationCompatibilitySystemPDEs(self, systemsOfPDEs, func, depth_max = 20, depth = 0):

		for k in range(len(systemsOfPDEs)):
			try:
				print(systemsOfPDEs[k])
				self.string_storer = str(systemsOfPDEs[k])
				self.string_storer = self.string_storer[-100:]
				Clock.schedule_once(self.scheduleOnceTester, 0)
			except:
				raise
		orders_list = [max_order_term(equation, func)[3] for equation in systemsOfPDEs]
		if len(systemsOfPDEs) == 2:
			maximum_order = max(orders_list)
			maximum_order_index = orders_list.index(maximum_order)
			equation_2_new = expand(self.simplificationPDE(systemsOfPDEs[0], systemsOfPDEs[1], func))
			if equation_2_new == 0:
				return [depth, True]
			else:
				if depth >= depth_max: return [depth, False]
				try:
					print([systemsOfPDEs[maximum_order_index], equation_2_new])
					self.string_storer = str([systemsOfPDEs[maximum_order_index], equation_2_new])
					self.string_storer = self.string_storer[-100:]
					Clock.schedule_once(self.scheduleOnceTester, 0)
					return self.ContextFreeDifferentialCommutationCompatibilitySystemPDEs([systemsOfPDEs[maximum_order_index], equation_2_new], func, depth_max, depth + 1)
				except Exception:
					return [depth, False]
		else:
			systemsOfPDEs_F = [None for k in range(len(systemsOfPDEs))]
			for k in range(len(systemsOfPDEs)-1):
				systemsOfPDEs_F[k] = expand(self.simplificationPDE(systemsOfPDEs[k], systemsOfPDEs[k+1], func))
			systemsOfPDEs_F[-1] = expand(self.simplificationPDE(systemsOfPDEs[-1], systemsOfPDEs[0], func))
			flag = False
			for k in range(len(systemsOfPDEs_F)):
				print(systemsOfPDEs_F[k])
				self.string_storer = str(systemsOfPDEs_F[k])
				self.string_storer = self.string_storer[-100:]
				Clock.schedule_once(self.scheduleOnceTester, 0)
				flag = flag or systemsOfPDEs_F[k].equals(0)
			if flag: return [depth, flag]
			else:
				if depth >= depth_max: return [depth, False]
				try:
					return self.ContextFreeDifferentialCommutationCompatibilitySystemPDEs(systemsOfPDEs_F, func, depth_max, depth + 1)
				except Exception:
					return [depth, False]

	def ContextSensitiveDifferentialCommutationCompatibilitySystemPDEs(self, systemsOfPDEs, contextualPDEs, phi, list_of_u_functions):
		print('context sensitive check')
		print('Entered system of PDEs =>', systemsOfPDEs)
		print('Number of entered system of PDEs =>', len(systemsOfPDEs))
		print('ContextualPDEs =>', contextualPDEs)
		print('len(ContextualPDEs) =>', len(contextualPDEs))
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
			print(solution_dictionary_for_all_higher_derivs)
			self.string_storer = str(solution_dictionary_for_all_higher_derivs)
			self.string_storer = self.string_storer[-100:]
			Clock.schedule_once(self.scheduleOnceTester, 0)
			if highest_order_term_of_u_function_k not in list_of_all_highest_derivative_terms_in_uk:
				list_of_all_highest_derivative_terms_in_uk.append(highest_order_term_of_u_function_k)
		for k in range(num_kappa):
			print('k = ',k,' num_kappa = ',num_kappa)
			equation_k_from_the_phi_system = systemsOfPDEs[k]
			list_of_terms = equation_k_from_the_phi_system.args
			list_of_differential_terms_from_phi = []
			for m in range(len(list_of_terms)):
				list_of_sub_terms_in_a_term = list_of_terms[m].args
				for differential_term in list_of_sub_terms_in_a_term:
					if differential_term not in list_of_differential_terms_from_phi:
						list_of_differential_terms_from_phi.append(differential_term)
			if k != num_kappa - 1:
				print('equation_k_from_the_phi_system = ', equation_k_from_the_phi_system)
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
				print('highest_order_term_of_phi_from_equation_k = ',highest_order_term_of_phi_from_equation_k)
				if highest_order_term_of_phi_from_equation_k in solution_dictionary_for_all_higher_derivs:
					while True:
						new_equation = expand(equation_k_from_the_phi_system.subs(highest_order_term_of_phi_from_equation_k,
							solution_dictionary_for_all_higher_derivs[highest_order_term_of_phi_from_equation_k]).doit())
						new_equation, _ = fraction(together(new_equation).doit())
						new_equation = expand(new_equation)
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
						print('new highest_order_term_of_phi_from_equation_k -> ', highest_order_term_of_phi_from_equation_k)
						if highest_order_term_of_phi_from_equation_k in solution_dictionary_for_all_higher_derivs: pass
						else:
							new_equation, _ = fraction(together(new_equation).doit())
							new_equation = expand(new_equation)
							solution_of_highest_order_term_of_phi_from_equation_k = expand(solve(new_equation, highest_order_term_of_phi_from_equation_k)[0])
							solution_of_highest_order_term_of_phi_from_equation_k = expand(solution_of_highest_order_term_of_phi_from_equation_k)
							solution_dictionary_for_all_higher_derivs[highest_order_term_of_phi_from_equation_k] = solution_of_highest_order_term_of_phi_from_equation_k
							print(solution_dictionary_for_all_higher_derivs)
							self.string_storer = str(solution_dictionary_for_all_higher_derivs)
							self.string_storer = self.string_storer[-100:]
							Clock.schedule_once(self.scheduleOnceTester, 0)
							if highest_order_term_of_phi_from_equation_k not in list_of_all_highest_derivative_terms_in_phi:
								list_of_all_highest_derivative_terms_in_phi.append(highest_order_term_of_phi_from_equation_k)
							break
				else:
					print('equation_k_from_the_phi_system = ', equation_k_from_the_phi_system)
					solution_of_highest_order_term_of_phi_from_equation_k = solve(equation_k_from_the_phi_system, highest_order_term_of_phi_from_equation_k)[0]
					solution_of_highest_order_term_of_phi_from_equation_k = expand(solution_of_highest_order_term_of_phi_from_equation_k)
					solution_dictionary_for_all_higher_derivs[highest_order_term_of_phi_from_equation_k] = solution_of_highest_order_term_of_phi_from_equation_k
					print(solution_dictionary_for_all_higher_derivs)
					self.string_storer = str(solution_dictionary_for_all_higher_derivs)
					self.string_storer = self.string_storer[-100:]
					Clock.schedule_once(self.scheduleOnceTester, 0)
					if highest_order_term_of_phi_from_equation_k not in list_of_all_highest_derivative_terms_in_phi:
						list_of_all_highest_derivative_terms_in_phi.append(highest_order_term_of_phi_from_equation_k)
			else:
				print('equation_k_from_the_phi_system = ', equation_k_from_the_phi_system)
				order_of_phi_from_equation_k = ode_order(equation_k_from_the_phi_system, phi)
				for term in list_of_differential_terms_from_phi:
					if ode_order(term, phi) == order_of_phi_from_equation_k:
						highest_order_term_of_phi_from_equation_k = term
						break
				if highest_order_term_of_phi_from_equation_k in solution_dictionary_for_all_higher_derivs:
					while True:
						new_equation = expand(equation_k_from_the_phi_system.subs(highest_order_term_of_phi_from_equation_k,
							solution_dictionary_for_all_higher_derivs[highest_order_term_of_phi_from_equation_k]).doit())
						new_equation, _ = fraction(together(new_equation).doit())
						new_equation = expand(new_equation)
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
							print('equation_k_from_the_phi_system = ', equation_k_from_the_phi_system)
							print('new_equation = ', new_equation)
							print('highest_order_term_of_phi_from_equation_k = ', highest_order_term_of_phi_from_equation_k)
							self.string_storer = str(highest_order_term_of_phi_from_equation_k)
							self.string_storer = self.string_storer[-100:]
							Clock.schedule_once(self.scheduleOnceTester, 0)
							print('new_equation = ', new_equation)
							self.string_storer = str(new_equation)
							self.string_storer = self.string_storer[-100:]
							Clock.schedule_once(self.scheduleOnceTester, 0)
							solution_of_highest_order_term_of_phi_from_equation_k = solve(new_equation, highest_order_term_of_phi_from_equation_k)[0]
							solution_of_highest_order_term_of_phi_from_equation_k = expand(solution_of_highest_order_term_of_phi_from_equation_k)
							solution_dictionary_for_all_higher_derivs[highest_order_term_of_phi_from_equation_k] = solution_of_highest_order_term_of_phi_from_equation_k
							self.string_storer = str(solution_dictionary_for_all_higher_derivs)
							self.string_storer = self.string_storer[-100:]
							Clock.schedule_once(self.scheduleOnceTester, 0)
							if highest_order_term_of_phi_from_equation_k not in list_of_all_highest_derivative_terms_in_phi:
								list_of_all_highest_derivative_terms_in_phi.append(highest_order_term_of_phi_from_equation_k)
							break
				else:
					print('equation_k_from_the_phi_system = ', equation_k_from_the_phi_system)
					solution_of_highest_order_term_of_phi_from_equation_k = solve(equation_k_from_the_phi_system, highest_order_term_of_phi_from_equation_k)[0]
					solution_of_highest_order_term_of_phi_from_equation_k = expand(solution_of_highest_order_term_of_phi_from_equation_k).doit()
					solution_dictionary_for_all_higher_derivs[highest_order_term_of_phi_from_equation_k] = solution_of_highest_order_term_of_phi_from_equation_k
					print('highest_order_term_of_phi_from_equation_k ->', highest_order_term_of_phi_from_equation_k)
					print('solution_dictionary_for_all_higher_derivs[highest_order_term_of_phi_from_equation_k]', solution_dictionary_for_all_higher_derivs[highest_order_term_of_phi_from_equation_k])
					self.string_storer = str(solution_dictionary_for_all_higher_derivs)
					self.string_storer = self.string_storer[-100:]
					Clock.schedule_once(self.scheduleOnceTester, 0)
					if highest_order_term_of_phi_from_equation_k not in list_of_all_highest_derivative_terms_in_phi:
						list_of_all_highest_derivative_terms_in_phi.append(highest_order_term_of_phi_from_equation_k)
		# From here on I have gotten the basic preparatory lists maintained, now I am following the rest of it from the pseudocode of the paper
		S_equations, S_variables = [], []
		for variable, substitution in solution_dictionary_for_all_higher_derivs.items():
			equation = expand(variable - substitution)
			if variable not in S_variables:
				S_variables.append(variable)
				S_equations.append(equation)
		self.string_storer = str(list_of_all_highest_derivative_terms_in_phi)
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		for k in range(len(S_equations)):
			self.string_storer = str(S_equations[k])
			self.string_storer = self.string_storer[-100:]
			Clock.schedule_once(self.scheduleOnceTester, 0)
		for k in range(len(S_variables)):
			self.string_storer = str(S_variables[k])
			self.string_storer = self.string_storer[-100:]
			Clock.schedule_once(self.scheduleOnceTester, 0)
		for k in range(len(S_variables)):
			variable = S_variables[k]
			equation = S_equations[k]
			solution = expand(variable-equation)
			solution_dictionary_for_all_higher_derivs[variable] = solution
			self.string_storer = str(S_equations[k])
			self.string_storer = self.string_storer[-100:]
			Clock.schedule_once(self.scheduleOnceTester, 0)
			self.string_storer = str(S_variables[k])
			self.string_storer = self.string_storer[-100:]
			Clock.schedule_once(self.scheduleOnceTester, 0)

		print('S_variables = ', S_variables)
		print('S_equations = ', S_equations)
		for variable in solution_dictionary_for_all_higher_derivs:
			print('solution_dictionary_for_all_higher_derivs[',variable,'] = ',solution_dictionary_for_all_higher_derivs[variable])
		print('list_of_all_highest_derivative_terms_in_phi = ', list_of_all_highest_derivative_terms_in_phi)

		list_of_A_orders, list_of_B_orders = [], []
		for term in list_of_all_highest_derivative_terms_in_phi:
			list_of_partial_derivatives = term.variables
			print('term = ', term)
			print('list_of_partial_derivatives = ',list_of_partial_derivatives)
			count_of_t, count_of_x = 0, 0
			for variable in list_of_partial_derivatives:
				if variable == t: count_of_t = count_of_t + 1
				elif variable == x: count_of_x = count_of_x + 1
			list_of_A_orders.append(count_of_x)
			list_of_B_orders.append(count_of_t)
		A_max, A_min, B_max, B_min = max(list_of_A_orders), min(list_of_A_orders), max(list_of_B_orders), min(list_of_B_orders)
		print('list_of_A_orders -> ', list_of_A_orders, 'A_max = ',A_max, 'A_min = ',A_min)
		print('list_of_B_orders -> ', list_of_B_orders, 'B_max = ',B_max, 'B_min = ',B_min)
		list_of_C_orders, list_of_D_orders = [], []
		for term in list_of_all_highest_derivative_terms_in_uk:
			list_of_partial_derivatives = term.variables
			count_of_t, count_of_x = 0, 0
			for variable in list_of_partial_derivatives:
				if variable == t: count_of_t = count_of_t + 1
				elif variable == x: count_of_x = count_of_x + 1
			list_of_C_orders.append(count_of_x)
			list_of_D_orders.append(count_of_t)
		print('list_of_C_orders -> ', list_of_C_orders)
		print('list_of_D_orders -> ', list_of_D_orders)
		for v in range(len(list_of_all_highest_derivative_terms_in_uk)):
			Cv = list_of_C_orders[v]
			Dv = list_of_D_orders[v]
			uv = list_of_all_highest_derivative_terms_in_uk[v]
			gv = solution_dictionary_for_all_higher_derivs[uv]
			for i in range(15+max([A_max, B_max])):
				for j in range(15+max([A_max, B_max])):
					equation = expand(diff(uv,(x,i+Cv),(t,j+Dv))-diff(gv,(x,i),(t,j)))
					variable = diff(uv,(x,i+Cv),(t,j+Dv))
					for other_variable in solution_dictionary_for_all_higher_derivs:
						if other_variable != variable:
							equation = expand(equation.subs(other_variable, solution_dictionary_for_all_higher_derivs[other_variable]).doit())
							equation = expand(equation)
					if equation not in S_equations: S_equations.append(equation)
					if variable not in S_variables: S_variables.append(variable)
					if variable not in solution_dictionary_for_all_higher_derivs:
						solution_dictionary_for_all_higher_derivs[variable] = expand(equation-variable).doit()
					for k in range(len(S_equations)):
						self.string_storer = str(S_equations[k])
						self.string_storer = self.string_storer[-100:]
						Clock.schedule_once(self.scheduleOnceTester, 0)
					for k in range(len(S_variables)):
						self.string_storer = str(S_variables[k])
						self.string_storer = self.string_storer[-100:]
						Clock.schedule_once(self.scheduleOnceTester, 0)
		for w in range(len(list_of_A_orders)):
			Aw = list_of_A_orders[w]
			Bw = list_of_B_orders[w]
			phi_term_w = list_of_all_highest_derivative_terms_in_phi[w]
			fw = solution_dictionary_for_all_higher_derivs[phi_term_w]
			for derivsInX in range(14):
				print('derivsInX = ', derivsInX)
				for derivsInT in range(14):
					print('derivsInT = ', derivsInT)
					equation = expand(diff(phi,(x,Aw+derivsInX),(t,Bw+derivsInT))-diff(fw,(x,derivsInX),(t,derivsInT)))
					print('equation = ', equation)
					variable = diff(phi,(x,Aw+derivsInX),(t,Bw+derivsInT))
					print('variable = ', variable)
					for other_variable in solution_dictionary_for_all_higher_derivs:
						print('other_variable = ', other_variable)
						if other_variable != variable:
							equation = expand(equation.subs(other_variable, solution_dictionary_for_all_higher_derivs[other_variable]).doit())
							equation = expand(equation)
					print('equation not in S_equations = ', equation not in S_equations)
					print('variable not in S_variables = ', variable not in S_variables)
					print('variable not in solution_dictionary_for_all_higher_derivs = ', variable not in solution_dictionary_for_all_higher_derivs)
					if equation not in S_equations: S_equations.append(equation)
					if variable not in S_variables: S_variables.append(variable)
					if variable not in solution_dictionary_for_all_higher_derivs:
						solution_dictionary_for_all_higher_derivs[variable] = expand(equation-variable).doit()
					for k in range(len(S_equations)):
						self.string_storer = str(S_equations[k])
						self.string_storer = self.string_storer[-100:]
						Clock.schedule_once(self.scheduleOnceTester, 0)
					for k in range(len(S_variables)):
						self.string_storer = str(S_variables[k])
						self.string_storer = self.string_storer[-100:]
						Clock.schedule_once(self.scheduleOnceTester, 0)
		# Now having set all the equations and the variables I will have to be a little careful
		# setting up the equation solver both in terms of the data structure and also
		# the discontinuation of the old solutions.
		print('S_variables = ', S_variables)
		print('S_equations = ', S_equations)
		print('solution_dictionary_for_all_higher_derivs = ', solution_dictionary_for_all_higher_derivs)
		while True:
			S_new_equations = []
			solution_to_data_structure = list(nonlinsolve(S_equations, S_variables))
			print('solution_to_data_structure = ', solution_to_data_structure)
			try:
				for variable, solution in solution_to_data_structure.items():
					if variable in S_variables:
						equation = expand(variable-solution)
						if equation not in S_new_equations:
							S_new_equations.append(equation)
			except:
				print('S_new_equations = ',S_new_equations)
				print('len(solution_to_data_structure) = ', len(solution_to_data_structure))
				print('len(solution_to_data_structure[0]) = ', len(solution_to_data_structure[0]))
				new_solution_set = solution_to_data_structure[0]
				for k in range(len(S_variables)):
					variable = S_variables[k]
					solution = new_solution_set[k]
					equation = expand(variable-solution)
					if equation not in S_new_equations:
						S_new_equations.append(equation)
			flag = len(list(set(S_new_equations) - set(S_equations))) == 0
			S_equations = S_new_equations
			for k in range(len(S_equations)):
				self.string_storer = str(S_equations[k])
				self.string_storer = self.string_storer[-100:]
				Clock.schedule_once(self.scheduleOnceTester, 0)
			for k in range(len(S_variables)):
				self.string_storer = str(S_variables[k])
				self.string_storer = self.string_storer[-100:]
				Clock.schedule_once(self.scheduleOnceTester, 0)
			print('S_variables = ', S_variables)
			print('S_equations = ', S_equations)
			print('len(S_variables) = ', len(S_variables))
			print('len(S_equations) = ', len(S_equations))
			print(flag)
			if flag: break

		print('len(S_variables) = ', len(S_variables))
		print('len(S_equations) = ', len(S_equations))
		for k in range(len(S_variables)):
			solution_dictionary_for_all_higher_derivs[S_variables[k]] = expand(S_equations[k]-S_variables[k])

		print(solution_dictionary_for_all_higher_derivs)
		print('S_variables ->')
		for k in range(len(S_variables)):
			print('k -> ', k, ' ->', S_variables[k])
		print('len(S_variables) -> ', len(S_variables))
		print('S_equations ->')
		for k in range(len(S_equations)):
			print('k -> ', k, ' ->', S_equations[k])
		print('len(S_equations) -> ', len(S_equations))

		n_0 = len(list_of_all_highest_derivative_terms_in_phi)
		W_set, R_set = [None for k in range(n_0)], [None for k in range(n_0)]
		flag1 = True
		for k in range(n_0-1):
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
				self.string_storer = str(target_root_order)
				self.string_storer = self.string_storer[-100:]
				Clock.schedule_once(self.scheduleOnceTester, 0)
				nth_root_instances = [term for term in W_set[k].find(lambda x: isinstance(x, Pow) and x.exp == 1/target_root_order)]
				has_nth_root = len(nth_root_instances) > 0
				if has_nth_root:
					reference_expression = expand(W_set[k]**target_root_order)
					radical = nth_root_instances[0]
					W_set[k] = resultant(W_set[k], reference_expression, radical)
					break
				target_root_order = target_root_order + 1
			print('W_set['+str(k)+'] = ',W_set[k])
			flag1 = flag1 and W_set[k] == 0
			self.string_storer = str(W_set[k])
			self.string_storer = self.string_storer[-100:]
			Clock.schedule_once(self.scheduleOnceTester, 0)
		A_n_minus_1, A_0 = list_of_A_orders[n_0-1], list_of_A_orders[0]
		B_n_minus_1, B_0 = list_of_B_orders[n_0-1], list_of_B_orders[0]
		phi_term_n_minus_1 = list_of_all_highest_derivative_terms_in_phi[n_0-1]
		f_n_minus_1 = solution_dictionary_for_all_higher_derivs[phi_term_n_minus_1]
		phi_term_0 = list_of_all_highest_derivative_terms_in_phi[0]
		f_0 = solution_dictionary_for_all_higher_derivs[phi_term_0]
		W_set[n_0-1] = expand(diff(f_n_minus_1, (x,A_max-A_n_minus_1), (t, B_max-B_n_minus_1))-diff(f_0, (x,A_max-A_0), (t, B_max-B_0)))
		W_set[n_0-1], _ = fraction(together(W_set[n_0-1]).doit())
		W_set[n_0-1] = expand(W_set[n_0-1])
		for variable, solution in solution_dictionary_for_all_higher_derivs.items():
			W_set[n_0-1] = expand(W_set[n_0-1].subs(variable, solution))
			W_set[n_0-1], _ = fraction(together(W_set[n_0-1]).doit())
			W_set[n_0-1] = expand(W_set[n_0-1])
		target_root_order = 2
		while target_root_order < 20:
			nth_root_instances = [term for term in W_set[n_0-1].find(lambda x: isinstance(x, Pow) and x.exp == 1/target_root_order)]
			has_nth_root = len(nth_root_instances) > 0
			if has_nth_root:
				reference_expression = expand(W_set[n_0-1]**target_root_order)
				radical = nth_root_instances[0]
				W_set[n_0-1] = resultant(W_set[n_0-1], reference_expression, radical)
				break
			target_root_order = target_root_order + 1
		print('W_set['+str(n_0-1)+'] = ',W_set[n_0-1])
		self.string_storer = str(W_set[n_0-1])
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		flag1 = flag1 and W_set[n_0-1] == 0
		flag2 = True
		expression_terms = W_set[-1].as_coefficients_dict()
		list_of_differential_terms_from_one_of_the_last_terms = []
		for term, coefficient in expression_terms.items():
			if isinstance(term, Derivative) and term not in list_of_differential_terms_from_one_of_the_last_terms:
				list_of_differential_terms_from_one_of_the_last_terms.append(term)
			elif isinstance(term, Mul):
				for sub_term in term.args:
					if isinstance(sub_term, Derivative) and sub_term not in list_of_differential_terms_from_one_of_the_last_terms:
						list_of_differential_terms_from_one_of_the_last_terms.append(sub_term)
		phi_term_random = choice(list_of_differential_terms_from_one_of_the_last_terms)
		print('phi_term_random = ',phi_term_random)
		for k in range(n_0-1):
			R_set[k] = expand(resultant(
				W_set[k],W_set[k+1],\
				phi_term_random
				))
			print('R_set['+str(k)+'] = ',R_set[k])
			flag2 = flag2 and R_set[k] == 0
			self.string_storer = str(R_set[k])
			self.string_storer = self.string_storer[-100:]
			Clock.schedule_once(self.scheduleOnceTester, 0)
		R_set[n_0-1] = expand(resultant(
			W_set[n_0-1],W_set[0],\
			phi_term_random
			))
		print('R_set['+str(n_0-1)+'] = ',R_set[n_0-1])
		flag2 = flag2 and R_set[n_0-1] == 0
		self.string_storer = str(R_set[n_0-1])
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		return (flag1 or flag2)

	def BacklundContextFreeAndContextSensitiveCompatibility(self):
		x, t, sigma, b, kappa = symbols('x t sigma b kappa')
		phi = Function('phi')(x, t)
		f = Function('f')
		input_pde_expression = parse_expr(self.root.ids.input_pde.text)
		alpha, _, _ = self.ultimatePainleve(input_pde_expression, x, t)
		beta_new = -alpha
		u_coeff = [Function('u'+str(k))(x, t) for k in range(beta_new+1)]
		u = 0
		for k in range(beta_new+1):
			u = u + u_coeff[k]/phi**(beta_new-k)
		backlund_transform_store = u
		print(u)
		generator_list = [phi, 1/phi]
		for k in range(10):
			for m in range(10):
				if diff(phi,(x,k),(t,m)) not in generator_list:
					generator_list.append(diff(phi,(x,k),(t,m)))
		painleve_expansion = expand(input_pde_expression.subs(f(x,t),u).doit())
		painleve_expansion_poly_in_terms_of_phi = Poly(painleve_expansion, gens = generator_list)
		painleve_expansion_poly_in_terms_of_phi = painleve_expansion_poly_in_terms_of_phi.as_expr()
		mu_min = degree(painleve_expansion_poly_in_terms_of_phi, 1/phi)
		print(mu_min)
		equation_list = [0 for k in range(mu_min+1)]
		u_solution = [0 for k in range(beta_new+1)]
		for k in range(beta_new):
			equation_list[k] = simplify(painleve_expansion_poly_in_terms_of_phi.coeff(1/phi,mu_min-k))
			print('k -> ', k, ' equation_list[k] = ', equation_list[k])
			for m in range(k):
				equation_list[k] = expand(equation_list[k].subs(u_coeff[m], u_solution[m]).doit())
			if equation_list[k] == 0:
				k_max_with_non_zero_coeffs = k
				break
			factors_list_of_equationk = factor_list(equation_list[k])
			try: equation_list[k] = factors_list_of_equationk[-1][-1][0]
			except: pass
			u_solution[k] = solve(equation_list[k],u_coeff[k])[0]
			print('k -> ', k, ' equation_list[k] = ', equation_list[k])
			print('k -> ', k, ' u_solution[k] = ', u_solution[k])
			backlund_transform_store = backlund_transform_store.subs(u_coeff[k], u_solution[k])
			assert expand(equation_list[k].subs(u_coeff[k], u_solution[k]).doit()) == 0
		if k == beta_new-1: k_max_with_non_zero_coeffs = beta_new
		print('k_max_with_non_zero_coeffs = ',k_max_with_non_zero_coeffs)
		self.backlund_transform = backlund_transform_store
		k_lower = k_max_with_non_zero_coeffs if k_max_with_non_zero_coeffs == beta_new else k_max_with_non_zero_coeffs+1
		ListOfEquation = [[None for j in range(mu_min+1-k_lower)] for i in range(k_max_with_non_zero_coeffs)]
		for k in range(k_lower, mu_min+1):
			if k!=mu_min:
				equation_list[k] = simplify(painleve_expansion_poly_in_terms_of_phi.coeff(1/phi,mu_min-k))
			else:
				equation_list[k] = simplify(painleve_expansion_poly_in_terms_of_phi.coeff(phi,0))
			for m in range(k_max_with_non_zero_coeffs):
				ListOfEquation[m][k-k_lower] = expand(equation_list[k].subs(u_coeff[m], u_solution[m]).doit())
				equation_list[k] = ListOfEquation[m][k-k_lower]
				ListOfEquation[m][k-k_lower], _ = fraction(together(ListOfEquation[m][k-k_lower]).doit())
				ListOfEquation[m][k-k_lower] = expand(ListOfEquation[m][k-k_lower])
		print('ListOfEquation = ', ListOfEquation)
		print('ListOfEquation[-1] = ', ListOfEquation[-1])
		self.string_storer = str(ListOfEquation[-1])
		self.string_storer = self.string_storer[-100:]
		Clock.schedule_once(self.scheduleOnceTester, 0)
		TotalSystemOfPDEs = ListOfEquation[-1]
		while True:
			try:
				TotalSystemOfPDEs.remove(0)
			except:
				break
		contextualPDEs = [TotalSystemOfPDEs[-1]]
		systemOfPDEs = TotalSystemOfPDEs
		systemOfPDEs.remove(TotalSystemOfPDEs[-1])
		filteredSystemOfPDEs = systemOfPDEs
		print('filteredSystemOfPDEs =',filteredSystemOfPDEs)
		print('len(filteredSystemOfPDEs) =',len(filteredSystemOfPDEs))
		flag_context_free_compatibility = self.ContextFreeDifferentialCommutationCompatibilitySystemPDEs(filteredSystemOfPDEs+contextualPDEs, phi)
		stored_filteredSystemOfPDEs = [filteredSystemOfPDEs[k] for k in range(len(filteredSystemOfPDEs))]
		if flag_context_free_compatibility[1]:
			self.context_free_compatibility_test = True
			while len(filteredSystemOfPDEs)+len(contextualPDEs) >= 2:
				switch_track = False
				for f1 in filteredSystemOfPDEs:
					print('f1 = ', f1)
					newSystemOfPDEs = filteredSystemOfPDEs
					newSystemOfPDEs.remove(f1)
					print('newSystemOfPDEs+contextualPDEs =',newSystemOfPDEs+contextualPDEs)
					print('len(newSystemOfPDEs)+len(contextualPDEs) =',len(newSystemOfPDEs)+len(contextualPDEs))
					self.string_storer = str(newSystemOfPDEs)
					self.string_storer = self.string_storer[-100:]
					Clock.schedule_once(self.scheduleOnceTester, 0)
					if len(newSystemOfPDEs)+len(contextualPDEs) == 2:
						self.context_free_compatibility_test = True
						self.context_sensitive_compatibility_test = True
						self.compatibility_test = True
						self.associated_conditions = [f1] # newSystemOfPDEs
						return
					else:
						print('>> newSystemOfPDEs+contextualPDEs =',newSystemOfPDEs+contextualPDEs)
						print('>> len(newSystemOfPDEs)+len(contextualPDEs) =',len(newSystemOfPDEs)+len(contextualPDEs))
						flag = self.ContextFreeDifferentialCommutationCompatibilitySystemPDEs(newSystemOfPDEs+contextualPDEs, phi)
						if flag[1]:
							filteredSystemOfPDEs = newSystemOfPDEs
							continue
						else:
							filteredSystemOfPDEs = newSystemOfPDEs
							switch_track = True
							break
				print('switch_track = ',switch_track)
				if switch_track: break
			print('filteredSystemOfPDEs =', filteredSystemOfPDEs)
			print('newSystemOfPDEs = ', newSystemOfPDEs)
			print('stored_filteredSystemOfPDEs =', stored_filteredSystemOfPDEs)
			print('stored_filteredSystemOfPDEs\\filteredSystemOfPDEs = ',list(set(stored_filteredSystemOfPDEs)-set(filteredSystemOfPDEs)))
			print('stored_filteredSystemOfPDEs\\newSystemOfPDEs =', list(set(stored_filteredSystemOfPDEs)-set(newSystemOfPDEs)))
			print('filteredSystemOfPDEs\\newSystemOfPDEs =', list(set(filteredSystemOfPDEs)-set(newSystemOfPDEs)))
			print('contextualPDEs =', contextualPDEs)
			if self.ContextSensitiveDifferentialCommutationCompatibilitySystemPDEs(filteredSystemOfPDEs, contextualPDEs, phi, [u_coeff[beta_new]]):
				self.context_free_compatibility_test = True
				self.context_sensitive_compatibility_test = True
				self.compatibility_test = True
				self.associated_conditions = filteredSystemOfPDEs
				return
			else:
				self.context_free_compatibility_test = True
				self.context_sensitive_compatibility_test = False
				self.compatibility_test = False
				self.associated_conditions = filteredSystemOfPDEs
				return
		elif self.ContextSensitiveDifferentialCommutationCompatibilitySystemPDEs(filteredSystemOfPDEs, contextualPDEs, phi, [u_coeff[beta_new]]):
			self.context_free_compatibility_test = False
			self.context_sensitive_compatibility_test = True
			self.compatibility_test = True
			self.associated_conditions = filteredSystemOfPDEs
			return
		else:
			self.context_free_compatibility_test = True
			self.context_sensitive_compatibility_test = False
			self.compatibility_test = False
			self.associated_conditions = filteredSystemOfPDEs
			return
	def show_dialog(self):
		if not self.dialog:
			self.dialog = MDDialog(
				title = 'Some calculation steps',
				type = 'custom',
				content_cls = Content()
				)
		# Introduce string storers and clock schedulers in appropriate places ....
		thread = threading.Thread(target=self.BacklundContextFreeAndContextSensitiveCompatibility)
		thread.start()
		self.dialog.open()

	def close_handler(self):
		if self.dialog:
			self.dialog.dismiss()
			if self.compatibility_test == 0:
				self.root.ids.integrability_handler.text = 'The equation '+self.root.ids.input_pde.text+' is non-integrable'
			else:
				self.root.ids.integrability_handler.text = 'The equation '+self.root.ids.input_pde.text+' is integrable with the auto-Backlund transform '+str(self.backlund_transform)
				self.root.ids.integrability_handler.text = self.root.ids.integrability_handler.text + ' with the following compatible set of '+str(self.associated_conditions)
		if hasattr(self, 'schedule_interval'):
			Clock.unschedule(self.scheduleOnceTester)
		self.dialog = None

PainleveBacklundCheckApp().run()
