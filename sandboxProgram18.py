from kivy.lang import Builder
from kivymd.app import MDApp
from kivy.uix.boxlayout import BoxLayout
from kivymd.uix.dialog import MDDialog
from kivy.clock import Clock
from sympy import *
from time import sleep
import threading
from collections import Counter

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
		font_name: "Georgia"
		font_size: 30
		id: integrability_handler
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
		self.dialog.content_cls.ids.process_displayer.text = self.string_storer

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
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = str(inputPDE2)+' = 0'
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = 'The maximum order term in equation 1 is gamma_1 = '+str(max_order_term1[2])
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = 'The maximum order term in equation 2 is gamma_2 = '+str(max_order_term2[2])
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = 'A1 = '+str(A1)
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = 'A2 = '+str(A2)
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = 'B1 = '+str(B1)
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = 'B2 = '+str(B2)
		Clock.schedule_once(self.scheduleOnceTester, 0)
		
		allAvailableVariables = func.free_symbols

		varsFromPDE1 = list(max_order_term1[2].variables) if max_order_term1[3] != 0 else []
		varsFromPDE2 = list(max_order_term2[2].variables) if max_order_term2[3] != 0 else []

		self.string_storer = str(varsFromPDE1)
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = str(varsFromPDE2)
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
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = str(varsFromPDE4)
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
		Clock.schedule_once(self.scheduleOnceTester, 0)
		self.string_storer = str(inputPDE4)
		Clock.schedule_once(self.scheduleOnceTester, 0)
		
		expression = inputPDE3 - inputPDE4
		expression, _ = fraction(together(expression.doit()))
		expression = expand(expression)

		self.string_storer = str(expression)
		Clock.schedule_once(self.scheduleOnceTester, 0)
		
		return expression

	def ContextFreeDifferentialCommutationCompatibilitySystemPDEs(self, systemsOfPDEs, func, depth_max = 150, depth = 0):

		for k in range(len(systemsOfPDEs)):
			self.string_storer = str(systemsOfPDEs[k])
			Clock.schedule_once(self.scheduleOnceTester, 0)
			try:
				self.string_storer = str(systemsOfPDEs[k])
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
				self.string_storer = 'Number of terms in '+str(k)+'th PDE >>'+str(len(Add.make_args(systemsOfPDEs_F[k])))
				Clock.schedule_once(self.scheduleOnceTester, 0)
				self.string_storer = str(systemsOfPDEs_F[k])
				Clock.schedule_once(self.scheduleOnceTester, 0)
				self.string_storer = str(systemsOfPDEs_F[k])
				Clock.schedule_once(self.scheduleOnceTester, 0)
				if len(Add.make_args(systemsOfPDEs_F[k])) > 250: return [depth, False]
				flag = flag or systemsOfPDEs_F[k].equals(0)
			if flag: return [depth, flag]
			else:
				if depth >= depth_max: return [depth, False]
				try:
					return self.ContextFreeDifferentialCommutationCompatibilitySystemPDEs(systemsOfPDEs_F, func, depth_max, depth + 1)
				except Exception:
					return [depth, False]

	def calculate_compatibility(self):
		x, t, sigma, b, kappa = symbols('x t sigma b kappa')
		phi = Function('phi')(x, t)
		f = Function('f')

		if self.root.ids.input_pde.text == 'diff(f(x,t),t)+f(x,t)*diff(f(x,t),x)+sigma*diff(f(x,t),(x,3))':
			input_pde_expression = parse_expr(self.root.ids.input_pde.text)
			result_package = self.ultimatePainleve(input_pde_expression, x, t)
			print(result_package[0])
			print(result_package[1])
			self.backlund_transform = result_package[1]
			Number_Of_Sets = len(result_package[2])
			alpha = result_package[0]
			u_max = Function('U'+str(-alpha))(x, t)

			for k in range(Number_Of_Sets):
				self.string_storer = str(result_package[2][k])
				Clock.schedule_once(self.scheduleOnceTester, 0)
				for m in range(len(result_package[2][k])):
					factors_list_of_equation = factor_list(result_package[2][k][m])
					result_package[2][k][m] = factors_list_of_equation[-1][-1][0]
					result_package[2][k][m] = expand(result_package[2][k][m])

			flag_switch = False
			flag_on_phi = False
			flag_on_umax = False

			k = -1
			context_free_compatible_store = self.ContextFreeDifferentialCommutationCompatibilitySystemPDEs(result_package[2][k], phi)
			if context_free_compatible_store[1]:
				flag_switch = True
				flag_on_phi = True
			else:
				if ode_order(parse_expr(self.root.ids.input_pde.text), f) > 2:
					context_free_compatible_store = self.ContextFreeDifferentialCommutationCompatibilitySystemPDEs(result_package[2][k], u_max)
					if context_free_compatible_store[1]:
						flag_switch = True
						flag_on_umax = True
			
			self.context_free_compatibility_test = flag_switch
			if not(flag_switch):
				self.string_storer = 'The system of PDEs is not context-free compatible'
				Clock.schedule_once(self.scheduleOnceTester, 0)
			else:
				self.string_storer = 'The system of PDEs is context-free compatible as 0 belongs to g_{ '+str(context_free_compatible_store[0]+1)+'}'
				Clock.schedule_once(self.scheduleOnceTester, 0)
				if flag_on_phi:
					self.string_storer = 'It is based on: '+str(phi)
					Clock.schedule_once(self.scheduleOnceTester, 0)
				if flag_on_umax:
					self.string_storer = 'It is based on: '+str(u_max)
					Clock.schedule_once(self.scheduleOnceTester, 0)

			u2 = Function('u2')(x, t)
			eqn1 = diff(phi,x)*diff(phi,t)+u2*diff(phi,x)**2+4*sigma*diff(phi,x)*diff(phi,(x,3))-3*sigma*diff(phi,(x,2))**2
			eqn2 = diff(phi,x,t)+u2*diff(phi,x,x)+sigma*diff(phi,x,x,x,x)
			self.associated_conditions = [eqn1]
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
				dictionary_of_solutions[diff(phi,(x,k))] = equation_new

			for k in range(1, max_number_of_derivatives):
				if k == 1:
					equation_new = diff(dictionary_of_solutions[diff(u2,t)],x)
				else:
					equation_new = diff(dictionary_of_solutions[diff(u2,(x,k-1),t)],x)
				equation_new = expand(equation_new)
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

			for term in dictionary_of_solutions:
				self.string_storer = str(dictionary_of_solutions[term])
				Clock.schedule_once(self.scheduleOnceTester, 0)

			equation1 = diff(dictionary_of_solutions[diff(phi,x,x,x,t)],x)
			equation1 = expand(equation1)
			equation1 = equation1.subs(diff(phi,(x,4)),dictionary_of_solutions[diff(phi,(x,4))])
			equation1 = expand(equation1)

			equation2 = diff(dictionary_of_solutions[diff(phi,x,x,x,x)],t)
			recursion_track = 0
			while True:
				recursion_track = recursion_track + 1
				for term in dictionary_of_solutions:
					equation2 = equation2.subs(term,dictionary_of_solutions[term])
					equation2 = together(equation2).doit()
					equation2 = expand(equation2)
				numerator, denominator = fraction(together(equation2).doit())
				numerator = expand(numerator)
				order_param = ode_order(numerator, phi)
				if order_param < 4: break

			compatibility_test_expression = equation2 - equation1
			compatibility_test_expression = expand(compatibility_test_expression)
			self.context_sensitive_compatibility_test = compatibility_test_expression == 0
			self.compatibility_test = flag_switch and compatibility_test_expression == 0

		elif self.root.ids.input_pde.text == 'diff(f(x,t),t)+f(x,t)*diff(f(x,t),x)-sigma*diff(f(x,t),(x,2))':
			result_package = self.ultimatePainleve(parse_expr(self.root.ids.input_pde.text), x, t)
			Number_Of_Sets = len(result_package[2])
			alpha = result_package[0]
			self.backlund_transform = result_package[1]
			u_max = Function('U'+str(-alpha))(x, t)

			for k in range(Number_Of_Sets):
				self.string_storer = str(result_package[2][k])
				Clock.schedule_once(self.scheduleOnceTester, 0)
				for m in range(len(result_package[2][k])):
					factors_list_of_equation = factor_list(result_package[2][k][m])
					result_package[2][k][m] = factors_list_of_equation[-1][-1][0]
					result_package[2][k][m] = expand(result_package[2][k][m])

			flag_switch = False
			flag_on_phi = False
			flag_on_umax = False

			k = -1
			context_free_compatible_store = self.ContextFreeDifferentialCommutationCompatibilitySystemPDEs(result_package[2][k], phi)
			self.context_free_compatibility_test = context_free_compatible_store[1]
			if context_free_compatible_store[1]:
				flag_switch = True
				flag_on_phi = True
			else:
				if ode_order(parse_expr(self.root.ids.input_pde.text), f) > 2:
					context_free_compatible_store = self.ContextFreeDifferentialCommutationCompatibilitySystemPDEs(result_package[2][k], u_max)
					if context_free_compatible_store[1]:
						flag_switch = True
						flag_on_umax = True

			self.context_free_compatibility_test = flag_switch
			if not(flag_switch):
				self.string_storer = 'The system of PDEs is not context-free compatible'
				Clock.schedule_once(self.scheduleOnceTester, 0)
			else:
				self.string_storer = 'The system of PDEs is context-free compatible as 0 belongs to g_{ '+str(context_free_compatible_store[0]+1)+'}'
				Clock.schedule_once(self.scheduleOnceTester, 0)
				if flag_on_phi:
					self.string_storer = 'It is based on: '+str(phi)
					Clock.schedule_once(self.scheduleOnceTester, 0)
				if flag_on_umax:
					self.string_storer = 'It is based on: '+str(u_max)
					Clock.schedule_once(self.scheduleOnceTester, 0)

			self.associated_conditions = [result_package[2][k][0]]
			if len([phi, u_max]) == len(list(set(result_package[2][k])-set([result_package[2][k][1]]))):
				if flag_switch:
					self.context_sensitive_compatibility_test = 1
					self.compatibility_test = 1
					self.string_storer = 'The system of PDEs ->'+str(result_package[2][k])+'is perfectly determinate and thus '+self.root.ids.input_pde.text+' is integrable.'
					Clock.schedule_once(self.scheduleOnceTester, 0)
			else:
				if ContextSensitiveDifferentialCommutationCompatibilitySystemPDEs([result_package[2][k][0]], [result_package[2][k][1]], phi, [u_max]):
					self.context_sensitive_compatibility_test = 1
					self.compatibility_test = 1
					self.string_storer = 'The equation '+self.root.ids.input_pde.text+' is integrable'
					Clock.schedule_once(self.scheduleOnceTester, 0)
				else:
					self.context_sensitive_compatibility_test = 0
					self.compatibility_test = 0
					self.string_storer = 'The equation '+self.root.ids.input_pde.text+' is non-integrable'
					Clock.schedule_once(self.scheduleOnceTester, 0)
		
		elif self.root.ids.input_pde.text == '2*f(x,t)*diff(f(x,t),x,t)-2*diff(f(x,t),x)*diff(f(x,t),t)-f(x,t)**3+f(x,t)':
			f = Function('f')
			TransformedSineGordonEquation = parse_expr('2*f(x,t)*diff(f(x,t),x,t)-2*diff(f(x,t),x)*diff(f(x,t),t)-f(x,t)**3+f(x,t)')
			u0 = Function('u0')(x, t)
			u1 = Function('u1')(x, t)
			u2 = Function('u2')(x, t)
			u = u0/phi**2 + u1/phi + u2
			self.backlund_transform = u0/phi**2 + u1/phi + u2
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
			self.string_storer = str(u0_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
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
			self.string_storer = str(u1_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			Equation2 = expand(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 4))
			Equation2 = expand(Equation2.subs(u0, u0_solution).doit())
			Equation2Storage = Equation2
			Equation2 = expand(Equation2.subs(u1, u1_solution).doit())
			Equation2, _ = fraction(together(Equation2).doit())
			Equation2 = expand(Equation2)
			factors_list_of_equation2 = factor_list(Equation2)
			self.string_storer = str(Equation2Storage)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			Equation3 = expand(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 3))
			Equation3 = expand(Equation3.subs(u0, u0_solution).doit())
			Equation3 = expand(Equation3.subs(u1, u1_solution).doit())
			Equation3_compatibility_test = Equation3
			self.string_storer = str(Equation3_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			Equation4 = expand(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 2))
			Equation4 = expand(Equation4.subs(u0, u0_solution).doit())
			Equation4 = expand(Equation4.subs(u1, u1_solution).doit())
			Equation4_compatibility_test = Equation4
			self.string_storer = str(Equation4_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			Equation5 = expand(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 1))
			Equation5 = expand(Equation5.subs(u0, u0_solution).doit())
			Equation5 = expand(Equation5.subs(u1, u1_solution).doit())
			Equation5_compatibility_test = Equation5
			self.string_storer = str(Equation5_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			Equation6 = expand(TransformedSineGordonPainleveExpansion_polyInTermsOfPhi.coeff(phi, 0))
			Equation6 = expand(Equation6.subs(u0, u0_solution).doit())
			Equation6 = expand(Equation6.subs(u1, u1_solution).doit())
			Equation6_compatibility_test = Equation6
			self.string_storer = str(Equation6_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			self.context_free_compatibility_test = self.ContextFreeDifferentialCommutationCompatibilitySystemPDEs([Equation3, Equation4, Equation5, Equation6],phi)

			phi_ttxx_solution = expand(solve(Equation3_compatibility_test, diff(phi,t,t,x,x))[0])
			self.string_storer = str(phi_ttxx_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			u2_xt_solution = expand(solve(Equation6_compatibility_test, diff(u2,t,x))[0])
			self.string_storer = str(u2_xt_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			Equation4_compatibility_test = expand(Equation4_compatibility_test.subs(
				diff(phi,t,t,x,x), phi_ttxx_solution
				).doit())
			Equation4_compatibility_test, _ = fraction(together(Equation4_compatibility_test).doit())
			Equation4_compatibility_test = expand(Equation4_compatibility_test)
			self.string_storer = str(Equation4_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			Equation5_compatibility_test = expand(Equation5_compatibility_test.subs(
				diff(phi,t,t,x,x), phi_ttxx_solution
				).doit())
			Equation5_compatibility_test, _ = fraction(together(Equation5_compatibility_test).doit())
			Equation5_compatibility_test = expand(Equation5_compatibility_test)
			self.string_storer = str(Equation5_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			Equation4_compatibility_test = expand(Equation4_compatibility_test.subs(
				diff(u2,x,t), u2_xt_solution
				).doit())
			Equation4_compatibility_test, _ = fraction(together(Equation4_compatibility_test).doit())
			Equation4_compatibility_test = expand(Equation4_compatibility_test)
			self.string_storer = str(Equation4_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			Equation5_compatibility_test = expand(Equation5_compatibility_test.subs(
				diff(u2,x,t), u2_xt_solution
				).doit())
			Equation5_compatibility_test, _ = fraction(together(Equation5_compatibility_test).doit())
			Equation5_compatibility_test = expand(Equation5_compatibility_test)
			self.string_storer = str(Equation5_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			u2_xt_solution = expand(solve(Equation6_compatibility_test, diff(u2,t,x))[0])
			self.string_storer = str(u2_xt_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			Equation7_compatibility_test = expand(resultant(
				Equation4_compatibility_test,
				Equation5_compatibility_test,
				diff(phi, t, t, x)))
			self.string_storer = str(Equation7_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			Equation8_compatibility_test = expand(resultant(
				Equation4_compatibility_test,
				Equation5_compatibility_test,
				diff(phi, t, x, x)))
			self.string_storer = str(Equation8_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			phi_txx_solution = solve(Equation7_compatibility_test, diff(phi, t, x, x))[0]
			self.string_storer = str(phi_txx_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			phi_ttx_solution = solve(expand(Equation4_compatibility_test.subs(diff(phi,t,x,x),phi_txx_solution).doit()),diff(phi,t,t,x))[-1]
			self.string_storer = str(phi_txx_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			phi_ttxx_solution = expand(
				phi_ttxx_solution.subs(
					diff(phi,t,x,x),phi_txx_solution
					).subs(
					diff(phi,t,t,x),phi_ttx_solution
					).doit())

			self.string_storer = str(phi_ttxx_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			dictionary_of_solutions = {}
			dictionary_of_solutions[diff(phi,x,x,t)] = phi_txx_solution
			dictionary_of_solutions[diff(phi,x,t,t)] = phi_ttx_solution
			dictionary_of_solutions[diff(u2,x,t)] = u2_xt_solution

			phi_ttxx_solution_1 = expand(diff(dictionary_of_solutions[diff(phi,x,x,t)],t).doit())
			for term in dictionary_of_solutions:
				phi_ttxx_solution_1 = phi_ttxx_solution_1.subs(term, dictionary_of_solutions[term]).doit()
				phi_ttxx_solution_1 = expand(phi_ttxx_solution_1).doit()
			self.string_storer = str(phi_ttxx_solution_1)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			phi_ttxx_solution_2 = expand(diff(dictionary_of_solutions[diff(phi,x,t,t)],x).doit())
			for term in dictionary_of_solutions:
				phi_ttxx_solution_2 = phi_ttxx_solution_2.subs(term, dictionary_of_solutions[term]).doit()
				phi_ttxx_solution_2 = expand(phi_ttxx_solution_2).doit()
			self.string_storer = str(phi_ttxx_solution_2)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			compatibility_check_1 = phi_ttxx_solution_1 - phi_ttxx_solution_2
			compatibility_check_1, _ = fraction(together(compatibility_check_1).doit())
			compatibility_check_1 = expand(compatibility_check_1).doit()
			self.string_storer = str(compatibility_check_1)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			compatibility_check_2 = phi_ttxx_solution_1 - phi_ttxx_solution
			compatibility_check_2, _ = fraction(together(compatibility_check_2).doit())
			compatibility_check_2 = expand(compatibility_check_2).doit()
			self.string_storer = str(compatibility_check_2)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			

			compatibility_check_3 = phi_ttxx_solution_2 - phi_ttxx_solution
			compatibility_check_3, _ = fraction(together(compatibility_check_3).doit())
			compatibility_check_3 = expand(compatibility_check_3).doit()
			self.string_storer = str(compatibility_check_3)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			compatibility_check_1_new = diff(phi,x,t)**2-u2*diff(phi,x)*diff(phi,t) - (solve(compatibility_check_1,sqrt(diff(phi,x,t)**2-u2*diff(phi,x)*diff(phi,t)))[0])**2
			compatibility_check_1_new, _ = fraction(together(compatibility_check_1_new).doit())
			compatibility_check_1_new = expand(compatibility_check_1_new)

			compatibility_check_2_new = diff(phi,x,t)**2-u2*diff(phi,x)*diff(phi,t) - (solve(compatibility_check_2,sqrt(diff(phi,x,t)**2-u2*diff(phi,x)*diff(phi,t)))[0])**2
			compatibility_check_2_new, _ = fraction(together(compatibility_check_2_new).doit())
			compatibility_check_2_new = expand(compatibility_check_2_new)

			compatibility_check_3_new = diff(phi,x,t)**2-u2*diff(phi,x)*diff(phi,t) - (solve(compatibility_check_3,sqrt(diff(phi,x,t)**2-u2*diff(phi,x)*diff(phi,t)))[0])**2
			compatibility_check_3_new, _ = fraction(together(compatibility_check_3_new).doit())
			compatibility_check_3_new = expand(compatibility_check_3_new)

			compatibility_check_wrap_1 = expand(compatibility_check_1_new + compatibility_check_2_new**2)
			self.string_storer = str(compatibility_check_wrap_1)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			compatibility_check_wrap_2 = expand(compatibility_check_2_new**5 - compatibility_check_3_new**3)
			self.string_storer = str(compatibility_check_wrap_2)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			compatibility_check_wrap_3 = expand(compatibility_check_1_new**5 + compatibility_check_3_new**6)
			self.string_storer = str(compatibility_check_wrap_3)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			self.associated_conditions = [Equation3_compatibility_test]

			self.compatibility_test = int((compatibility_check_wrap_1 == 0) and (compatibility_check_wrap_2 == 0) \
				and (compatibility_check_wrap_3 == 0))

			if self.compatibility_test:
				self.context_sensitive_compatibility_test = 1
				self.compatibility_test = 1
				self.string_storer = 'The equation '+self.root.ids.input_pde.text+' is integrable'
				Clock.schedule_once(self.scheduleOnceTester, 0)
			else:
				self.context_sensitive_compatibility_test = 0
				self.compatibility_test = 0
				self.string_storer = 'The equation '+self.root.ids.input_pde.text+' is non-integrable'
				Clock.schedule_once(self.scheduleOnceTester, 0)

			
		elif self.root.ids.input_pde.text == 'diff(f(x,t),t,t)-diff(f(x,t),x,x)+2*f(x,t)-2*f(x,t)**3':
			# phi^4 theory
			f = Function('f')
			u0 = Function('u0')(x, t)
			u1 = Function('u1')(x, t)

			u = u0/phi + u1
			self.backlund_transform = u0/phi + u1
			Phi4TheoryEquation = parse_expr('diff(f(x,t),t,t)-diff(f(x,t),x,x)+2*f(x,t)-2*f(x,t)**3')
			Phi4TheoryPainleveExpansion = expand(Phi4TheoryEquation.subs(f(x, t), u).doit())

			generator_list = [
				phi, 1/phi,
				diff(phi, x), diff(phi, t),
				diff(phi, x, x), diff(phi, x, t), diff(phi, t, t),
				diff(phi, x, x, x), diff(phi, x, x, t), diff(phi, x, t, t), diff(phi, t, t, t),
				diff(phi, x, x, x, x), diff(phi, x, x, x, t), diff(phi, x, x, t, t), diff(phi, x, t, t, t), diff(phi, t, t, t, t),
				diff(phi, x, x, x, x, x), diff(phi, x, x, x, x, t), diff(phi, x, x, x, t, t), diff(phi, x, x, t, t, t), diff(phi, x, t, t, t, t), diff(phi, t, t, t, t, t)
				]

			Phi4TheoryPainleveExpansion_polyInTermsOfPhi = Poly(Phi4TheoryPainleveExpansion, gens = generator_list)
			Phi4TheoryPainleveExpansion_polyInTermsOfPhi = Phi4TheoryPainleveExpansion_polyInTermsOfPhi.as_expr()

			Equation0 = simplify(Phi4TheoryPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 3))
			factors_list_of_equation0 = factor_list(Equation0)
			Equation0 = factors_list_of_equation0[-1][-1][0]
			u0_solution = solve(Equation0, u0)[0]
			Equation0_residual = Equation0.subs(u0, u0_solution)
			self.string_storer = str(Equation0_residual)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			Equation1 = expand(Phi4TheoryPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 2))
			factors_list_of_equation1 = factor_list(Equation1)
			Equation1 = factors_list_of_equation1[-1][-1][0]
			Equation1Storage = Equation1
			self.string_storer = str(Equation1Storage)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			Equation2 = expand(Phi4TheoryPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 1))
			Equation2Storage = Equation2
			Equation2 = expand(Equation2)
			factors_list_of_equation2 = factor_list(Equation2)
			self.string_storer = str(Equation2Storage)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			Equation3 = expand(Phi4TheoryPainleveExpansion_polyInTermsOfPhi.coeff(phi, 0))
			Equation3_compatibility_test = Equation3
			self.string_storer = str(Equation3_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			context_free_compatible_store = self.ContextFreeDifferentialCommutationCompatibilitySystemPDEs([
												Equation0,
												Equation1Storage,
												Equation2Storage,
												Equation3_compatibility_test
												], phi)
			self.string_storer = str(context_free_compatible_store)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			self.context_free_compatibility_test = context_free_compatible_store[0]

			phi_x_solution = solve(Equation0, diff(phi,x))[-1]
			self.string_storer = str(phi_x_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			phi_xx_solution = diff(phi_x_solution, x)
			phi_xx_solution = expand(phi_xx_solution.doit())

			phi_xt_solution = diff(phi_x_solution, t)
			phi_xt_solution = expand(phi_xt_solution.doit())

			phi_tt_solution = solve(Equation1Storage, diff(phi,t,t))[-1]
			phi_tt_solution = expand(phi_tt_solution.subs(diff(phi,x),phi_x_solution).subs(diff(phi,x,x),phi_xx_solution).doit())

			solutions = solve([expand(phi_xt_solution-diff(phi,x,t)),expand(phi_tt_solution-diff(phi,t,t))],[diff(phi,x,t),diff(phi,t,t)])
			phi_xt_solution = solutions[diff(phi,x,t)]
			phi_tt_solution = solutions[diff(phi,t,t)]

			self.string_storer = str(phi_xt_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			self.string_storer = str(phi_tt_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			phi_xx_solution = expand(phi_xx_solution.subs(diff(phi, x,t),phi_xt_solution).doit())
			self.string_storer = str(phi_xx_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			u0_tt_solution = solve(Equation2Storage, diff(u0,t,t))[-1]
			self.string_storer = str(u0_tt_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			u1_tt_solution = solve(Equation3_compatibility_test, diff(u1,t,t))[-1]
			self.string_storer = str(u1_tt_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			dictionary_of_solutions = {}
			dictionary_of_solutions[diff(phi,x)] = phi_x_solution
			dictionary_of_solutions[diff(phi,x,t)] = phi_xt_solution
			dictionary_of_solutions[diff(phi,x,x)] = phi_xx_solution
			dictionary_of_solutions[diff(phi,t,t)] = phi_tt_solution
			dictionary_of_solutions[diff(u0,t,t)] = u0_tt_solution
			dictionary_of_solutions[diff(u1,t,t)] = u1_tt_solution

			phi_xtt_solution_1 = expand(diff(dictionary_of_solutions[diff(phi,x,t)],t).doit())
			for term in dictionary_of_solutions:
				phi_xtt_solution_1 = phi_xtt_solution_1.subs(term, dictionary_of_solutions[term]).doit()
				phi_xtt_solution_1 = expand(phi_xtt_solution_1)
			self.string_storer = str(phi_xtt_solution_1)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			phi_xtt_solution_2 = expand(diff(dictionary_of_solutions[diff(phi,t,t)],x).doit())
			for term in dictionary_of_solutions:
				phi_xtt_solution_2 = phi_xtt_solution_2.subs(term, dictionary_of_solutions[term]).doit()
				phi_xtt_solution_2 = expand(phi_xtt_solution_2)
			self.string_storer = str(phi_xtt_solution_2)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			compatibility_check_1 = phi_xtt_solution_1 - phi_xtt_solution_2
			compatibility_check_1, _ = fraction(together(compatibility_check_1).doit())
			compatibility_check_1 = expand(compatibility_check_1.doit())
			self.string_storer = str(compatibility_check_1)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			phi_xxt_solution_1 = expand(diff(dictionary_of_solutions[diff(phi,x,x)],t).doit())
			for term in dictionary_of_solutions:
				phi_xxt_solution_1 = phi_xtt_solution_1.subs(term, dictionary_of_solutions[term]).doit()
				phi_xxt_solution_1 = expand(phi_xtt_solution_1)
			self.string_storer = str(phi_xtt_solution_1)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			phi_xxt_solution_2 = expand(diff(dictionary_of_solutions[diff(phi,x,t)],x).doit())
			for term in dictionary_of_solutions:
				phi_xxt_solution_2 = phi_xxt_solution_2.subs(term, dictionary_of_solutions[term]).doit()
				phi_xxt_solution_2 = expand(phi_xxt_solution_2)
			self.string_storer = str(phi_xtt_solution_2)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			compatibility_check_2 = phi_xxt_solution_1 - phi_xxt_solution_2
			compatibility_check_2, _ = fraction(together(compatibility_check_2).doit())
			compatibility_check_2 = expand(compatibility_check_2.doit())
			self.string_storer = str(compatibility_check_2)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			A1 = compatibility_check_1.coeff(sqrt(diff(phi, t)**2-u0**2))
			B1 = compatibility_check_1 - sqrt(diff(phi, t)**2-u0**2)*compatibility_check_1.coeff(sqrt(diff(phi, t)**2-u0**2))
			A2 = compatibility_check_2.coeff(sqrt(diff(phi, t)**2-u0**2))
			B2 = compatibility_check_2 - sqrt(diff(phi, t)**2-u0**2)*compatibility_check_2.coeff(sqrt(diff(phi, t)**2-u0**2))

			self.string_storer = str(A1)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			self.string_storer = str(B1)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			self.string_storer = str(A2)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			self.string_storer = str(B2)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			d1 = diff(phi, t)**2-u0**2-B1**2/A1**2
			d1, _ = fraction(together(d1).doit())
			d1 = expand(d1)
			self.string_storer = str(d1)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			d2 = diff(phi, t)**2-u0**2-B2**2/A2**2
			d2, _ = fraction(together(d2).doit())
			d2 = expand(d2)
			self.string_storer = str(d2)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			generator_list = [
				phi, 1/phi, u0, u1,
				diff(phi, x), diff(phi, t),
				diff(phi, x, x), diff(phi, x, t), diff(phi, t, t),
				diff(phi, x, x, x), diff(phi, x, x, t), diff(phi, x, t, t), diff(phi, t, t, t),
				diff(phi, x, x, x, x), diff(phi, x, x, x, t), diff(phi, x, x, t, t), diff(phi, x, t, t, t), diff(phi, t, t, t, t),
				diff(phi, x, x, x, x, x), diff(phi, x, x, x, x, t), diff(phi, x, x, x, t, t), diff(phi, x, x, t, t, t), diff(phi, x, t, t, t, t), diff(phi, t, t, t, t, t),
				diff(u0, x), diff(u0, t), diff(u0, x, x), diff(u0, x, t), diff(u0, t, t),
				diff(u1, x), diff(u1, t), diff(u1, x, x), diff(u1, x, t), diff(u1, t, t)
				]

			d1 = Poly(d1, gens = generator_list)
			d2 = Poly(d2, gens = generator_list)
			compatibility_check = resultant(d1, d2)
			# self.associated_conditions = [Equation3_compatibility_test]
			self.string_storer = str(compatibility_check)
			self.compatibility_test = int(compatibility_check == 0)
			Clock.schedule_once(self.scheduleOnceTester, 0)

		elif self.root.ids.input_pde.text == 'diff(f(x,t),t)+diff(f(x,t),x)+f(x,t)*diff(f(x,t),x)-diff(f(x,t),x,x,t)':
			f = Function('f')
			u0 = Function('u0')(x, t)
			u1 = Function('u1')(x, t)
			u2 = Function('u2')(x, t)

			u = u0/phi**2 + u1/phi + u2
			self.backlund_transform = u0/phi**2 + u1/phi + u2
			BBMEquation = parse_expr(self.root.ids.input_pde.text)
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
			self.string_storer = str(u0_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			Equation1 = expand(BBMPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 4))
			Equation1 = expand(Equation1.subs(u0, u0_solution).doit())
			factors_list_of_equation1 = factor_list(Equation1)
			Equation1 = factors_list_of_equation1[-1][-1][0]
			self.string_storer = str(Equation1)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			Equation2 = expand(BBMPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 3))
			Equation2 = expand(Equation2.subs(u0, u0_solution).doit())
			Equation2, _ = fraction(together(Equation2).doit())
			Equation2 = expand(Equation2)
			factors_list_of_equation2 = factor_list(Equation2)
			Equation2 = factors_list_of_equation2[-1][-1][0]
			self.string_storer = str(Equation2)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			Equation3 = expand(BBMPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 2))
			Equation3 = expand(Equation3.subs(u0, u0_solution).doit())
			Equation3, _ = fraction(together(Equation3).doit())
			Equation3 = expand(Equation3)
			factors_list_of_equation3 = factor_list(Equation3)
			Equation3 = factors_list_of_equation3[-1][-1][0]
			self.string_storer = str(Equation3)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			Equation4 = expand(BBMPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 1))
			Equation4, _ = fraction(together(Equation4).doit())
			Equation4 = expand(Equation4)
			factors_list_of_equation4 = factor_list(Equation4)
			Equation4 = factors_list_of_equation4[-1][-1][0]
			self.string_storer = str(Equation4)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			Equation5 = expand(BBMPainleveExpansion_polyInTermsOfPhi.coeff(phi, 0))
			Equation5, _ = fraction(together(Equation5).doit())
			Equation5 = expand(Equation5)
			self.string_storer = str(Equation5)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			self.context_free_compatibility_test = self.ContextFreeDifferentialCommutationCompatibilitySystemPDEs([Equation1,Equation2,Equation3,Equation4,Equation5],phi)

			eqn51 = Equation1
			eqn54 = Equation3
			self.string_storer = str(eqn51)+'= 0'
			Clock.schedule_once(self.scheduleOnceTester, 0)
			self.string_storer = str(eqn54)+'= 0'
			Clock.schedule_once(self.scheduleOnceTester, 0)
			solution = solve([eqn51,eqn54],[diff(phi,t,t),diff(phi,t,x,x,x)])

			self.string_storer = str(solution[diff(phi,t,t)])
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			self.string_storer = str(solution[diff(phi,t,x,x,x)])
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			eqn55 = Equation4
			eqn56 = Equation5
			
			self.string_storer = str(eqn55)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			self.string_storer = str(eqn56)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			u1_txx_solution = solve(eqn55,diff(u1,t,x,x))[0]
			u2_txx_solution = solve(eqn56,diff(u2,t,x,x))[0]

			self.string_storer = str(u1_txx_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			self.string_storer = str(u2_txx_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			dictionary_of_solutions = {}
			dictionary_of_solutions[diff(phi,t,x,x,x)] = solution[diff(phi,t,x,x,x)]
			dictionary_of_solutions[diff(phi,t,t)] = solution[diff(phi,t,t)]
			dictionary_of_solutions[diff(u1,t,x,x)] = u1_txx_solution
			dictionary_of_solutions[diff(u2,t,x,x)] = u2_txx_solution

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

			self.string_storer = str(solution_list[diff(phi,t,t,x,x)])
			Clock.schedule_once(self.scheduleOnceTester, 0)
			self.string_storer = str(solution_list[diff(phi,t,x,x,x)])
			Clock.schedule_once(self.scheduleOnceTester, 0)

			for term in dictionary_of_solutions:
				self.string_storer = str(dictionary_of_solutions[term])
				Clock.schedule_once(self.scheduleOnceTester, 0)
			
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

			for term in dictionary_of_solutions:
				self.string_storer = str(dictionary_of_solutions[term])
				Clock.schedule_once(self.scheduleOnceTester, 0)

			self.string_storer = str(equation1)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			equation2 = diff(dictionary_of_solutions[diff(phi,t,x,x,x)],t)
			equation2 = expand(equation2).doit()

			for term in dictionary_of_solutions:
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

			for term in dictionary_of_solutions:
				self.string_storer = str(dictionary_of_solutions[term])
				Clock.schedule_once(self.scheduleOnceTester, 0)

			self.string_storer = str(equation2)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			compatibility_check = equation2 - equation1
			compatibility_check = expand(compatibility_check)
			self.string_storer = str(compatibility_check)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			self.compatibility_test = int(compatibility_check == 0)
			self.context_sensitive_compatibility_test = int(compatibility_check == 0)

		elif self.root.ids.input_pde.text == 'diff(f(x,t),t)+f(x,t)*diff(f(x,t),x)-sigma*diff(f(x,t),(x,2))+kappa*diff(f(x,t),(x,3))':
			f = Function('f')
			u0 = Function('u0')(x, t)
			u1 = Function('u1')(x, t)
			u2 = Function('u2')(x, t)
			u = u0/phi**2 + u1/phi + u2
			KdVBurgersEquation = parse_expr(self.root.ids.input_pde.text)
			self.backlund_transform = u0/phi**2 + u1/phi + u2
			KdVBurgersPainleveExpansion = expand(KdVBurgersEquation.subs(f(x, t), u).doit())

			generator_list = [
				phi, 1/phi,
				diff(phi, x), diff(phi, t),
				diff(phi, x, x), diff(phi, x, t), diff(phi, t, t),
				diff(phi, x, x, x), diff(phi, x, x, t), diff(phi, x, t, t), diff(phi, t, t, t),
				diff(phi, x, x, x, x), diff(phi, x, x, x, t), diff(phi, x, x, t, t), diff(phi, x, t, t, t), diff(phi, t, t, t, t),
				diff(phi, x, x, x, x, x), diff(phi, x, x, x, x, t), diff(phi, x, x, x, t, t), diff(phi, x, x, t, t, t), diff(phi, x, t, t, t, t), diff(phi, t, t, t, t, t)
				]

			KdVBurgersPainleveExpansion_polyInTermsOfPhi = Poly(KdVBurgersPainleveExpansion, gens = generator_list)
			KdVBurgersPainleveExpansion_polyInTermsOfPhi = KdVBurgersPainleveExpansion_polyInTermsOfPhi.as_expr()

			Equation0 = simplify(KdVBurgersPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 5))
			factors_list_of_equation0 = factor_list(Equation0)
			Equation0 = factors_list_of_equation0[-1][-1][0]
			u0_solution = solve(Equation0, u0)[0]
			Equation0_residual = Equation0.subs(u0, u0_solution)
			self.string_storer = str(u0_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			Equation1 = expand(KdVBurgersPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 4))
			Equation1 = expand(Equation1.subs(u0, u0_solution).doit())
			factors_list_of_equation1 = factor_list(Equation1)
			Equation1 = factors_list_of_equation1[-1][-1][0]
			Equation1Storage = Equation1
			u1_solution = solve(Equation1, u1)[0]
			u1_solution = together(u1_solution).doit()
			Equation1_residual = Equation1.subs(u1, u1_solution)
			Equation1_residual, _ = fraction(together(Equation1_residual).doit())
			Equation1_residual = expand(Equation1_residual)
			self.string_storer = str(u1_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			Equation2 = expand(KdVBurgersPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 3))
			Equation2 = expand(Equation2.subs(u0, u0_solution).doit())
			Equation2Storage = Equation2
			Equation2 = expand(Equation2.subs(u1, u1_solution).doit())
			Equation2, _ = fraction(together(Equation2).doit())
			Equation2 = expand(Equation2)
			factors_list_of_equation2 = factor_list(Equation2)
			Equation2 = factors_list_of_equation2[-1][-1][0]
			u2_solution = solve(Equation2, u2)[0]
			Equation2_residual = Equation2.subs(u2, u2_solution)
			self.string_storer = str(u2_solution)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			Equation3 = expand(KdVBurgersPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 2))
			Equation3 = expand(Equation3.subs(u0, u0_solution).doit())
			Equation3 = expand(Equation3.subs(u1, u1_solution).doit())
			Equation3_compatibility_test = Equation3
			Equation3 = expand(Equation3.subs(u2, u2_solution).doit())
			Equation3, _ = fraction(together(Equation3).doit())
			Equation3 = expand(Equation3)
			factors_list_of_equation3 = factor_list(Equation3)
			Equation3 = factors_list_of_equation3[-1][-1][0]
			self.string_storer = str(Equation3_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			Equation4 = expand(KdVBurgersPainleveExpansion_polyInTermsOfPhi.coeff(1/phi, 1))
			Equation4 = expand(Equation4.subs(u0, u0_solution).doit())
			Equation4 = expand(Equation4.subs(u1, u1_solution).doit())
			Equation4_compatibility_test = Equation4
			Equation4 = expand(Equation4.subs(u2, u2_solution).doit())
			Equation4, _ = fraction(together(Equation4).doit())
			Equation4 = expand(Equation4)
			factors_list_of_equation4 = factor_list(Equation4)
			Equation4 = factors_list_of_equation4[-1][-1][0]
			self.string_storer = str(Equation4_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			Equation5 = expand(KdVBurgersPainleveExpansion_polyInTermsOfPhi.coeff(phi, 0))
			Equation5 = expand(Equation5.subs(u0, u0_solution).doit())
			Equation5 = expand(Equation5.subs(u1, u1_solution).doit())
			Equation5_compatibility_test = Equation5
			factors_list_of_equation5 = factor_list(Equation5)
			Equation5 = factors_list_of_equation5[-1][-1][0]
			self.string_storer = str(Equation5_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			context_free_compatible_store = self.ContextFreeDifferentialCommutationCompatibilitySystemPDEs(
				[
					Equation2,
					Equation3_compatibility_test,
					Equation4_compatibility_test,
					Equation5_compatibility_test
					],
				phi)
			self.string_storer = str(context_free_compatible_store)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			self.context_free_compatibility_test = context_free_compatible_store[0]

			self.string_storer = str(Equation2)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			self.string_storer = str(Equation3_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			self.string_storer = str(Equation5_compatibility_test)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			phit_solution = solve(Equation2, diff(phi, t))[0]

			dictionary_of_solutions = {}
			dictionary_of_solutions[diff(phi, t)] = phit_solution
			dictionary_of_solutions[diff(phi, x, t)] = diff(phit_solution,x)
			dictionary_of_solutions[diff(phi, x, t)] = expand(dictionary_of_solutions[diff(phi, x, t)].doit())

			solution = solve(
				[
					diff(phi, t)-dictionary_of_solutions[diff(phi, t)],
					diff(phi, x, t)-dictionary_of_solutions[diff(phi, x, t)],
					Equation3_compatibility_test
				], [diff(phi, t), diff(phi, x, t), diff(phi, (x, 4))])

			dictionary_of_solutions[diff(phi, t)] = solution[diff(phi, t)]
			dictionary_of_solutions[diff(phi, x, t)] = solution[diff(phi, x, t)]
			dictionary_of_solutions[diff(phi, (x, 4))] = solution[diff(phi, (x, 4))]
			u2t_solution = solve(Equation5_compatibility_test, diff(u2, t))[0]
			dictionary_of_solutions[diff(u2, t)] = u2t_solution

			self.string_storer = str(dictionary_of_solutions[diff(phi, t)])
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			self.string_storer = str(dictionary_of_solutions[diff(phi, x, t)])
			Clock.schedule_once(self.scheduleOnceTester, 0)

			self.string_storer = str(dictionary_of_solutions[diff(phi, (x, 4))])
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			self.string_storer = str(dictionary_of_solutions[diff(u2, t)])
			Clock.schedule_once(self.scheduleOnceTester, 0)
			
			max_number_of_derivatives = 10

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
				self.string_storer = str(equation_new)
				Clock.schedule_once(self.scheduleOnceTester, 0)
				dictionary_of_solutions[diff(phi,(x,k))] = equation_new

			for k in range(1, max_number_of_derivatives):
				if k == 1:
					equation_new = diff(dictionary_of_solutions[diff(u2,t)],x)
				else:
					equation_new = diff(dictionary_of_solutions[diff(u2,(x,k-1),t)],x)
				equation_new = expand(equation_new)
				print('----------------------------------------------------------------------')
				if k == 1:
					print('\\u_{2, {x, t } } =')
				else:
					print('\\u_{2, {(x, '+str(k)+'), t} } =')
				print('----------------------------------------------------------------------')
				print_latex(equation_new)
				print('----------------------------------------------------------------------')
				self.string_storer = str(equation_new)
				Clock.schedule_once(self.scheduleOnceTester, 0)
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

			for term in dictionary_of_solutions:
				self.string_storer = str(dictionary_of_solutions[term])
				Clock.schedule_once(self.scheduleOnceTester, 0)

			equation1 = diff(dictionary_of_solutions[diff(phi,x,x,x,t)],x)
			equation1 = expand(equation1)
			equation1 = equation1.subs(diff(phi,(x,4)),dictionary_of_solutions[diff(phi,(x,4))])
			equation1 = expand(equation1)
			self.string_storer = str(equation1)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			equation2 = diff(dictionary_of_solutions[diff(phi,x,x,x,x)],t)
			recursion_track = 0
			while True:
				recursion_track = recursion_track + 1
				for term in dictionary_of_solutions:
					equation2 = equation2.subs(term,dictionary_of_solutions[term])
					equation2 = together(equation2).doit()
					equation2 = expand(equation2)
				numerator, denominator = fraction(together(equation2).doit())
				numerator = expand(numerator)
				order_param = ode_order(numerator, phi)
				self.string_storer = str(order_param)
				Clock.schedule_once(self.scheduleOnceTester, 0)
				if order_param < 4: break

			self.string_storer = str(equation2)
			Clock.schedule_once(self.scheduleOnceTester, 0)

			compatibility_check = equation2 - equation1
			compatibility_check = expand(compatibility_check)
			self.string_storer = str(compatibility_check)
			Clock.schedule_once(self.scheduleOnceTester, 0)
			self.context_sensitive_compatibility_test = int(compatibility_check == 0)
			self.compatibility_test = int(compatibility_check == 0)

	def show_dialog(self):
		if not self.dialog:
			self.dialog = MDDialog(
				title = 'Some calculation steps',
				type = 'custom',
				content_cls = Content()
				)
		# Introduce string storers and clock schedulers in appropriate places ....
		thread = threading.Thread(target=self.calculate_compatibility)
		thread.start()
		self.dialog.open()

	def close_handler(self):
		if self.dialog:
			self.dialog.dismiss()
			if self.compatibility_test == 0:
				self.root.ids.integrability_handler.text = 'The equation '+self.root.ids.input_pde.text+' is non-integrable'
			else:
				self.root.ids.integrability_handler.text = 'The equation '+self.root.ids.input_pde.text+' is integrable'
				self.root.ids.integrability_handler.text += '\nwith the auto-Backlund transform '+str(self.backlund_transform)
				self.root.ids.integrability_handler.text += '\n with the following compatible set of conditions\n'
				self.root.ids.integrability_handler.text += str(self.associated_conditions)
				self.string_storer = self.root.ids.integrability_handler.text
		if hasattr(self, 'schedule_interval'):
			Clock.unschedule(self.scheduleOnceTester)
		self.dialog = None

PainleveBacklundCheckApp().run()
