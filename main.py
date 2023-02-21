from kivy.lang import Builder
from kivymd.app import MDApp
from kivymd.uix.floatlayout import MDFloatLayout
from kivymd.uix.button import MDRoundFlatButton
from kivy.properties import ObjectProperty
from kivy.uix.gridlayout import GridLayout
from kivy.uix.label import Label
from kivy.uix.popup import Popup
from sympy import *
from pprint import pprint as pretty_printer
init_printing()
x, t, sigma, b = symbols('x t sigma b')
f = Function('f')(x, t)
# Painleve property and Backlund transform
def totalTableOfCoeffsForPoly(expr, zeta):
	generator_list = [1/zeta, zeta]
	for m in range(10):
		for n in range(10):
			if m == 0 and n == 0: pass
			elif m == 0 and n == 1:
				generator_list.append(Derivative(zeta, t))
			elif m == 1 and n == 0:
				generator_list.append(Derivative(zeta, x))
			elif m == 1 and n == 1:
				generator_list.append(Derivative(zeta, x, t))
				generator_list.append(Derivative(zeta, t, x))
			elif m == 1 and n > 1:
				generator_list.append(Derivative(zeta, x, (t, n)))
			elif n == 1 and m > 1:
				generator_list.append(Derivative(zeta, t, (x, m)))
			else:
				generator_list.append(
					Derivative(zeta, (x, m), (t, n))
					)
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
	function_PDE, # Input PDE
	x, # The input variable used for the spatial variable
	t, # The input variable used for the temporal variable
	alpha # The power balancing exponent
	):
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
	totalTable = totalTableOfCoeffsForPoly(orig_eqn, phi(x, t))
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
		simplify(totalTable[k][1].subs([
		(U[m](x,t), U_final[m]) if m<(-alpha) else \
		((U[m](x,t), U[m](x,t)) if m==(-alpha) else (U[m](x,t), 0))\
		for m in range(M+1)
		]).doit()) \
		for k in range(len(totalTable))
	]
	refined_compatibility_Conditions = []
	for k in range(len(totalTable)):
		if compatibility_Conditions[k] != 0:
			refined_compatibility_Conditions.append(
				compatibility_Conditions[k]
				)
	# refined_compatibility_Conditions.append(
		# Eq(U[-alpha](x, t), phi(x, t)))
	return [backlund_transform, refined_compatibility_Conditions]
def ultimatePainleve(
	function_PDE, # Input PDE
	x, # The input variable used for the spatial variable
	t, # The input variable used for the temporal variable
	):
	range_of_exponents = [-1, -2, -3, -4]
	flag = False
	for alpha in range_of_exponents:
		M = int(-alpha + 3) # Maximum number of terms
		U = [Function('U'+str(k)) for k in range(M+1)]
		backlund_transform, compatibility_Conditions = PainlevePDE_withBacklundTransform(function_PDE, x, t, alpha)
		if len(compatibility_Conditions) < 2: continue
		else:
			flag = compatibility_Conditions[-1].subs(U[-alpha](x, t),f).doit() == expand(function_PDE)
			if not flag: continue
			else: break
	if alpha == -4 and len(compatibility_Conditions) == 1 and backlund_transform == U[4](x,t):
		return (None, backlund_transform, compatibility_Conditions)
	for k in range(len(compatibility_Conditions)):
		if str(compatibility_Conditions[k]).find(str(U[-alpha](x, t)))==-1:
			return (None, backlund_transform, [])				
	if not flag:
		alpha, compatibility_Conditions = None, []
	return (alpha, backlund_transform, compatibility_Conditions)
kivyBuilder = Builder.load_string(
"""
<PainleveBacklundCheck>
	MDLabel:
		id: painlevebacklund
		text: 'Painleve property and Backlund transforms'
		size_hint: (0.4, 0.1)
		pos_hint: {'x': 0.05, 'y': 0.85}
		font_size: '20sp'
		hint_font_size: '30sp'
	MDLabel:
		id: inputPDElabel
		text: 'Enter the PDE here =>'
		size_hint: (0.2, 0.1)
		pos_hint: {'x': 0.05, 'y': 0.7}
		font_size: '20sp'
		hint_font_size: '30sp'
	MDTextField:
		id: inputPDEText
		hint_text: 'Enter the PDE in the sympy-compatible format'
		size_hint: (0.6, 0.1)
		pos_hint: {'x': 0.35, 'y': 0.7}
		multiline: True
		font_size: '20sp'
		hint_font_size: 50
	MDRoundFlatButton:
		id: evaluateButton
		text: 'Evaluate'
		text_color: 'yellow'
		size_hint: (0.25, 0.05)
		pos_hint: {'x': 0.2, 'y': 0.575}
		on_release: app.evaluate_button()
		font_size: '20sp'
	MDRoundFlatButton:
		id: resetButton
		text: 'Reset'
		text_color: 'yellow'
		size_hint: (0.25, 0.05)
		pos_hint: {'x': 0.55, 'y': 0.575}
		on_release: app.reset_button()
		font_size: '20sp'
	MDTextField:
		id: outputText
		hint_text: 'You will see the output here'
		text: ''
		size_hint: (0.9, 0.5)
		pos_hint: {'x': 0.05, 'y': 0.05}
		multiline: True
		font_size: '20sp'
		hint_font_size: 50
""")
class PainleveBacklundCheck(MDFloatLayout):
	painlevebacklund = ObjectProperty(None)
class PainleveBacklundCheckApp(MDApp):
	def build(self):
		self.theme_cls.theme_style = "Dark"
		self.theme_cls.primary_palette = "Orange"
		return PainleveBacklundCheck()
	def evaluate_button(self):
		layout = GridLayout(cols = 1)
		popupLabel = Label(text = 'Calculations underway ....')
		layout.add_widget(popupLabel)
		popup = Popup(
			title='Let\'s wait for the result',
			content = layout, size_hint = (None, None),
			size = (400, 400)
			)
		popup.open()
		counter = 10**6
		while counter > 0: counter -= 1
		popup.dismiss()
		self.root.ids.outputText.text = self.root.ids.inputPDEText.text
		inputPDE = sympify(self.root.ids.inputPDEText.text)
		alpha, transform, compatibility = ultimatePainleve(
			inputPDE, # Input PDE
			x, # The input variable used for the spatial variable
			t # The input variable used for the temporal variable
		)
		if alpha == None:
			self.root.ids.outputText.text = 'The equation is non-integrable due to the failure of WTC test.'
		else:
			final_output = 'The power balancing exponent is alpha = '+\
							str(alpha)+'\n'
			final_output += 'The transform is given by U = '+\
							str(transform)+'\n'
			final_output += 'where the additional conditions are'+\
							' given by\n'
			for equation in compatibility:
				if equation != 0:
					final_output += '>'+str(Eq(equation, 0))+'\n'
			final_output += '\n Note that you might have to check the compatibility of the equations (except the last one)'
			self.root.ids.outputText.text = final_output
	def reset_button(self):
		self.root.ids.outputText.text = ''
		self.root.ids.inputPDEText.text = ''
if __name__ == "__main__":
    PainleveBacklundCheckApp().run()
