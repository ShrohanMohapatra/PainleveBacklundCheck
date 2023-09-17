from sympy import *
from collections import Counter
import unittest
x, t, sigma = symbols('x t sigma')
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

def simplificationPDE(inputPDE1, inputPDE2, func):
	max_order_term1 = max_order_term(inputPDE1, func)
	max_order_term2 = max_order_term(inputPDE2, func)
	coefficient1 = max_order_term1[0] # max_order_term1[0]/(max_order_term1[1]**max_order_term1[2])
	coefficient2 = max_order_term2[0] # max_order_term2[0]/(max_order_term2[1]**max_order_term2[2])
	inputPDE1 = inputPDE1/coefficient1 - max_order_term1[1] # (max_order_term1[1]**max_order_term1[2])
	inputPDE2 = inputPDE2/coefficient2 - max_order_term2[1] # (max_order_term2[1]**max_order_term2[2])
	allAvailableVariables = func.free_symbols
	varsFromPDE1 = list(max_order_term1[2].variables) if max_order_term1[3] != 0 else []
	varsFromPDE2 = list(max_order_term2[2].variables) if max_order_term2[3] != 0 else []

	# print('-------------------------------')
	# print('varsFromPDE1')
	# print('-------------------------------')
	# print(varsFromPDE1)
	# print('-------------------------------')
	# print('-------------------------------')
	# print('varsFromPDE2')
	# print('-------------------------------')
	# print(varsFromPDE2)
	# print('-------------------------------')

	dictFromPDE1 = Counter(varsFromPDE1)
	dictFromPDE2 = Counter(varsFromPDE2)

	varsFromPDE3 = list((dictFromPDE2-dictFromPDE1).keys())
	varsFromPDE4 = list((dictFromPDE1-dictFromPDE2).keys())

	# print('-------------------------------')
	# print('varsFromPDE3')
	# print('-------------------------------')
	# print(varsFromPDE3)
	# print('-------------------------------')
	# print('-------------------------------')
	# print('varsFromPDE4')
	# print('-------------------------------')
	# print(varsFromPDE4)
	# print('-------------------------------')

	inputPDE3 = inputPDE1
	for variable in varsFromPDE3:
		inputPDE3 = diff(inputPDE3, variable)

	inputPDE4 = inputPDE2
	for variable in varsFromPDE4:
		inputPDE4 = diff(inputPDE4, variable)

	expression = inputPDE3 - inputPDE4
	expression, _ = fraction(together(expression.doit()))
	expression = expand(expression)
	return expression

def CompatibilitySystemPDEs(systemsOfPDEs, func, depth_max, depth = 0):

	print('------------------------------------------------------------------------------------------------')
	print('Welcome to recursion number = ',depth+1)
	print('------------------------------------------------------------------------------------------------')

	for k in range(len(systemsOfPDEs)):
		print('--------------------------------------------------------------------------------------------')
		print('>> PDE Number '+str(k+1)+' -->')
		print('--------------------------------------------------------------------------------------------')
		print(systemsOfPDEs[k])
		print('--------------------------------------------------------------------------------------------')

	# for equation in systemsOfPDEs:
	# 	print('--------------------------------------------')
	# 	print('>> max_order_term(equation, func)')
	# 	print('--------------------------------------------')
	# 	print([equation], max_order_term(equation, func))
	# 	print('--------------------------------------------')

	# print('Depth = ', depth)

	orders_list = [max_order_term(equation, func)[3] for equation in systemsOfPDEs]
	if len(systemsOfPDEs) == 2:
		maximum_order = max(orders_list)
		maximum_order_index = orders_list.index(maximum_order)
		equation_2_new = simplificationPDE(systemsOfPDEs[0], systemsOfPDEs[1], func)

		# print('------------------------------------------------')
		# print('>> systemsOfPDEs[maximum_order_index]')
		# print('------------------------------------------------')
		# print(systemsOfPDEs[maximum_order_index])
		# print('------------------------------------------------')
		# print('------------------------------------------------')
		# print('>> equation_2_new')
		# print('------------------------------------------------')
		# print(equation_2_new)
		# print('------------------------------------------------')

		if equation_2_new == 0:
			return [depth, True]
		else:
			if depth >= depth_max: return [depth, False]
			try:
				return CompatibilitySystemPDEs([systemsOfPDEs[maximum_order_index], equation_2_new], func, depth_max, depth + 1)
			except Exception:
				return [depth, False]
	else:
		systemsOfPDEs_F = [None for k in range(len(systemsOfPDEs))]
		for k in range(len(systemsOfPDEs)-1):
			systemsOfPDEs_F[k] = simplificationPDE(systemsOfPDEs[k], systemsOfPDEs[k+1], func)
		systemsOfPDEs_F[-1] = simplificationPDE(systemsOfPDEs[-1], systemsOfPDEs[0], func)

		# print('-------------------')
		# print(systemsOfPDEs_F)
		# print('-------------------')

		flag = False
		for k in range(len(systemsOfPDEs_F)):
			flag = flag or systemsOfPDEs_F[k] == 0
		if flag: return [depth, flag]
		else:
			if depth >= depth_max: return [depth, False]
			try:
				return CompatibilitySystemPDEs(systemsOfPDEs_F, func, depth_max, depth + 1)
			except Exception:
				return [depth, False]

class TestCompatibilityPDE(unittest.TestCase):

	def test_TestCaseMultivariateCauchyRiemannVersion1(self):
		# I am just trying if the multivariate Cauchy Reimann works
		# with two variables or not
		# Tested and it works
		print('------------------------------------------------------------------------------------------------------------------------------------')
		print('Here I am testing for the compatibility of the Cauchy-Reimann system of PDEs')
		print('in functions <u(x, y), v(x, y)> by recursing through the function u ....')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		x, y = symbols('x y')
		u, v = symbols('u v', cls = Function)
		equations = [
			diff(u(x, y), x, x) + diff(u(x, y), y, y),
			diff(u(x, y), x)-diff(v(x, y), y),
			diff(u(x, y), y)+diff(v(x, y), x),
			diff(v(x, y), x, x) + diff(v(x, y), y, y)
			]
		compatibility_checker = CompatibilitySystemPDEs(equations, u(x, y), 10)
		depth = compatibility_checker[0]
		status = compatibility_checker[1]
		if status:
			print('The Cauchy-Riemann system 1 renders compatible')
		else:
			print('The Cauchy-Riemann system 1 renders incompatible')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		# print('Number of cycles of reduction in Cauchy Riemann system 1', depth)
		self.assertTrue(status)

	def test_TestCaseMultivariateCauchyRiemannVersion2(self):
		# I am just trying if the multivariate Cauchy Reimann works
		# with two variables (w.r.t. the other variables)
		print('------------------------------------------------------------------------------------------------------------------------------------')
		print('Here I am testing for the compatibility of the Cauchy-Reimann system of PDEs')
		print('in functions <u(x, y), v(x, y)> by recursing through the function v ....')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		x, y = symbols('x y')
		u, v = symbols('u v', cls = Function)
		equations = [
			diff(u(x, y), x, x) + diff(u(x, y), y, y),
			diff(u(x, y), x)-diff(v(x, y), y),
			diff(u(x, y), y)+diff(v(x, y), x),
			diff(v(x, y), x, x) + diff(v(x, y), y, y)
			]
		compatibility_checker = CompatibilitySystemPDEs(equations, u(x, y), 10)
		depth = compatibility_checker[0]
		status = compatibility_checker[1]
		if status:
			print('The Cauchy-Riemann system 2 renders compatible')
		else:
			print('The Cauchy-Riemann system 2 renders incompatible')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		# print('Number of cycles of reduction in Cauchy Riemann system 2', depth)
		self.assertTrue(status)

	def test_TestCaseBurgers(self):
		print('------------------------------------------------------------------------------------------------------------------------------------')
		print('Here I am testing for the compatibility of the system of PDEs')
		print('that come from the WTC Painleve analysis of the Burgers equation ....')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		x, t, sigma = symbols('x t sigma')
		phi, u0 = symbols('phi u0', cls = Function)
		u = (phi(x,t)*(sigma*diff(phi(x, t), (x, 2)) - diff(phi(x, t), t)) - 2*sigma*diff(phi(x, t), x)**2)/phi(x, t)/diff(phi(x, t), x)
		u2 = (sigma*diff(phi(x, t), (x, 2)) - diff(phi(x, t), t))/diff(phi(x, t), x)
		original_equation = diff(u, t) + u*diff(u, x) - sigma*diff(u, (x, 2))
		original_equation, _ = fraction(together(original_equation.doit()))
		original_equation = expand(original_equation)
		backlund_equation = diff(u2, t) + u2*diff(u2, x) - sigma*diff(u2, (x, 2))
		backlund_equation, _ = fraction(together(backlund_equation.doit()))
		backlund_equation = expand(backlund_equation)
		compatibility_checker = CompatibilitySystemPDEs([original_equation, backlund_equation], phi(x, t), 15)
		depth = compatibility_checker[0]
		status = compatibility_checker[1]
		if status:
			print('The PDE system from WTC Painleve analysis of the Burgers equation renders compatible')
		else:
			print('The PDE system from WTC Painleve analysis of the Burgers equation renders incompatible')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		# print('Number of cycles of reduction in Burgers equation', depth)
		self.assertTrue(status)

	def test_TestCaseKdV(self):
		print('------------------------------------------------------------------------------------------------------------------------------------')
		print('Here I am testing for the compatibility of the system of PDEs')
		print('that come from the WTC Painleve analysis of the KdV equation ....')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		x, t, sigma = symbols('x t sigma')
		phi, u2 = symbols('phi u2', cls = Function)
		u = 12*sigma*(phi(x, t)*diff(phi(x, t), (x, 2)) - diff(phi(x, t), x)**2)/phi(x, t)/diff(phi(x, t), x) + u2(x, t)
		equations = [None]*4
		equations[0] = diff(phi(x, t), x)*diff(phi(x, t), t) + diff(phi(x, t), x)**2*u2(x, t) + 4*sigma*diff(phi(x, t), x)*diff(phi(x, t), (x, 3)) - 3*sigma*diff(phi(x, t), (x, 2))**2
		equations[1] = diff(phi(x, t), x, t) + diff(phi(x, t), x, x)*u2(x, t) + sigma*diff(phi(x, t), (x, 4))
		equations[2] = diff(u2(x, t), t) + u2(x, t)*diff(u2(x, t), x) + sigma*diff(u2(x, t), (x, 3))
		equations[3] = diff(u, t) + u*diff(u, x) + sigma*diff(u, (x, 3))
		equations[3], _ = fraction(together(equations[3].doit()))
		equations[3] = expand(equations[3])
		compatibility_checker = CompatibilitySystemPDEs(equations, phi(x, t), 10)
		depth = compatibility_checker[0]
		status = compatibility_checker[1]
		if status:
			print('The PDE system from WTC Painleve analysis of the KdV equation renders compatible')
		else:
			print('The PDE system from WTC Painleve analysis of the KdV equation renders incompatible')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		# print('Number of cycles of reduction in KdV equation', depth)
		self.assertTrue(status)

	def test_TestCaseBBM(self):
		print('------------------------------------------------------------------------------------------------------------------------------------')
		print('Here I am testing for the compatibility of the system of PDEs')
		print('that come from the WTC Painleve analysis of the BBM (Benjamin-Bona-Mahony) equation ....')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		u0, u1, u2, phi = symbols('u0 u1 u2 phi', cls = Function)
		x, t, sigma = symbols('x t sigma')
		u = u0(x, t)/phi(x, t)**2 + u1(x, t)/phi(x, t) + u2(x, t)
		bbm_equation = phi(x, t)**5*(diff(u, t) + diff(u, x) + u*diff(u, x) - diff(u, x, x, t))
		bbm_equation, _ = fraction(together(bbm_equation.doit()))
		bbm_equation = expand(bbm_equation)
		bbm_equation_coeffs = [bbm_equation.coeff(phi(x, t), k+1) for k in range(5)]
		bbm_equation_with_all_phi_powers = 0
		for k in range(5):
			bbm_equation_with_all_phi_powers = bbm_equation_with_all_phi_powers + bbm_equation_coeffs[k] * phi(x, t)**(k+1)
		bbm_equation_coeffs.insert(0, simplify(bbm_equation - bbm_equation_with_all_phi_powers))
		compatibility_checker = CompatibilitySystemPDEs(bbm_equation_coeffs, phi(x, t), 4)
		depth = compatibility_checker[0]
		status = compatibility_checker[1]
		if status:
			print('The PDE system from WTC Painleve analysis of the BBM equation renders compatible')
		else:
			print('The PDE system from WTC Painleve analysis of the BBM equation renders incompatible')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		# print('Number of cycles of reduction in BBM equation', depth)
		self.assertFalse(status)

	def test_TestCaseKuramotoSivashinsky(self):
		print('------------------------------------------------------------------------------------------------------------------------------------')
		print('Here I am testing for the compatibility of the system of PDEs')
		print('that come from the WTC Painleve analysis of the KS (Kuramoto-Sivashinsky) equation ....')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		u0, u1, u2, u3, phi = symbols('u0 u1 u2 u3 phi', cls = Function)
		x, t, sigma = symbols('x t sigma')
		u = u0(x, t)/phi(x, t)**3 + u1(x, t)/phi(x, t)**2 + u2(x, t)/phi(x, t) + u3(x, t)
		ks_equation = phi(x, t)**7*(diff(u, t) + u*diff(u, x) + diff(u, x, x) + diff(u, x, x, x, x))
		ks_equation, _ = fraction(together(ks_equation.doit()))
		ks_equation = expand(ks_equation)
		ks_equation_coeffs = [ks_equation.coeff(phi(x, t), k+1) for k in range(7)]
		ks_equation_with_all_phi_powers = 0
		for k in range(7):
			ks_equation_with_all_phi_powers = ks_equation_with_all_phi_powers + ks_equation_coeffs[k] * phi(x, t)**(k+1)
		ks_equation_coeffs.insert(0, simplify(ks_equation - ks_equation_with_all_phi_powers))
		compatibility_checker = CompatibilitySystemPDEs(ks_equation_coeffs, phi(x, t), 4)
		depth = compatibility_checker[0]
		status = compatibility_checker[1]
		if status:
			print('The PDE system from WTC Painleve analysis of the KS equation renders compatible')
		else:
			print('The PDE system from WTC Painleve analysis of the KS equation renders incompatible')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		# print('Number of cycles of reduction in BBM equation', depth)
		self.assertFalse(status)

	def test_TestCaseKdVBurgers(self):
		print('------------------------------------------------------------------------------------------------------------------------------------')
		print('Here I am testing for the compatibility of the system of PDEs')
		print('that come from the WTC Painleve analysis of the KdV-Burgers equation ....')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		u0, u1, u2, u3, phi = symbols('u0 u1 u2 u3 phi', cls = Function)
		x, t, sigma = symbols('x t sigma')
		u = u0(x, t)/phi(x, t)**2 + u1(x, t)/phi(x, t) + u2(x, t)
		kdv_burgers_equation = phi(x, t)**5*(diff(u, t) + u*diff(u, x) + diff(u, x, x, x) - sigma**2 * diff(u, x, x))
		kdv_burgers_equation, _ = fraction(together(kdv_burgers_equation.doit()))
		kdv_burgers_equation = expand(kdv_burgers_equation)
		kdv_burgers_equation_coeffs = [kdv_burgers_equation.coeff(phi(x, t), k+1) for k in range(5)]
		kdv_burgers_equation_with_all_phi_powers = 0
		for k in range(5):
			kdv_burgers_equation_with_all_phi_powers = kdv_burgers_equation_with_all_phi_powers + kdv_burgers_equation_coeffs[k] * phi(x, t)**(k+1)
		kdv_burgers_equation_coeffs.insert(0, simplify(kdv_burgers_equation - kdv_burgers_equation_with_all_phi_powers))
		compatibility_checker = CompatibilitySystemPDEs(kdv_burgers_equation_coeffs, phi(x, t), 7)
		depth = compatibility_checker[0]
		status = compatibility_checker[1]
		if status:
			print('The PDE system from WTC Painleve analysis of the KdV-Burgers equation renders compatible')
		else:
			print('The PDE system from WTC Painleve analysis of the KdV-Burgers equation renders incompatible')
		print('------------------------------------------------------------------------------------------------------------------------------------')
		# print('Number of cycles of reduction in BBM equation', depth)
		self.assertFalse(status)

if __name__ == '__main__':
	unittest.main()
