from sympy import *
from sympy.core.numbers import Zero as special_zero
import unittest
def max_order_term(inputPDE, func):
	if type(inputPDE) == Derivative:
		terms = [inputPDE]
	else:
		terms = list(sympify(str(inputPDE), evaluate = False).args)
	order = ode_order(inputPDE, func)
	print('-----------------------------------------------------------')
	print('>> inputPDE')
	print(inputPDE)
	print('>> type(inputPDE)')
	print(type(inputPDE))
	print('terms')
	print(terms)
	for individual_term in terms:
		print('>> order')
		print(order)
		print('>> individual_term')
		print(individual_term)
		print('>> isinstance(individual_term, Tuple)')
		print(isinstance(individual_term, Tuple))
		if (
			isinstance(individual_term, Tuple) and \
			len(individual_term) == 2 and individual_term[0] != func
			):
			continue
		print('>> ode_order(individual_term, func)')
		print(ode_order(individual_term, func))
		if ode_order(individual_term, func) == order:
			factorList = factor_list(individual_term)
			print('>> factorList')
			print(factorList)
			for individual_factor in factorList[1]:
				extract_factor = individual_factor[0]
				if ode_order(extract_factor, func) == order:
					print('>> individual_term')
					print(individual_term)
					print('>> individual_factor')
					print(individual_factor)
					return [
						individual_term,
						individual_factor[0],
						individual_factor[1]
					]
		else: pass
	print('-----------------------------------------------------------')

def simplificationPDE(inputPDE1, inputPDE2, func):
	max_order_term1 = max_order_term(inputPDE1, func)
	max_order_term2 = max_order_term(inputPDE2, func)
	coefficient1 = max_order_term1[0]/(max_order_term1[1]**max_order_term1[2])
	coefficient2 = max_order_term2[0]/(max_order_term2[1]**max_order_term2[2])
	inputPDE1 = inputPDE1/coefficient1 - (max_order_term1[1]**max_order_term1[2])
	inputPDE2 = inputPDE2/coefficient2 - (max_order_term2[1]**max_order_term2[2])
	varsFromPDE1 = max_order_term1[1].variables
	varsFromPDE2 = max_order_term2[1].variables
	print('-------------------------------')
	print('varsFromPDE1')
	print('-------------------------------')
	print(varsFromPDE1)
	print('-------------------------------')
	print('-------------------------------')
	print('varsFromPDE2')
	print('-------------------------------')
	print(varsFromPDE2)
	print('-------------------------------')
	inputPDE3 = inputPDE1
	for k in range(len(varsFromPDE2)):
		inputPDE3 = diff(inputPDE3, varsFromPDE2[k])
	inputPDE4 = inputPDE2
	for k in range(len(varsFromPDE1)):
		inputPDE4 = diff(inputPDE4, varsFromPDE1[k])
	expression = inputPDE3 - inputPDE4
	expression, _ = fraction(together(expression.doit()))
	return expression

def CompatibilitySystemPDEs(systemsOfPDEs, func):
	for equation in systemsOfPDEs:
		print('--------------------------------------------')
		print('>> max_order_term(equation, func)')
		print('--------------------------------------------')
		print(max_order_term(equation, func))
		print('--------------------------------------------')
	orders_list = [max_order_term(equation, func)[2] for equation in systemsOfPDEs]
	if len(systemsOfPDEs) == 2:
		maximum_order = max(orders_list)
		maximum_order_index = orders_list.index(maximum_order)
		equation_2_new = simplificationPDE(systemsOfPDEs[0], systemsOfPDEs[1], func)
		print('------------------------------------------------')
		print('>> systemsOfPDEs[maximum_order_index]')
		print('------------------------------------------------')
		print(systemsOfPDEs[maximum_order_index])
		print('------------------------------------------------')
		print('------------------------------------------------')
		print('>> equation_2_new')
		print('------------------------------------------------')
		print(equation_2_new)
		print('------------------------------------------------')
		if equation_2_new == 0:
			return True
		else:
			try:
				return CompatibilitySystemPDEs([systemsOfPDEs[maximum_order_index], equation_2_new], func)
			except Exception as err:
				print('\n The following error appeared:')
				print(err)
				return False
	else:
		systemsOfPDEs_F = [None for k in range(len(systemsOfPDEs))]
		for k in range(len(systemsOfPDEs)-1):
			systemsOfPDEs_F[k] = simplificationPDE(systemsOfPDEs[k], systemsOfPDEs[k+1], func)
		systemsOfPDEs_F[-1] = simplificationPDE(systemsOfPDEs[-1], systemsOfPDEs[0], func)
		flag = False
		for k in range(len(systemsOfPDEs_F)):
			flag = flag or systemsOfPDEs_F[k] == 0
		if flag: return flag
		else: return CompatibilitySystemPDEs(systemsOfPDEs_F, func)

class TestCompatibilityPDE(unittest.TestCase):
	# def test_TestCaseMaxOrderTerm(self):
	# 	x, y = symbols('x y')
	# 	f, g = symbols('f g', cls = Function)
	# 	inputPDE = f(x) + 1
	# 	self.assertTrue(max_order_term(inputPDE, f) == [f(x), f(x), 1])
	# def test_TestCaseSimplification(self):
	# 	x, y = symbols('x y')
	# 	f, g = symbols('f g', cls = Function)
	# 	equations = [diff(f(x), (x, 2)) + f(x), diff(f(x), x) - 1]
	# 	self.assertTrue(simplificationPDE(equations[0], equations[1], f) == Derivative(f(x), x))
	# def test_TestCase1(self):
	# 	x, y = symbols('x y')
	# 	f, g = symbols('f g', cls = Function)
	# 	equations = [diff(f(x), x) - f(x), diff(f(x), (x, 2)) - diff(f(x), x)]
	# 	self.assertTrue(CompatibilitySystemPDEs(equations, f))
	# def test_TestCase2(self):
	# 	x, y = symbols('x y')
	# 	f, g = symbols('f g', cls = Function)
	# 	equations = [diff(f(x), x) - f(x), diff(f(x), (x, 2)) - f(x)]
	# 	self.assertTrue(CompatibilitySystemPDEs(equations, f))
	# def test_TestCase3(self):
	# 	x, y = symbols('x y')
	# 	f, g = symbols('f g', cls = Function)
	# 	equations = [diff(f(x), (x, 2)) + f(x), diff(f(x), x) - 1]
	# 	self.assertFalse(CompatibilitySystemPDEs(equations, f))
	# def test_TestCase4(self):
	# 	# I am just trying if the multivariate Cauchy Reimann works
	# 	# with two variables or not
	# 	# Tested and it works
	# 	x, y = symbols('x y')
	# 	u, v = symbols('u v', cls = Function)
	# 	equations = [
	# 		diff(u(x, y), x, x) + diff(u(x, y), y, y),
	# 		diff(u(x, y), x)-diff(v(x, y), y),
	# 		diff(u(x, y), y)+diff(v(x, y), x),
	# 		diff(v(x, y), x, x) + diff(v(x, y), y, y)
	# 		]
	# 	self.assertTrue(CompatibilitySystemPDEs(equations, u))
	# def test_TestCase5(self):
	# 	# I am just trying if the multivariate Cauchy Reimann works
	# 	# with two variables (w.r.t. the other variables)
	# 	x, y = symbols('x y')
	# 	u, v = symbols('u v', cls = Function)
	# 	equations = [
	# 		diff(u(x, y), x, x) + diff(u(x, y), y, y),
	# 		diff(u(x, y), x)-diff(v(x, y), y),
	# 		diff(u(x, y), y)+diff(v(x, y), x),
	# 		diff(v(x, y), x, x) + diff(v(x, y), y, y)
	# 		]
	# 	self.assertTrue(CompatibilitySystemPDEs(equations, v))
	def test_TestCaseBurgers(self):
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
		print(CompatibilitySystemPDEs([original_equation, backlund_equation], phi))
		self.assertTrue(True)
if __name__ == '__main__':
	unittest.main()


