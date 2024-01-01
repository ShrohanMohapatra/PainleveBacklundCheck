from sympy import *
from TransformExistence import *
x, t, sigma = symbols('x t sigma')
u = Function('f')(x, t)
BurgersEquation = parse_expr('diff(f(x,t),t)+f(x,t)*diff(f(x,t),x)-sigma*diff(f(x,t),(x,2))')
result_package = ultimatePainleve(BurgersEquation, x, t)
print('------------------------------------------------------------------')
print('Final answer of the sandbox program with Burgers equation')
print('------------------------------------------------------------------')
print(result_package[0])
print(result_package[1])
for k in range(len(result_package[2])):
	print(result_package[2][k])
print('------------------------------------------------------------------')