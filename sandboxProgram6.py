from sympy import *
from TransformExistence import *
from ContextFreeDifferentialCommutationCompatibility import *
from time import time

start = time()

x, t, sigma = symbols('x t sigma')
f = Function('f')(x, t)
try:
	phi = Function('phi')(x, t)
except: pass

Burgers_Equation = parse_expr('diff(f(x,t),t)+f(x,t)*diff(f(x,t),x)-sigma*diff(f(x,t),(x,2))')
result_package = ultimatePainleve(Burgers_Equation, x, t)
print('------------------------------------------------------------------')
print('Final answer of the sandbox program with Burgers equation')
print('------------------------------------------------------------------')

print(result_package[0])
print(result_package[1])
Number_Of_Sets = len(result_package[2])
alpha = result_package[0]
u_max = Function('U'+str(-alpha))(x, t)

for k in range(Number_Of_Sets):
	print('------------------------------------------------------------------')
	print(result_package[2][k])
	print('------------------------------------------------------------------')
	for m in range(len(result_package[2][k])):
		factors_list_of_equation = factor_list(result_package[2][k][m])
		result_package[2][k][m] = factors_list_of_equation[-1][-1][0]
		result_package[2][k][m] = expand(result_package[2][k][m])

flag_switch = False
flag_on_phi = False
flag_on_umax = False

# for k in range(len(result_package[2])-1,-1,-1):
# 	context_free_compatible_store = ContextFreeDifferentialCommutationCompatibilitySystemPDEs(result_package[2][k], phi)
# 	if context_free_compatible_store[1]:
# 		flag_switch = True
# 		flag_on_phi = True
# 		break
# 	context_free_compatible_store = ContextFreeDifferentialCommutationCompatibilitySystemPDEs(result_package[2][k], u_max)
# 	if context_free_compatible_store[1]:
# 		flag_switch = True
# 		flag_on_umax = True
# 		break

k = -1
context_free_compatible_store = ContextFreeDifferentialCommutationCompatibilitySystemPDEs(result_package[2][k], phi)
if context_free_compatible_store[1]:
	flag_switch = True
	flag_on_phi = True
else:
	if ode_order(Burgers_Equation, f) > 2:
		context_free_compatible_store = ContextFreeDifferentialCommutationCompatibilitySystemPDEs(result_package[2][k], u_max)
		if context_free_compatible_store[1]:
			flag_switch = True
			flag_on_umax = True

if not(flag_switch):
	print('The system of PDEs is not context-free compatible')
else:
	print('The system of PDEs is context-free compatible as 0 belongs to g_{ '+str(context_free_compatible_store[0]+1)+'}')
	print('Equation -> '+str(result_package[2][k]))
	if flag_on_phi: print('It is based on:', str(phi))
	if flag_on_umax: print('It is based on:', str(u_max))

end = time()

print(' Time consumed => ', end-start)

