from sympy import *
from TransformExistence import *
from ContextFreeDifferentialCommutationCompatibility import *
import threading

def ContextFreeCompatibilityThread(systemOfPDE, func, result_holder):
	result_holder.append(ContextFreeDifferentialCommutationCompatibilitySystemPDEs(systemOfPDE, func))

x, t, sigma = symbols('x t sigma')
f = Function('f')(x, t)
try:
	phi = Function('phi')(x, t)
except: pass
KdVEquation = parse_expr('diff(f(x,t),t)+f(x,t)*diff(f(x,t),x)+sigma*diff(f(x,t),(x,3))')
result_package = ultimatePainleve(KdVEquation, x, t)
print('------------------------------------------------------------------')
print('Final answer of the sandbox program with KdV equation')
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

result_holder_for_context_free_checks = {phi:[], u_max: []}
result_holder_from_the_threads_for_phi = [[] for k in range(len(result_package[2]))]
result_holder_from_the_threads_for_umax = [[] for k in range(len(result_package[2]))]

list_of_threads = [None for k in range(2*len(result_package[2]))]

for k in range(len(result_package[2])):
	list_of_threads[k] = threading.Thread(target = ContextFreeCompatibilityThread, args = (result_package[2][k], phi, result_holder_from_the_threads_for_phi[k]))
	list_of_threads[k+1] = threading.Thread(target = ContextFreeCompatibilityThread, args = (result_package[2][k], u_max, result_holder_from_the_threads_for_umax[k]))

for k in range(2*len(result_package[2])): list_of_threads[k].start()
for k in range(2*len(result_package[2])): list_of_threads[k].join()

for k in range(len(result_package[2])):
	try:
		if result_holder_from_the_threads_for_phi[k][0] == True:
			result_holder_for_context_free_checks[phi].append(result_holder_from_the_threads_for_phi[k][0])
			break
	except:
		try:
			if True in result_holder_from_the_threads_for_phi[k][0]:
				result_holder_for_context_free_checks[phi].append(result_holder_from_the_threads_for_phi[k][0])
				break
		except:
			result_holder_for_context_free_checks[phi].append(result_holder_from_the_threads_for_phi[k][0])
	try:
		if result_holder_from_the_threads_for_umax[k][0] == True:
			result_holder_for_context_free_checks[u_max].append(result_holder_from_the_threads_for_umax[k][0])
			break
	except:
		try:
			if True in result_holder_from_the_threads_for_umax[k][0]:
				result_holder_for_context_free_checks[u_max].append(result_holder_from_the_threads_for_umax[k][0])
				break
		except:
			result_holder_for_context_free_checks[u_max].append(result_holder_from_the_threads_for_umax[k][0])

for k in range(2*len(result_package[2])):
	if list_of_threads[k].is_alive():
		list_of_threads[k].join(timeout=0)
