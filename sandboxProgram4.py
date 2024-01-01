from sympy import *
from TransformExistence import *
from ContextFreeDifferentialCommutationCompatibility import *
import threading
import random
from time import time

def ContextFreeCompatibilityThread(systemOfPDE, func, semaphore, lock, event, result_holder):
	with semaphore:
		final_answer = ContextFreeDifferentialCommutationCompatibilitySystemPDEs(systemOfPDE, func)
		with lock:
			result_holder.append(final_answer)
		event.set()

start = time()

x, t, sigma = symbols('x t sigma')
f = Function('f')(x, t)
try:
	phi = Function('phi')(x, t)
except: pass
BurgersEquation = parse_expr('diff(f(x,t),t)+f(x,t)*diff(f(x,t),x)+sigma*diff(f(x,t),(x,3))')
result_package = ultimatePainleve(BurgersEquation, x, t)
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

semaphore = threading.Semaphore(value = 1)
event = threading.Event()
lock = threading.Lock()

result_holder_for_context_free_checks = {phi:[], u_max: []}
result_holder_from_the_threads_for_phi = [[] for k in range(len(result_package[2]))]
result_holder_from_the_threads_for_umax = [[] for k in range(len(result_package[2]))]

list_of_threads = [None for k in range(2*len(result_package[2]))]

for k in range(len(result_package[2])):
	list_of_threads[2*k] = threading.Thread(target = ContextFreeCompatibilityThread, args = (result_package[2][k], phi, semaphore, lock, event, result_holder_from_the_threads_for_phi[k]))
	list_of_threads[2*k+1] = threading.Thread(target = ContextFreeCompatibilityThread, args = (result_package[2][k], u_max, semaphore, lock, event, result_holder_from_the_threads_for_umax[k]))

random.shuffle(list_of_threads)

for k in range(2*len(result_package[2])): list_of_threads[k].start()
event.wait()
semaphore.release()
for k in range(2*len(result_package[2])): list_of_threads[k].join()

for k in range(len(result_package[2])):
	with lock:
		result_holder_for_context_free_checks[phi].append(result_holder_from_the_threads_for_phi[k][0])
		result_holder_for_context_free_checks[u_max].append(result_holder_from_the_threads_for_umax[k][0])

all_threads_finished = all(not thread.is_alive() for thread in list_of_threads)

# Now you can use 'all_threads_finished' to determine if all threads have finished
if all_threads_finished:
	print("All threads have finished.")
else:
	print("Some threads are still running.")

print(result_holder_for_context_free_checks)

end = time()

print(' Time consumed => ', end-start)

# fastest_thread = min(list_of_threads, key = lambda x: x._tstate_lock)

# list_of_threads.remove(fastest_thread)

# for k in range(len(list_of_threads)):
# 	if list_of_threads[k].is_alive():
# 		list_of_threads[k].join(timeout=0)
