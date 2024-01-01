# I am just trying to see if the multithreading strategy
# works well with a dummy Collatz conjecture verification
import threading
def multithreadedRuns(array,x_variable):
	print(x_variable)
	if array[x_variable] is None:
		if x_variable == 0:
			array[x_variable] = 1
		if x_variable%2 == 0:
			array[x_variable] = x_variable//2
		else:
			array[x_variable] = 3*x_variable+1
N_total = 2**16
LinkStores = [None for k in range(N_total)]
list_of_threads = [threading.Thread(target=multithreadedRuns,args=(LinkStores,k,)) for k in range((N_total-1)//3+1)]
for k in range((N_total-1)//3+1): list_of_threads[k].start()
for k in range((N_total-1)//3+1): list_of_threads[k].join()
for k in range((N_total-1)//3-1):
	link_runner = k
	while link_runner > 1:
		print(' -> '+str(link_runner),end='')
		next_link = LinkStores[link_runner]
		if next_link is None: break
		elif next_link == 0: break
		else:
			link_runner = LinkStores[link_runner]
	print(' -> ', link_runner)
