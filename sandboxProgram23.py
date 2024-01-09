def increase_by_one(a, b, c):
	return [a+1, b+1, c+1]
a, b, c = 1, 2, 3
print('a = ', a, 'b = ', b, 'c = ', c)
[a, b, c] = increase_by_one(a, b, c)
print('a = ', a, 'b = ', b, 'c = ', c)