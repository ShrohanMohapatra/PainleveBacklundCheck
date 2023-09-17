from sympy import *
u0, u1, u2, phi = symbols('u0 u1 u2 phi', cls = Function)
x, t, sigma = symbols('x t sigma')
u = u0(x, t)/phi(x, t)**2 + u1(x, t)/phi(x, t) + u2(x, t)
bbm_equation = phi(x, t)**5*(diff(u, t) + diff(u, x) + u*diff(u, x) - diff(u, x, x, t))
bbm_equation, _ = fraction(together(bbm_equation.doit()))
bbm_equation = expand(bbm_equation)
print('-------------------------------')
print('bbm_equation')
print('-------------------------------')
print(bbm_equation)
print('-------------------------------')
bbm_equation_coeffs = [bbm_equation.coeff(phi(x, t), k+1) for k in range(6)]
for k in range(len(bbm_equation_coeffs)):
	print('-------------------------------')
	print('bbm_equation_coeffs['+str(k)+']')
	print('-------------------------------')
	print(bbm_equation_coeffs[k])
	print('-------------------------------')
bbm_equation_with_all_phi_powers = 0
for k in range(6):
	bbm_equation_with_all_phi_powers = bbm_equation_with_all_phi_powers + bbm_equation_coeffs[k] * phi(x, t)**(k+1)
bbm_equation_coeffs.insert(0, simplify(bbm_equation - bbm_equation_with_all_phi_powers))
print('-------------------------------')
print('New coefficient lists')
print('-------------------------------')
for k in range(len(bbm_equation_coeffs)):
	print('-------------------------------')
	print('bbm_equation_coeffs['+str(k)+']')
	print('-------------------------------')
	print(bbm_equation_coeffs[k])
	print('-------------------------------')
