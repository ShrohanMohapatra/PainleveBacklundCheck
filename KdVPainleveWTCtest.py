from sympy import *
u, u0, u1, u2, phi = symbols('u u0 u1 u2 phi', cls = Function)
x, t, sigma = symbols('x t sigma')
KdV_Equation = diff(u(x, t), t)+u(x, t)*diff(u(x, t), x)+\
				sigma*diff(u(x, t), (x, 3))
print('----------------------------------------------------------------')
print('KdV equation is as follows .....')
print('----------------------------------------------------------------')
print(KdV_Equation)
print('----------------------------------------------------------------')
print('Painlev\'e expansion')
print('----------------------------------------------------------------')
PainleveExpansion = u0(x, t)/phi(x, t)**2 + u1(x, t)/phi(x, t) + u2(x, t)
WTCSubstitution = phi(x, t)**5*KdV_Equation.subs(
	u(x, t), PainleveExpansion)
WTCSubstitution = expand(WTCSubstitution.doit())
print(WTCSubstitution)
print('----------------------------------------------------------------')
print('Coefficients')
print('----------------------------------------------------------------')
coeff_list = [WTCSubstitution.coeff(phi(x, t), k) for k in range(1, 6)]
sum1 = 0
for k in range(1, 6):
	sum1 = sum1 + coeff_list[k-1] * phi(x, t)**k
coeff_list.insert(0, simplify(WTCSubstitution-sum1))
for k in range(len(coeff_list)):
	print('k = ',k,' Coefficient = ', coeff_list[k])
print('----------------------------------------------------------------')
u0_1 = solve(Eq(coeff_list[0], 0), u0(x, t))[1]
u1_1 = solve(Eq(coeff_list[1].subs(u0(x, t), u0_1), 0), u1(x, t))[0]
u2_1 = solve(
	Eq(coeff_list[2].subs(u0(x, t), u0_1).\
		subs(u1(x, t), u1_1), 0), u2(x, t))[0]
print('----------------------------------------------------------------')
print('u0(x, t) = ',u0_1)
print('----------------------------------------------------------------')
print('u1(x, t) = ',u1_1)
print('----------------------------------------------------------------')
print('u2(x, t) = ',u2_1)
print('----------------------------------------------------------------')
init_printing()
u2_2 = expand(coeff_list[2].subs(u0(x, t), u0_1).subs(u1(x, t), u1_1))
u2_2, _ = fraction(together(u2_2.doit().as_independent([x, t, sigma])[0]))
u2_2 = expand(u2_2.expand()/(24*sigma*diff(phi(x, t), x)))
print('Equation 1')
print('----------------------------------------------------------------')
print(u2_2)
print('----------------------------------------------------------------')
print('----------------------------------------------------------------')
init_printing()
u3 = expand(coeff_list[3].subs(u0(x, t), u0_1).subs(u1(x, t), u1_1))
u3, _ = fraction(together(u3.doit().as_independent([x, t, sigma])[0]))
u3 = expand(u3.expand()/(-12*sigma))
print('Equation 2')
print('----------------------------------------------------------------')
print(u3)
print('----------------------------------------------------------------')
print('----------------------------------------------------------------')
print('Equation 3')
print('----------------------------------------------------------------')
u4 = expand(coeff_list[4].subs(u0(x, t), u0_1).subs(u1(x, t), u1_1))
u4, _ = fraction(together(u4.doit().as_independent([x, t, sigma])[0]))
u4 = expand(u4.expand()/(12*sigma))
print(u4)
print('----------------------------------------------------------------')
print('Equation 4')
print('----------------------------------------------------------------')
u5 = expand(coeff_list[5].subs(u0(x, t), u0_1).subs(u1(x, t), u1_1))
u5, _ = fraction(together(u5.doit().as_independent([x, t, sigma])[0]))
u5 = u5.expand()
print(u5)
print('----------------------------------------------------------------')

print('----------------------------------------------------------------')
print('1/diff(phi(x,t),x)*diff(Equation_2 - diff(Equation_1, x), x)'+
	' - Equation_3')
print('----------------------------------------------------------------')
test_1 = expand(diff(1/diff(phi(x,t), x)*(u3 - diff(u2, x)), x) - u4)
print(test_1)
print('----------------------------------------------------------------')

print('----------------------------------------------------------------')
equation3 = expand(coeff_list[3].subs(u0(x, t), u0_1).subs(u1(x, t), u1_1).subs(u2(x, t), u2_1))
equation3, _ = fraction(together(equation3.doit().as_independent([x, t, sigma])[0]))
equation3 = expand(equation3/(-12*sigma))
print('Equation number 3')
print('----------------------------------------------------------------')
print(equation3)
print('----------------------------------------------------------------')

# operator3 = lambdify([x, t, phi], equation3)

print('----------------------------------------------------------------')
equation4 = expand(coeff_list[4].subs(u0(x, t), u0_1).subs(u1(x, t), u1_1).subs(u2(x, t), u2_1))
equation4, _ = fraction(together(equation4.doit().as_independent([x, t, sigma])[0]))
equation4 = expand(equation4/(-12*sigma))
print('Equation number 4')
print('----------------------------------------------------------------')
print(equation4)
print('----------------------------------------------------------------')

# operator4 = lambdify((x, t, phi), equation4)

print('----------------------------------------------------------------')
equation5 = expand(coeff_list[5].subs(u0(x, t), u0_1).subs(u1(x, t), u1_1).subs(u2(x, t), u2_1))
equation5, _ = fraction(together(equation5.doit().as_independent([x, t, sigma])[0]))
equation5 = expand(equation5)
print('Equation number 5')
print('----------------------------------------------------------------')
print(equation5)
print('----------------------------------------------------------------')

print('----------------------------------------------------------------')
coefficient3 = equation3.coeff(Derivative(phi(x, t), (x, 4)), 1)
coefficient4 = equation4.coeff(Derivative(phi(x, t), (x, 5)), 1)
equation3_new = equation3/coefficient3 - Derivative(phi(x, t), (x, 4))
equation3_new, _ = fraction(together(equation3_new.doit().as_independent([x, t, sigma])[0]))
equation4_new = equation4/coefficient4 - Derivative(phi(x, t), (x, 5))
equation4_new, _ = fraction(together(equation4_new.doit().as_independent([x, t, sigma])[0]))
equation4_next = diff(equation3, (x, 5)) - diff(equation4, (x, 4))
equation4_next, _ = fraction(together(equation4_next.doit().as_independent([x, t, sigma])[0]))
print('----------------------------------------------------------------')
print('equation4_next = ')
print('----------------------------------------------------------------')
print(equation4_next)
print('----------------------------------------------------------------')

# print('----------------------------------------------------------------')
# commutator34 = expand(equation3.doit().subs({phi(x,t): equation4})-\
# 	equation4.doit().subs({phi(x,t): equation3}))
# print('Commutator attempt')
# print('----------------------------------------------------------------')
# print(commutator34)
# print('----------------------------------------------------------------')

# print('----------------------------------------------------------------')
# commutator45 = expand(equation4.doit().subs({phi(x,t): equation5})-\
# 	equation5.doit().subs({phi(x,t): equation4}))
# print('Commutator attempt')
# print('----------------------------------------------------------------')
# print(commutator45)
# print('----------------------------------------------------------------')

# print('----------------------------------------------------------------')
# commutator35 = expand(equation3.doit().subs({phi(x,t): equation5})-\
# 	equation5.doit().subs({phi(x,t): equation3}))
# print('Commutator attempt')
# print('----------------------------------------------------------------')
# print(commutator35)
# print('----------------------------------------------------------------')

# print('----------------------------------------------------------------')
# commutator55 = expand(equation5.doit().subs({phi(x,t): equation5})-\
# 	equation5.doit().subs({phi(x,t): equation5}))
# print('Commutator attempt')
# print('----------------------------------------------------------------')
# print(commutator55)
# print('----------------------------------------------------------------')

# operator5 = lambdify((x, t, phi), equation5)

# print(operator3(x, t, phi(x, t)))

