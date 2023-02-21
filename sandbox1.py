from AlgorithmSkeleton1 import *

BurgerPDE = parse_expr('diff(f(x,t),t)+f(x,t)*diff(f(x,t),x)-sigma*diff(f(x,t),x,x)')
print()
print(
	'Let\'s check the Backlund transform from the '+
	'Painleve property of the Burger\'s equation'
	)
alpha, transform, compatibility = ultimatePainleve(
		BurgerPDE, # Input PDE
		x, # The input variable used for the spatial variable
		t # The input variable used for the temporal variable
		# alpha # The power balancing exponent
	)
print()
if alpha == None:
	print('The system is certainly non-integrable!')
else:
	print('The power balancing exponent is ',alpha)
	print('The transform is given by U =',  transform)
	print('where the additional conditions are given by')
	for equation in compatibility:
		if equation != 0: print('>',Eq(equation, 0))
print()

KdVPDE = parse_expr('diff(f(x,t),t)+f(x,t)*diff(f(x,t),x)+sigma*diff(f(x,t),x,x,x)')
print()
print(
	'Let\'s check the Backlund transform from the '+
	'Painleve property of the KdV equation'
	)
alpha, transform, compatibility = ultimatePainleve(
		KdVPDE, # Input PDE
		x, # The input variable used for the spatial variable
		t # The input variable used for the temporal variable
		# alpha # The power balancing exponent
	)
print()
if alpha == None:
	print('The system is certainly non-integrable!')
else:
	print('The power balancing exponent is ',alpha)
	print('The transform is given by U =',  transform)
	print('where the additional conditions are given by')
	for equation in compatibility:
		if equation != 0: print('>',Eq(equation, 0))
print()


TransformSineGordon = parse_expr('2*f(x,t)*diff(f(x,t),x,t)'+
	'-2*diff(f(x,t),x)*diff(f(x,t),t)-f(x,t)**3+f(x,t)')
print()
print(
	'Let\'s check the Backlund transform from the '+\
	'Painleve property of the sine-Gordon equation'+\
	' (with the transformation V = e^\\{i u\\})'
	)
alpha, transform, compatibility = ultimatePainleve(
		TransformSineGordon, # Input PDE
		x, # The input variable used for the spatial variable
		t # The input variable used for the temporal variable
		# alpha # The power balancing exponent
	)
print()
if alpha == None:
	print('The system is certainly non-integrable!')
else:
	print('The power balancing exponent is ',alpha)
	print('The transform is given by U =',  transform)
	print('where the additional conditions are given by')
	for equation in compatibility:
		if equation != 0: print('>',Eq(equation, 0))
print()


BoussinesqEquation = parse_expr('diff(f(x,t),t,t)+2*f(x,t)*diff(f(x,t),x,x)+2*diff(f(x,t),x)**2+1/3*diff(f(x,t),x,x,x,x)')
print()
print(
	'Let\'s check the Backlund transform from the '+\
	'Painleve property of the Boussinesq equation'
	)
alpha, transform, compatibility = ultimatePainleve(
		BoussinesqEquation, # Input PDE
		x, # The input variable used for the spatial variable
		t # The input variable used for the temporal variable
		# alpha # The power balancing exponent
	)
print()
if alpha == None:
	print('The system is certainly non-integrable!')
else:
	print('The power balancing exponent is ',alpha)
	print('The transform is given by U =',  transform)
	print('where the additional conditions are given by')
	for equation in compatibility:
		if equation != 0: print('>',Eq(equation, 0))
print()


ModifiedKdVPDE = parse_expr('diff(f(x, t), t)-diff(f(x, t)**3'+
	'-2*sigma**2*diff(f(x, t),x,x),x)')
print()
print(
	'Let\'s check the Backlund transform from the '+
	'Painleve property of the modified KdV equation'
	)
alpha, transform, compatibility = ultimatePainleve(
		ModifiedKdVPDE, # Input PDE
		x, # The input variable used for the spatial variable
		t # The input variable used for the temporal variable
		# alpha # The power balancing exponent
	)
print()
if alpha == None:
	print('The system is certainly non-integrable!')
else:
	print('The power balancing exponent is ',alpha)
	print('The transform is given by U =',  transform)
	print('where the additional conditions are given by')
	for equation in compatibility:
		if equation != 0: print('>',Eq(equation, 0))
print()


PowerLawNonlinearKG = parse_expr('diff(f(x,t),x,t)-sigma*f(x,t)**4')
print()
print(
	'Let\'s check the Backlund transform from the '+\
	'Painleve property of the power law nonlinear Klein-Gordon'
	)
alpha, transform, compatibility = ultimatePainleve(
		PowerLawNonlinearKG, # Input PDE
		x, # The input variable used for the spatial variable
		t # The input variable used for the temporal variable
		# alpha # The power balancing exponent
	)
print()
if alpha == None:
	print('The system is certainly non-integrable!')
else:
	print('The power balancing exponent is ',alpha)
	print('The transform is given by U =',  transform)
	print('where the additional conditions are given by')
	for equation in compatibility:
		if equation != 0: print('>',Eq(equation, 0))
print()


BBM_PDE = parse_expr('diff(f(x, t), t)+diff(f(x, t), x)'
	+'+f(x, t)*diff(f(x, t), x)-diff(f(x, t), x, x, t)')
print()
print(
	'Let\'s check the Backlund transform from the '+
	'Painleve property of the Benjamin-Bona-Mahoney equation'
	)
alpha, transform, compatibility = ultimatePainleve(
		BBM_PDE, # Input PDE
		x, # The input variable used for the spatial variable
		t # The input variable used for the temporal variable
		# alpha # The power balancing exponent
	)
print()
if alpha == None:
	print('The system is certainly non-integrable!')
else:
	print('The power balancing exponent is ',alpha)
	print('The transform is given by U =',  transform)
	print('where the additional conditions are given by')
	for equation in compatibility:
		if equation != 0: print('>',Eq(equation, 0))


BBM_potential = parse_expr('diff(f(x,t),t)+diff(f(x,t),x)+1/2*diff(f(x,t),x)**2-diff(f(x,t),x,x,t)')
print()
print(
	'Let\'s check the Backlund transform from the '+\
	'Painleve property of the Benjamin-Bona-Mahoney equation'+\
	' in the potential form'
	)
alpha, transform, compatibility = ultimatePainleve(
		BBM_potential, # Input PDE
		x, # The input variable used for the spatial variable
		t # The input variable used for the temporal variable
		# alpha # The power balancing exponent
	)
print()
if alpha == None:
	print('The system is certainly non-integrable!')
else:
	print('The power balancing exponent is ',alpha)
	print('The transform is given by U =',  transform)
	print('where the additional conditions are given by')
	for equation in compatibility:
		if equation != 0: print('>',Eq(equation, 0))
print()

KuramotoSivashinsky = parse_expr('diff(f(x,t),t)+1/2*diff(f(x,t),x)**2+diff(f(x,t),x,x)+diff(f(x,t),x,x,x,x)')
print()
print(
	'Let\'s check the Backlund transform from the '+\
	'Painleve property of the Kuramoto-Sivashinsky equation'+\
	' in the potential form'
	)
alpha, transform, compatibility = ultimatePainleve(
		KuramotoSivashinsky, # Input PDE
		x, # The input variable used for the spatial variable
		t # The input variable used for the temporal variable
		# alpha # The power balancing exponent
	)
print()
if alpha == None:
	print('The system is certainly non-integrable!')
else:
	print('The power balancing exponent is ',alpha)
	print('The transform is given by U =',  transform)
	print('where the additional conditions are given by')
	for equation in compatibility:
		if equation != 0: print('>',Eq(equation, 0))


# higherPowerGenKdV = parse_expr('diff(f(x,t),t)+diff(f(x,t),x,x,x)+f(x,t)**14*diff(f(x,t),x)')
# print()
# print(
# 	'Let\'s check the Backlund transform from the '+\
# 	'Painleve property of the power-14 generalized KdV'
# 	)
# alpha, transform, compatibility = ultimatePainleve(
# 		higherPowerGenKdV, # Input PDE
# 		x, # The input variable used for the spatial variable
# 		t # The input variable used for the temporal variable
# 		# alpha # The power balancing exponent
# 	)
# print()
# if alpha == None:
# 	print('The system is certainly non-integrable!')
# else:
# 	print('The power balancing exponent is ',alpha)
# 	print('The transform is given by U =',  transform)
# 	print('where the additional conditions are given by')
# 	for equation in compatibility:
# 		if equation != 0: print('>',Eq(equation, 0))
# print()


nonlinearPseudoParabolicDiffusion = parse_expr('(diff(f(x,t),t)-diff(f(x,t),x,x,t))**2-9/4*sigma**2*diff(f(x,t),x)*diff(f(x,t),x,x)**2')
print()
print(
	'Let\'s check the Backlund transform from the '+\
	'Painleve property of the nonlinear pseudo-parabolic type diffusion equation'
	)
alpha, transform, compatibility = ultimatePainleve(
		nonlinearPseudoParabolicDiffusion, # Input PDE
		x, # The input variable used for the spatial variable
		t # The input variable used for the temporal variable
		# alpha # The power balancing exponent
	)
print()
if alpha == None:
	print('The system is certainly non-integrable!')
else:
	print('The power balancing exponent is ',alpha)
	print('The transform is given by U =',  transform)
	print('where the additional conditions are given by')
	for equation in compatibility:
		if equation != 0: print('>',Eq(equation, 0))
print()


# b = 2
# CamassaHolmEquation = parse_expr('2*f(x,t)**3*diff(f(x,t),x)'
# 	+'+diff(f(x,t),t)*(diff(f(x,t),x)**2-f(x,t)*diff(f(x,t),x,x)-1)'
# 	+'+f(x,t)*(f(x,t)*diff(f(x,t),x,x,t)-diff(f(x,t),x)*diff(f(x,t),x,t))'
# 	)
# print()
# print(
# 	'Let\'s check the Backlund transform from the '+
# 	'Painleve property of the Camassa-Holm equation'
# 	)
# alpha, transform, compatibility = ultimatePainleve(
# 		CamassaHolmEquation, # Input PDE
# 		x, # The input variable used for the spatial variable
# 		t, # The input variable used for the temporal variable
# 		# alpha # The power balancing exponent
# 	)
# print()
# print('The power balancing exponent is ',alpha)
# print('The transform is given by U = ',  transform)
# print('where the additional conditions are given by')
# for equation in compatibility:
# 	if equation != 0: print('>',Eq(equation, 0))


# b = 3
# DegasperisProcesiEquation = parse_expr('3*f(x,t)**4*diff(f(x,t),x)'
# 	+'+diff(f(x,t),t)*(diff(f(x,t),x)**2-f(x,t)*diff(f(x,t),x,x)-1)'
# 	+'+f(x,t)*(f(x,t)*diff(f(x,t),x,x,t)-diff(f(x,t),x)*diff(f(x,t),x,t))'
# 	)
# print()
# print(
# 	'Let\'s check the Backlund transform from the '+
# 	'Painleve property of the Degasperis-Procesi equation'
# 	)
# alpha, transform, compatibility = ultimatePainleve(
# 		DegasperisProcesiEquation, # Input PDE
# 		x, # The input variable used for the spatial variable
# 		t, # The input variable used for the temporal variable
# 		# alpha # The power balancing exponent
# 	)
# print()
# print('The power balancing exponent is ',alpha)
# print('The transform is given by U = ',  transform)
# print('where the additional conditions are given by')
# for equation in compatibility:
# 	if equation != 0: print('>',Eq(equation, 0))

# b = 2.5
# alpha = -2
# bFamilyEquation = parse_expr('2.5*f(x,t)**(3.5)*diff(f(x,t),x)'
# 	+'+diff(f(x,t),t)*(diff(f(x,t),x)**2-f(x,t)*diff(f(x,t),x,x)-1)'
# 	+'+f(x,t)*(f(x,t)*diff(f(x,t),x,x,t)-diff(f(x,t),x)*diff(f(x,t),x,t))'
# 	)
# print()
# print(
# 	'Let\'s check the Backlund transform from the '+
# 	'Painleve property of the b-Family equation with b = '+str(b)
# 	)
# transform, compatibility = PainlevePDE_withBacklundTransform(
# 		bFamilyEquation, # Input PDE
# 		x, # The input variable used for the spatial variable
# 		t, # The input variable used for the temporal variable
# 		alpha # The power balancing exponent
# 	)
# print('The transform is given by U = ',  transform)
# print('where the additional conditions are given by')
# for equation in compatibility:
# 	if equation != 0: print('>',Eq(equation, 0))
