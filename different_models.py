from matplotlib import pyplot as plt
import numpy as np
import math
import statistics

start = input('enter the starting keyword:: \n')

#tumor grwoth model Gomparison relation
def tgm(l, a ,v0, t_start, t_final, h):
	l = float(l); 	a = float(a); 	v0 = float(v0); 	t_start = float(t_start); 	t_final =  float(t_final); 	h = float(h)
	g = []
	y = []
	while t0 < tmax:
		k1 = f(t0 , v0)
		k2 = f(t0 + h/2.0 , v0 + h*k1*0.5)
		k3 = f(t0 + h/2.0 , v0 + h*k2*0.5)
		k4 = f(t0 + h , v0 + h*k3)
	
		v = v0	+ (k1 + 2*(k2 + k3) + k4)*h/6.0
		g.append(v)
		v0 = v
	
		t0 = t0 + h
		y.append(t0)

	plt.plot(y, g, 'b--')
	plt.show()

#non linear prey predator model for single species
def nlppm(c0, n, r, k):
	c = []
	for i in N:
		x = c0*(1 + r*(1 - c0 /k))
		c.append(x)
		c0 = x

	plt.plot(c)
	plt.show()

#non linear prey predator model for 2 species xn+1 = xn(1+a) - xnyn*b
def nlppm2(n, g, d, b, i1, i2):
	pre0 = np.random.random()
	prd0 = np.random.random()
	prey = np.zeros(n)
	predator = np.zeros(n)
	prey[0] = pre0
	predator[0] = prd0

	for i in range(1,n):
		prey[i] = prey[i-1]*(1 + g) - i1*prey[i-1]*predator[i-1] - b*prey[i-1]**2
		predator[i] = predator[i-1] * (1 - d)  + i2*predator[i-1]*prey[i-1]

	print(prey, predator)
	plt.plot(prey, 'k--', label='prey')
	plt.plot(predator, 'b--' , label = 'predator')
	plt.legend()
	plt.show()


if start == 'nlppm2':
	n = int(input('enter the value of n::'))
	g = float(input('enter the grwoth rate::'))
	d = float(input('enter the death rate::'))
	b = float(input('enter the value of beta::'))
	i1 = float(input('enter the interaction value for prey::'))
	i2 = float(input('enter the interaction value for predator::'))
	nlppm2(n, g, d,b, i1, i2)

## plot of Pearl-Verhulst logistic Equation
def pvle(r, n):
	x0 = np.random.random()*0.1
	x = np.zeros(n)
	x[0] = x0
	for i in range(1 , n):
		x[i] = r*x[i-1]*(1 - x[i-1])

	plt.plot(x , 'b.--')
	plt.xlabel('index')
	plt.ylabel('xn')
	plt.show()

###initiating parameters
if start == 'pvle':
	r = float(input('enter the value of r::\n'))
	n = int(input('enter the value of n::\n'))
	pvle(r,n)

##continuous verhulst logistic equation
def clem(p0 , r1 , k, t0, tmax, dt):
	p = []
	p.append(p0)
	t = []
	t.append(t0)
	dp_dt = lambda p: r1*p*(1 - p/k)
	
	while t0 < tmax:
		p1 = p0 + dp_dt(p0) * dt
		p.append(p1)
		p0 = p1

		t0 = t0 + dt
		t.append(t0)

	plt.plot(t, p, 'b-.')
	plt.xlabel('time ---->')
	plt.ylabel('populaton ---->')
	plt.show()

####initiating parameters 
if start == 'clem':
	p0 = float(input('enter the initial population:\n'))
	r1 = float(input('enter the value of r1:\n'))
	k = float(input('ente the value of k:\n'))
	t0 = float(input('enter the start time:\n'))
	tmax = float(input('enter the stopping time:\n'))
	dt = float(input('enter the step size for time:\n'))
	clem(p0, r1, k, t0, tmax, dt)
#birfurcation
	#saddle node bifurcation
	#Trans critical bifurcation
	#pitchfork bifurcation

#chemostat model
#using the dimensionless parameters  #coupled equations
def chemostat(F, V,  C1, C0, N0, N1, t1 , tmax, Kn, Kmax, a, h):
	Numb = []
	Conc = []
	t = []
	
	Numb.append(N0)
	Conc.append(C0)
	
	dN_dt = lambda C, N : t1*Kmax * (C/(Kn/C1 + C)) - t1*F * N/V 	#varitation of N*(number of microbes)
	dC_dt = lambda C, N: (-a*Kmax/C1) * t1 * C/((Kn/C1) + C*C1)*N*N1 - t1*F*C/V + t1*F*C0/(V*C1)
	
	t0 = 0
	t.append(t0)

	while t0 < tmax:
		k1 = dN_dt(C0 , N0)
		L1 = dC_dt(C0 , N0)
		k2 = dN_dt(C0 + 0.5*h, N0 + k1*h*0.5)
		L2 = dC_dt(C0 + 0.5*h , N0 + L1 * h * 0.5)
		k3 = dN_dt(C0 + 0.5*h , N0 + k2*h*0.5)
		L3 = dC_dt(C0 + 0.5*h , N0 + L2*h*0.5)
		k4 = dN_dt(C0 + h , N0 + k3*h)
		L4 = dC_dt(C0 + h , N0 + L3*h)

		C0 = C0 + (k1 + 2*(k2 + k3) + k4)/6.0
		N0 = N0 + (L1 + 2*(L2 + L3) + L4)/6.0
		Conc.append(C0)
		Numb.append(N0)

		t0 = t0 + h
		t.append(t0)

	plt.plot(t , C0 , 'b--')
	plt.show()

### starting parameters

if start == 'chemostat':
	N0 = float(input('enter the inital number of microbes::'))
	C0 = float(input('enter the inital concentration::'))
	N1 = float(input('enter the value for const N1::'))
	C1 = float(input('enter the value of const C1::'))
	t1 = float(input('enter the value of const t1::'))
	tmax = float(input('enter the maximum value of time::'))
	Kn = float(input('enter the value of Kn::'))
	Kmax = float(input('enter the value of value Kmax::'))
	a = float(input('enter the value for a ::'))
	h = float(input('enter the step size::'))
	F = float(input('enter the influx rate::'))
	V = float(input('enter the value of Volume of chemostat::'))
	chemostat(F, V, C1, C0, N0, N1, t1, tmax, Kn, Kmax, a, h)