import numpy as np
import matplotlib.pyplot as plt

def f1(chord, x):
	res = chord * ((x/chord)**(0.5)*(1-x/chord))/(np.exp(15*x/chord))
	return res

def f2(chord, x):
	res = chord * np.sin(np.pi*(x/chord)**0.25)**3
	return res

def f3(chord, x):
	res = chord * np.sin(np.pi*(x/chord)**0.757)**3
	return res

def f4(chord, x):
	res = chord * np.sin(np.pi*(x/chord)**1.357)**3
	return res

def defPoints(points,chord,coeffUp,coeffLo):
	x = points[:,0]
	y = points[:,1]
	for i in range(len(x)):
			if y[i] > 0:
				for j in range(len(coeffUp)):
					y[i] += eval("coeffUp["+str(j)+"]*f"+str(j+1)+"(chord,x["+str(i)+"])")
			else:
				for j in range(len(coeffLo)):
					y[i] += eval("coeffLo["+str(j)+"]*f"+str(j+1)+"(chord,x["+str(i)+"])")
	points_mod = points
	points_mod[:,1] = y
	return points_mod
	






chord = 1
x0 = -0.5

a = np.loadtxt("points")
a[:,0] = a[:,0]-x0
xtest = 0 #a[0,1]
coeffUp=[0.01,0.01,0.01,0.01]*0
coeffLo=[0,0,0,0]*0
pmod = defPoints(a,chord,coeffUp,coeffLo)
plt.plot(pmod[:,0],pmod[:,1],"*")
plt.plot([-1,1],[0,0])
plt.show()
print(f1(chord,xtest)+f2(chord,xtest)+f4(chord,xtest))