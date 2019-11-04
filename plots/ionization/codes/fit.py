from scipy.optimize import curve_fit
import numpy as np

def function(x,eps0,a,b,c,d):
	return eps0 + np.log10( (1+x)**a * np.exp(-b*x)/(np.exp(c*x)+d))

data = np.genfromtxt("emissivity.dat")
x,y = data[:,0],data[:,1]

args,cov = curve_fit(function, x, y, p0=np.array([24., 6,  1., 1., 25.]))#,
	#bounds=([20,0,0,0,20],[30,np.inf,np.inf,np.inf,30]))
print args

print cov

import matplotlib.pyplot as plt
plt.plot(x,y)

xf = np.linspace(0,7,100)
plt.plot(xf, function(xf,*args))
plt.show()

