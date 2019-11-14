import matplotlib.pyplot as plt 
import numpy as np 

f= np.genfromtxt("./FG19_mfp.dat")

z= np.linspace(0.1,7,100)

def func(z):
	return 50*((1+z)/(1+3.5))**(-(2.5+1.94))


plt.plot(f[:,0], f[:,1], 'b')
plt.plot(z, func(z), 'r')
#plt.xscale('log')
plt.yscale('log')
plt.xlim(0,7)
plt.show()