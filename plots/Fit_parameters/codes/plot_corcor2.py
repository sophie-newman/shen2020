from data import *
import numpy as np
from lf_shape import *
import matplotlib
import matplotlib.pyplot as plt
import corner
from scipy import ndimage
matplotlib.style.use('classic')
matplotlib.rc('xtick.major', size=10, width=2)
matplotlib.rc('xtick.minor', size=5, width=2)
matplotlib.rc('ytick.major', size=10, width=2)
matplotlib.rc('ytick.minor', size=5, width=2)
matplotlib.rc('lines',linewidth=2)
matplotlib.rc('axes', linewidth=2)
'''
def return_par_at_z(redshift,p):
	zref = 2.
	gamma1 = polynomial(redshift,(p[0],p[1],p[2]),2)
	gamma2 = doublepower(redshift,(p[3],zref, p[4], p[5]))
	logphi = polynomial(redshift,(p[6],p[7]),1)
	Lbreak = doublepower(redshift,(p[8],zref, p[9], p[10]))
	return np.array([gamma1,gamma2,logphi,Lbreak])

data = np.genfromtxt("../../../codes/lf_fit/output/chain.dat")
nburn = 1000
nwalker = 100
samples = data[nburn*nwalker+1:,2:]

redshift = 5
par_dpl = np.zeros((samples.shape[0],4)) 

print samples.shape[0]
for j in range(samples.shape[0]):
	if (j % 1000)==0: print j
	par_dpl[j,:] = return_par_at_z(redshift, samples[j,:])

np.save("chain", par_dpl)

'''
redshift = 5.
source = np.genfromtxt("../../Fit_parameters/codes/zevolution_fit_global.dat",names=True)
zref = 2.
p=source['value'][ source['paraid']==0 ]
gamma1 = polynomial(redshift,p,2)
p=source['value'][ source['paraid']==1 ]
gamma2 = doublepower(redshift,(p[0],zref, p[1], p[2]))
p=source['value'][ source['paraid']==2 ]
logphi = polynomial(redshift,p,1) 
p=source['value'][ source['paraid']==3 ]
Lbreak = doublepower(redshift,(p[0],zref, p[1], p[2]))
par = np.array([gamma1,gamma2,logphi,Lbreak])

par = (par[0],par[1],par[2],par[3]-12)

samples = np.load("chain.npy")
samples[:,-1] = samples[:,-1] - 12

dict={
	"fontsize":20,
	"labelpad":2.5
}

fig = corner.corner(samples,top_ticks=True,color='k', range=[(0.6,0.9),(1.8,2.05),(-6.1,-5.8),(0.6,0.9)],
	labels=(r"$\gamma_1$", r"$\gamma_2$", r"$\log{(\phi_{\ast})}$", r"$\log{(L_{\ast})} - 12$"), label_kwargs=dict,
	plot_contours=True, quantiles=(0.16,0.84), truths=par)

#plt.show()
plt.savefig("../figs/corner2.png",fmt='png')



