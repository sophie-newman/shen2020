from data import *
import numpy as np
from lf_shape import *
import matplotlib
import matplotlib.pyplot as plt
from scipy import ndimage
matplotlib.style.use('classic')
matplotlib.rc('xtick.major', size=15, width=3)
matplotlib.rc('xtick.minor', size=7.5, width=3)
matplotlib.rc('ytick.major', size=15, width=3)
matplotlib.rc('ytick.minor', size=7.5, width=3)
matplotlib.rc('lines',linewidth=4)
matplotlib.rc('axes', linewidth=4)

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

#redshifts = np.array([2,3,4,5,6,7])
redshifts = np.array([2,3,4,5,6])
cmaps = ("Oranges","Reds","Blues","Greens","Purples","Greys")

xmin, xmax = -6.5, -4.5
ymin, ymax = 12, 13.5

gamma1s = np.zeros((len(redshifts), samples.shape[0])) 
gamma2s = np.zeros((len(redshifts), samples.shape[0]))

fig=plt.figure(figsize = (15,10))
ax = fig.add_axes([0.11,0.12,0.79,0.83])

for i in range(len(redshifts)):
	for j in range(samples.shape[0]):
		_,_ ,gamma1s[i,j], gamma2s[i,j] = return_par_at_z(redshifts[i], samples[j,:])

	pixel = 1000
	image,_,_ = np.histogram2d(gamma1s[i,:], gamma2s[i,:] ,bins=[np.linspace(xmin,xmax,pixel+1), np.linspace(ymin,ymax,pixel+1)])

	image[np.invert(np.isfinite(image))] = 0
	image = ndimage.gaussian_filter(image, sigma=2+np.pi/10.)
	image = np.log10(image)
	mini= np.min(image[np.isfinite(image)])
	maxi= np.max(image[np.isfinite(image)])
	medi= np.median(image[np.isfinite(image)])
	image[np.invert(np.isfinite(image))] = np.nan
	x,y = np.meshgrid( np.linspace(xmin,xmax,pixel+1), np.linspace(ymin,ymax,pixel+1) )
	'''
	x = np.linspace(0,2,pixel+1)
	x = (x[:-1]+x[1:])/2.
	y = np.linspace(1,3,pixel+1)
	y = (y[:-1]+y[1:])/2.
	x,y = np.meshgrid(x,y)
	'''
	
	cmap = plt.get_cmap(cmaps[i])
	pos = ax.pcolormesh(x, y, np.transpose(image), cmap=cmap,
        norm=matplotlib.colors.Normalize(vmin=mini,vmax=maxi),alpha=1)
	#ax.contour(x,y,np.transpose(image), (0,), colors='g')

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=2,ncol=1,frameon=False)
ax.set_xlabel(r'$\log{(\phi_{\ast})}$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(L_{\ast})}$',fontsize=40,labelpad=5)

ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.show()
#plt.savefig("../figs/corner1.png",fmt='png')



