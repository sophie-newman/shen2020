from data import *
import numpy as np 
import astropy.constants as con
import matplotlib.pyplot as plt 
import matplotlib
matplotlib.style.use('classic')
matplotlib.rc('xtick.major', size=15, width=3)
matplotlib.rc('xtick.minor', size=7.5, width=3)
matplotlib.rc('ytick.major', size=15, width=3)
matplotlib.rc('ytick.minor', size=7.5, width=3)
matplotlib.rc('lines',linewidth=4)
matplotlib.rc('axes', linewidth=4)

T0 = np.polynomial.chebyshev.Chebyshev((1,0,0,0))
T1 = np.polynomial.chebyshev.Chebyshev((0,1,0,0))
T2 = np.polynomial.chebyshev.Chebyshev((0,0,1,0))
T3 = np.polynomial.chebyshev.Chebyshev((0,0,0,1))
def polynomial(z,p,n=3):
	xsi=1.+z
	if n==1: return p[0]*T0(xsi)+p[1]*T1(xsi)
	elif n==2: return p[0]*T0(xsi)+p[1]*T1(xsi)+p[2]*T2(xsi)
	elif n==3:
		return p[0]*T0(xsi)+p[1]*T1(xsi)+p[2]*T2(xsi)+p[3]*T3(xsi)
	else: return False

def doublepower(z,p):
	xsi=1.+z
	zref=p[1]
	return 2*p[0]/(np.power(xsi/(1+zref),p[2]) + np.power(xsi/(1+zref),p[3]))

def bestfit(z,field):
	source=np.genfromtxt("zevolution_fit.dat",names=['gamma1','gamma2','phi_s','Lbreak'])
	p=source[field]
	if (field=='gamma1'): 
		return polynomial(z,p,2)
	elif (field=='phi_s'):
		return polynomial(z,p,1)
	else: return doublepower(z,p)

def bestfit_global(z,paraid):
	source=np.genfromtxt("zevolution_fit_global.dat",names=True)
	zref = 2.0
	p=source['value'][ source['paraid']==paraid ]
	print p
	if (paraid==0): 
		return polynomial(z,p,2)
	elif (paraid==2):
		return polynomial(z,p,1)
	else: return doublepower(z,(p[0],zref,p[1],p[2]))

def Hopkins07(z):
	parameters_init = np.array([0.41698725, 2.17443860, -4.82506430, 13.03575300, 0.63150872, -11.76356000, -14.24983300, -0.62298947, 1.45993930, -0.79280099])
	xsi_log	= np.log10((1.+z)/(1.+2.))
	gamma1_0	= parameters_init[0]	#faint-end slope
	gamma2_0	= parameters_init[1] 	#bright-end slope
	P0			= parameters_init[2]	#normalization in log
	L0			= parameters_init[3]	#break in log
	k1, k2, k3 = parameters_init[4], parameters_init[5], parameters_init[6]
	k_gamma1 = parameters_init[7]
	k_gamma2_1 = parameters_init[8]
	k_gamma2_2 = parameters_init[9]
	gamma1   = gamma1_0 * np.power(10., k_gamma1*xsi_log)
	gamma2   = 2.*gamma2_0 / (np.power(10., xsi_log*k_gamma2_1) + np.power(10., xsi_log*k_gamma2_2))
	Lbreak  = L0 + k1*xsi_log + k2*xsi_log**2 + k3*xsi_log**3
	return gamma1,gamma2,P0*np.ones(len(z)),Lbreak

z_a=np.linspace(0.01,8,1000)
gamma1_a, gamma2_a, phi_s_a, Lbreak_a = Hopkins07(z_a)

fig=plt.figure(figsize = (15,10))
ax = fig.add_axes([0.11,0.12,0.79,0.83])

data=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_nofix.dat",names=True)
ax.errorbar(data["z"],data["gamma1"],yerr=data['err1'],linestyle='',marker='o',
	c='gray',mec='gray',ms=18,capsize=10,capthick=4,alpha=0.5,label=r'$\rm Local$ $\rm fits$ ($\rm free$)')

data=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_fix.dat",names=True)
ax.plot(data["z"],data["gamma1"],linestyle='',marker='o',
	c='royalblue',mec='royalblue',ms=18,label=r'$\rm Local$ $\rm fits$ ($\phi_{\ast}$ $\rm fixed$)')

ax.plot(z_a,gamma1_a,'--',dashes=(25,15),c='crimson',label=r'$\rm Hopkins+$ $\rm 2007$')
ax.plot(z_a,bestfit(z_a,'gamma1'),'-',c='seagreen',label=r'$\rm Fit$ $\rm on$ $\rm local$ $\rm fits$')
ax.plot(z_a,bestfit_global(z_a,0),'-',c='darkorchid',label=r'$\rm Global$ $\rm fit$')

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=2,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\rm Faint-end$ $\rm Slope$ $\rm \gamma_{\rm 1}$',fontsize=40,labelpad=5)

ax.set_xlim(0,7.)
ax.set_ylim(-0.1,1.7)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/gamma1.pdf",fmt='pdf')

fig=plt.figure(figsize = (15,10))
ax = fig.add_axes([0.11,0.12,0.79,0.83])

data=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_nofix.dat",names=True)
ax.errorbar(data["z"],data["gamma2"],yerr=data['err2'],linestyle='',marker='o',c='gray',mec='gray',ms=18,capsize=10,capthick=4,alpha=0.5)

data=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_fix.dat",names=True)
ax.plot(data["z"],data["gamma2"],linestyle='',marker='o',c='royalblue',mec='royalblue',ms=18)


ax.plot(z_a,gamma2_a,'--',dashes=(25,15),c='crimson')
ax.plot(z_a,bestfit(z_a,'gamma2'),'-',c='seagreen')
ax.plot(z_a,bestfit_global(z_a,1),'-',c='darkorchid')

ax.axhspan(0,1,color='yellow',alpha=0.5,label=r'$\rm Divergence$')
prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=2,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\rm Bright-end$ $\rm Slope$ $\rm \gamma_{\rm 2}$',fontsize=40,labelpad=5)

ax.set_xlim(0,7.)
ax.set_ylim(0.9,2.9)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/gamma2.pdf",fmt='pdf')

fig=plt.figure(figsize = (15,10))
ax = fig.add_axes([0.11,0.12,0.79,0.83])

data=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_nofix.dat",names=True)
ax.errorbar(data["z"],data["phi_s"],yerr=data['err3'],linestyle='',marker='o',c='gray',mec='gray',ms=18,capsize=10,capthick=4,alpha=0.5)

data=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_fix.dat",names=True)
fixid = data['err3']==0
unfix = data['err3']!=0
ax.plot(data["z"][fixid],data["phi_s"][fixid],linestyle='',marker='o',fillstyle='none', mew=3,c='black',mec='black',ms=18)
ax.plot(data["z"][unfix],data["phi_s"][unfix],linestyle='',marker='o',c='royalblue',mec='royalblue',ms=18)

ax.plot(z_a,phi_s_a,'--',dashes=(25,15),c='crimson')
ax.plot(z_a,bestfit(z_a,'phi_s'),'-',c='seagreen') 
ax.plot(z_a,bestfit_global(z_a,2),'-',c='darkorchid')

#prop = matplotlib.font_manager.FontProperties(size=25.0)
#ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=3,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(\phi_{\ast}\,[{\rm dex}^{-1}\,{\rm Mpc}^{-3}])}$',fontsize=40,labelpad=5)

ax.set_xlim(0,7.)
#ax.set_ylim(-6.3,-4.1)
ax.set_ylim(-6.3,-3.1)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/phi_s.pdf",fmt='pdf')

####################################################
fig=plt.figure(figsize = (15,10))
ax = fig.add_axes([0.11,0.12,0.79,0.83])

#load and plot observational data
####################################################
from lf_fitter_data import *
from convolve import *
from scipy import ndimage
def get_data(alldata,zmin,zmax,dset_name,dset_id):
        z_lis = 0.5*(zmin + zmax)
        for iz in range(len(z_lis)):
                redshift = z_lis[iz]
                L_data, PHI_data, DPHI_data = load_LF_data[dset_name](redshift)

		if (dset_id!=-5):
                	for i in range(len(L_data)):
                        	L_data[i] = bolometric_correction_inverse(L_data[i],dset_id)
        	elif (dset_id==-5):
                	L_1450 = -0.4*L_data + np.log10(Fab*(con.c.value/1450e-10)) - L_solar
			for i in range(len(L_data)):
                        	L_data[i] = bolometric_correction_inverse(L_1450[i],dset_id)
		
		DL_data = 0.0*L_data
		if len(L_data)!=1:
			for j in range(len(L_data)):
				if (j==0): DL_data[j] = (L_data[j+1]-L_data[j])/2.
				elif (j==len(L_data)-1): DL_data[j] = (L_data[j]-L_data[j-1])/2.
				else: DL_data[j] = (L_data[j+1]-L_data[j-1])/4.
		else: DL_data = np.array([0.2])

                alldata["L_OBS"]  = np.append(alldata["L_OBS"]  , L_data)
		alldata["DLOBS"]  = np.append(alldata["DLOBS"]  , DL_data)
                alldata["Z"]  = np.append(alldata["Z"]  , np.ones(len(L_data)) * redshift)
		alldata["Zlo"]= np.append(alldata["Zlo"], np.ones(len(L_data)) * zmin[iz])
		alldata["Zup"]= np.append(alldata["Zup"], np.ones(len(L_data)) * zmax[iz])
                return len(L_data)

def loop_over_all_dset():
        alldata={"L_OBS":np.array([]),"Zlo":np.array([]),"Zup":np.array([]),"Z":np.array([]),"DLOBS":np.array([])}
        for key in dset_ids.keys():
        	get_data(alldata,zmins[key],zmaxs[key],key,dset_ids[key])

	return alldata["Z"], alldata["Zlo"], alldata["Zup"] ,alldata["L_OBS"], alldata["DLOBS"]

pixel = 500
x,xlo,xup,y,dy = loop_over_all_dset()
r1,_,_ = np.histogram2d(x,y,bins=[np.linspace(0,7,pixel+1), np.linspace(10.3,13.8,pixel+1)])
r2,_,_ = np.histogram2d(xup,y,bins=[np.linspace(0,7,pixel+1), np.linspace(10.3,13.8,pixel+1)])
r3,_,_ = np.histogram2d(xlo,y,bins=[np.linspace(0,7,pixel+1), np.linspace(10.3,13.8,pixel+1)])
r4,_,_ = np.histogram2d(x,y+dy,bins=[np.linspace(0,7,pixel+1), np.linspace(10.3,13.8,pixel+1)])
r5,_,_ = np.histogram2d(x,y-dy,bins=[np.linspace(0,7,pixel+1), np.linspace(10.3,13.8,pixel+1)])

image = r1+r2+r3+r4+r5
image[np.invert(np.isfinite(image))] = 0
image = ndimage.gaussian_filter(image, sigma=30)
image = np.log10(image)
mini= np.min(image[np.isfinite(image)])
maxi= np.max(image[np.isfinite(image)])
medi= np.median(image[np.isfinite(image)])
print mini,maxi
image[np.invert(np.isfinite(image))] = mini
x,y = np.meshgrid( np.linspace(0,7,pixel+1), np.linspace(10.3,13.8,pixel+1) )

#cmap = plt.get_cmap('magma_r')
cmap = plt.get_cmap('Greys')
pos = ax.pcolormesh(x, y, np.transpose(image), cmap=cmap, 
	norm=matplotlib.colors.Normalize(vmin=medi-0.2*(medi-mini),vmax=maxi+3.*(maxi-medi)),alpha=1)
####################################################

data=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_nofix.dat",names=True)
ax.errorbar(data["z"],data["L_s"],yerr=data['err4'],linestyle='',marker='o',c='gray',mec='gray',ms=18,capsize=10,capthick=4,alpha=1)

data=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_fix.dat",names=True)
ax.plot(data["z"],data["L_s"],linestyle='',marker='o',c='royalblue',mec='royalblue',ms=18)

ax.plot(z_a,Lbreak_a,'--',dashes=(25,15),c='crimson')
ax.plot(z_a,bestfit(z_a,'Lbreak'),'-',c='seagreen')
ax.plot(z_a,bestfit_global(z_a,3),'-',c='darkorchid')

#limit = np.genfromtxt("../../fitresult/stat_limit.dat",names=True)
#ax.step(limit['z'],limit['Lmin3'],where='mid',linestyle=':',color='chocolate',label=r'$\rm 3$ $\rm data$ $\rm points$ $\rm limit$')
#ax.step(limit['z'],limit['Lmin5'],where='mid',linestyle=':',color='gray',label=r'$\rm 5$ $\rm data$ $\rm points$ $\rm limit$')

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=2,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\log{(L_{\ast}\,[{\rm L}_{\odot}])}$',fontsize=40,labelpad=5)

ax.set_xlim(0,7.)
ax.set_ylim(10.3,13.8)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/Lbreak.png",fmt='png')
####################################################

fig=plt.figure(figsize = (15,10))
ax = fig.add_axes([0.11,0.12,0.79,0.83])

data=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_nofix.dat",names=True)
ax.plot(data["z"],data["redchi"],linestyle='-',marker='o',c='gray',mec='gray',ms=18)

data=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_fix.dat",names=True)
ax.plot(data["z"],data["redchi"],linestyle='-',marker='o',c='royalblue',mec='royalblue',ms=18)

prop = matplotlib.font_manager.FontProperties(size=25.0)
ax.legend(prop=prop,numpoints=1, borderaxespad=0.5,loc=2,ncol=1,frameon=False)
ax.set_xlabel(r'$\rm z$',fontsize=40,labelpad=2.5)
ax.set_ylabel(r'$\chi_{\nu}$',fontsize=40,labelpad=5)

ax.set_xlim(0,7.)
#ax.set_ylim(0.,1.)
ax.tick_params(labelsize=30)
ax.tick_params(axis='x', pad=7.5)
ax.tick_params(axis='y', pad=2.5)
ax.minorticks_on()
plt.savefig("../figs/redchi.pdf",fmt='pdf')
