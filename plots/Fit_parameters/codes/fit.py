from scipy.optimize import curve_fit as cfit 
import numpy as np

T0 = np.polynomial.chebyshev.Chebyshev((1,0,0,0))
T1 = np.polynomial.chebyshev.Chebyshev((0,1,0,0))
T2 = np.polynomial.chebyshev.Chebyshev((0,0,1,0))
T3 = np.polynomial.chebyshev.Chebyshev((0,0,0,1))

data=np.genfromtxt("../../fitresult/fit_at_z.dat",names=True)

def func1(z,*p):
	xsi=1.+z
	return p[0]*T0(xsi)+p[1]*T1(xsi)#+p[2]*T2(xsi)+p[3]*T3(xsi)

def func2(z,*p):
	xsi=1.+z
	return p[0]*T0(xsi)+p[1]*T1(xsi)+p[2]*T2(xsi)#+p[3]*T3(xsi)

def func3(z,*p):
	xsi=1.+z
	return p[0]*T0(xsi)+p[1]*T1(xsi)+p[2]*T2(xsi)+p[3]*T3(xsi)

def doublepower(z,*p):
	xsi=1.+z
	zref=p[1]
	return 2*p[0]/(np.power(xsi/(1+zref),p[2]) + np.power(xsi/(1+zref),p[3]))

id= (data["z"]!=4) & (data["z"]!=4.5)
fit1=cfit(func3,data["z"][id],data["gamma1"][id],sigma=data["err1"][id],p0=np.array([0,0,0,0]))[0]
id= (data["z"]!=5) & (data["z"]!=6)
fit2=cfit(doublepower,data["z"][id],data["gamma2"][id],sigma=data["err2"][id],p0=np.array([0,2,0,0]))[0]
id= (data["z"]<=3)
fit3=cfit(func1,data["z"][id],data["phi_s"][id],sigma=data["err3"][id],p0=np.array([0,0,0,0]))[0]
id= (data["z"]<5)
fit4=cfit(doublepower,data["z"][id],data["L_s"][id],sigma=data["err4"][id],p0=np.array([0,2,0,0]))[0]

print "intepreted logphi_s:"
print func1(np.array([3.5,4,4.5,5,5.5,6,6.5]),*fit3)

print 'gamma1',fit1
print 'gamma2',fit2
print 'logphi_s',fit3
print 'Lbreak',fit4

np.savetxt("zevolution_fit.dat",np.c_[fit1,fit2,fit3,fit4],header='gamma1 gamma2 phi_s Lbreak')