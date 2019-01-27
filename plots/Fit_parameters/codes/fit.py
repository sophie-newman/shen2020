from scipy.optimize import curve_fit as cfit 
import numpy as np

T0 = np.polynomial.chebyshev.Chebyshev((1,0,0,0))
T1 = np.polynomial.chebyshev.Chebyshev((0,1,0,0))
T2 = np.polynomial.chebyshev.Chebyshev((0,0,1,0))
T3 = np.polynomial.chebyshev.Chebyshev((0,0,0,1))

data=np.genfromtxt("../../fitresult/fit_at_z.dat",names=True)

def func1(z,*p):
	xsi=1.+z
	return p[0]*T0(xsi)+p[1]*T1(xsi)+p[2]*T2(xsi)#+p[3]*T3(xsi)

def func2(z,*p):
	xsi=1.+z
	return p[0]*T0(xsi)+p[1]*T1(xsi)+p[2]*T2(xsi)+p[3]*T3(xsi)


id= (data["z"]!=4) & (data["z"]!=4.5) & (data["z"]!=5)

fit1=cfit(func2,data["z"][id],data["gamma1"][id],p0=np.array([0,0,0,0]))[0]
fit2=cfit(func2,data["z"][id],data["gamma2"][id],p0=np.array([1,1,1,1]))[0]
fit3=cfit(func1,data["z"][id],data["phi_s"][id],p0=np.array([1,1,1,0]))[0]
fit4=cfit(func1,data["z"][id],data["L_s"][id],p0=np.array([1,1,1,0]))[0]

print 'gamma1',fit1
print 'gamma2',fit2
print 'logphi_s',fit3
print 'Lbreak',fit4

np.savetxt("zevolution_fit.dat",np.c_[fit1,fit2,fit3,fit4],header='gamma1 gamma2 phi_s Lbreak')