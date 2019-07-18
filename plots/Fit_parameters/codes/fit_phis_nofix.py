from scipy.optimize import curve_fit as cfit 
import numpy as np

T0 = np.polynomial.chebyshev.Chebyshev((1,0,0,0))
T1 = np.polynomial.chebyshev.Chebyshev((0,1,0,0))
T2 = np.polynomial.chebyshev.Chebyshev((0,0,1,0))
T3 = np.polynomial.chebyshev.Chebyshev((0,0,0,1))

data=np.genfromtxt("../../../codes/lf_fit/output/fit_at_z_nofix.dat",names=True)

def func1(z,*p):
	xsi=1.+z
	return p[0]*T0(xsi)+p[1]*T1(xsi)#+p[2]*T2(xsi)+p[3]*T3(xsi)

id= (data["z"]>=0.4) & (data['z']<=2.8)
fit3, cov3=cfit(func1,data["z"][id],data["phi_s"][id],sigma=data["err3"][id],p0=np.array([-3.5,-0.35,0.,0.]))

print fit3
