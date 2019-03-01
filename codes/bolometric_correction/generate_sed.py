from data import *
import numpy as np 
import astropy.constants as con

data1=np.genfromtxt(datapath+"K13_SED.dat",names=["lognu","logall"],)
freq1=data1['lognu']
lamb1=con.c.value/(10**data1['lognu'])*1e10
sed1 =data1['logall']

#IR construction
IR_end=lamb1[0]
data2=np.genfromtxt(datapath+"R06_SED.dat",names=["lognu","logall","sigall","blue"],)
lamb2=con.c.value/(10**data2['lognu'])*1e10
sed2 =data2['blue']

lamb2_ir= lamb2[lamb2>IR_end]
sed2_ir = sed2[lamb2>IR_end]

scale = ((sed1[1]-sed1[0])/(lamb1[1]-lamb1[0]) * (lamb2_ir[-1]-lamb1[0]) + sed1[0])/sed2_ir[-1]

lamb1= np.append( lamb2_ir, lamb1)
sed1 = np.append( sed2_ir*scale, sed1)

#UV construction
UV_end  =912.
UVslope =1.70

lamb1,sed1 = lamb1[lamb1>912],sed1[lamb1>912]

lambUV= np.linspace(912.,500.,10)
sedUV = sed1[-1] + (UVslope-1) * (np.log10(lambUV) - np.log10(lamb1[-1]))

lamb1= np.append( lamb1, lambUV)
sed1 = np.append( sed1 , sedUV)
freq1= np.log10(con.c.value/(lamb1*1e-10))

#Xray construction and connection
f2500 = np.power(10.,sed1[np.abs(lamb1-2500)<=100])/np.power(10.,freq1[np.abs(lamb1-2500)<=100])
alphaox = -0.137*np.log10(f2500)+2.638
f2kev=10**(alphaox/0.384)*f2500
data3=np.genfromtxt(datapath+"XRAY_SED.dat",names=["lognu","logall"],)
lamb3=con.c.value/(10**data3['lognu'])*1e10
ratio= data3['logall'][np.abs(10**data3['lognu']*con.h.value/con.e.value/1000-2)<=0.05]-(-1.47406)
L_HX = np.log10( f2kev*(2.*1000.*con.e.value/con.h.value) / 10**ratio )
sed3 = data3['logall']-(-1.47406)+L_HX

lamb1= np.append( lamb1, lamb3[lamb3<50])
sed1 = np.append( sed1 , sed3[lamb3<50])

np.savetxt(datapath+"MySED.dat",np.c_[lamb1,sed1])
