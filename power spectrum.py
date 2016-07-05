import cosmolopy
import numpy as np
import time
tStart=time.time()
import matplotlib.pyplot as pl
from ed_functions_memo2 import sigma_v_function
from ed_functions_memo2 import sigma_v_nonlinear


z=np.arange(0.1,1.,0.01)
#omega_m=np.arange(0.2,0.4,0.01)
#sigma_8=np.arange(0.5,1.5,0.05)
#np_sigma_nl=np.frompyfunc(sigma_v_nonlinear,3,1)
np_sigma=np.frompyfunc(sigma_v_function,3,1)

def  sigma_approx(z,sigma_8,omega_m):
    c_z1=0.06582+4.06445*1E-5*np.exp(17.453*omega_m)
    c_z2=0.36981+33.76946*np.exp(-13.83038*omega_m)
    c_z3=-9.42233+27.8714*omega_m-31.1962*omega_m**2
    c_z4=1.22587+1.15218*omega_m-6.72929*omega_m**2
    sigma_v=0.006*sigma_8*np.exp(c_z1*z**-1+c_z2*z**1.5+c_z3*z+c_z4)
    return sigma_v

sigma_8=2.8
omega_m=0.36
sigma_v1 = -np_sigma(z, sigma_8,omega_m )
tEnd=time.time()
print(tEnd-tStart)
sigma_approx1 = sigma_approx(z, sigma_8, omega_m)

pl.ylabel('sigma_v1 & sigma_approx1')
pl.xlabel('redshift z')
pl.plot(z, sigma_v1, color='r')
pl.plot(z, sigma_approx1, "--")
pl.title('sigma_8:1.,omega_m:0.36')
pl.show()

"""omega_m100=20
sigma_8=1.
for i in range(omega_m100) :
    omega_m = (i+20)* 1. / 100
    print(omega_m)
    sigma_v1 = -np_sigma(z, sigma_8,omega_m )
    sigma_approx1 = sigma_approx(z, sigma_8, omega_m)
    pl.ylabel('sigma_v1 & sigma_approx1')
    pl.xlabel('redshift z')
    pl.plot(z, sigma_v1, color='r')
    pl.plot(z, sigma_approx1, "--")
    pl.title("sigma_8:1.,(omega_m-0.2)*100:{0:20}".format(i))
    pl.show()
"""




"""y1=0.006*np.exp(0.07*z**-1+z**1.5-4*z+0.6)
sigma_v2=-np_sigma(z,0.6,0.31)
y2=0.005*np.exp(0.09*z**-1+z**1.5-4*z+0.6)
sigma_v3=-np_sigma(z,0.6,0.27)
y3=0.0055*np.exp(0.08*z**-1+z**1.5-4*z+0.6)"""
#sigma_v4=-np_sigma(0.3,0.8,omega_m)
"""pl.plot(z,y1,color='b')
pl.plot(z,0.01*sigma_v1/y1,color='g')
pl.plot(z,sigma_v2,'--',color='r')
pl.plot(z,y2,'--',color='b')
pl.plot(z,0.01*sigma_v2/y2,'--',color='g')
pl.plot(z,sigma_v3,'.',color='r')
pl.plot(z,y3,'.',color='b')
pl.plot(z,0.01*sigma_v3/y3,'.',color='g')
#pl.plot(omega_m,sigma_v4,color='k')"""






# The first try to draw plot sigma_v--z.
# The next step is to vary sigma_8 and get the function [seems a quite hard job!


