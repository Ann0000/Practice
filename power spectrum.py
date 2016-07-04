import cosmolopy
import numpy as np
import matplotlib.pyplot as pl
import math
from ed_functions_memo2 import sigma_v_function
from ed_functions_memo2 import sigma_v_nonlinear
"""z=np.arange(0.01,1.,0.01)
np_sigma=np.frompyfunc(sigma_v_function,3,1)
sigma_v1=np_sigma(z,0.8,0.3)
sigma_v2=np_sigma(z,0.6,0.3)
pl.plot(z,sigma_v1,color='r')
pl.plot(z,sigma_v2,color='b')
pl.show()"""
"""sigma_v1=[]
sigma_v2=[]
for z_i in range(1,100,1):
    sigma_v1_i=(-sigma_v_function(z_i/100.,0.8,0.3))
    sigma_v1.append(sigma_v1_i)
    sigma_v2_i=(-sigma_v_function(z_i/100.,0.6,0.3))
    sigma_v2.append(sigma_v2_i)
pl.plot(z,sigma_v1,color='r')
pl.plot(z,sigma_v2,color='b')
pl.ylabel('sigma_v1 & sigma_v2')
pl.xlabel('redshift z')
pl.show()"""

# The first try to draw plot sigma_v--z.
# 21:52 write another function
# The next step is to vary sigma_8 and get the function [seems a quite hard job!

def f(x):
    return x-math.exp(x)+50*x**3+math.log1p(x)


z=np.arange(0.1,1.,0.01)
omega_m=np.arange(0.2,0.4,0.01)
sigma_8=np.arange(0.5,1.5,0.05)
np_sigma=np.frompyfunc(sigma_v_function,3,1)
sigma_v1=-np_sigma(z,0.6,0.25)
y1=0.006*np.exp(0.07*z**-1+z**1.5-4*z+0.6)
sigma_v2=-np_sigma(z,0.6,0.31)
y2=0.005*np.exp(0.09*z**-1+z**1.5-4*z+0.6)
sigma_v3=-np_sigma(z,0.6,0.27)
y3=0.0055*np.exp(0.08*z**-1+z**1.5-4*z+0.6)
#sigma_v4=-np_sigma(0.3,0.8,omega_m)
pl.plot(z,sigma_v1,color='r')
pl.plot(z,y1,color='b')
pl.plot(z,0.01*sigma_v1/y1,color='g')
pl.plot(z,sigma_v2,'--',color='r')
pl.plot(z,y2,'--',color='b')
pl.plot(z,0.01*sigma_v2/y2,'--',color='g')
pl.plot(z,sigma_v3,'.',color='r')
pl.plot(z,y3,'.',color='b')
pl.plot(z,0.01*sigma_v3/y3,'.',color='g')
#pl.plot(omega_m,sigma_v4,color='k')
pl.ylabel('sigma_v1 & sigma_v2')
pl.xlabel('redshift z')

pl.show()
# fitting sigma_8 and z