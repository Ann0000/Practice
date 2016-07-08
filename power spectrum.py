import cosmolopy
import numpy as np
import matplotlib.pyplot as pl
from ed_functions_memo2 import sigma_v_function


np_sigma=np.frompyfunc(sigma_v_function,3,1)

def  sigma_approx(z,sigma_8,omega_m):
    c_z1=0.12201-0.43413*omega_m+1.55192*omega_m**2-1.64078*omega_m**3
    c_z2=-4.4354+2.32415*np.exp(-11.335*omega_m)
    c_z3=2.68
    c_z4=1.1608-0.8429*omega_m

    sigma_v=0.006*sigma_8*np.exp(c_z1*z**-1+c_z2*z+c_z3*z**1.3+c_z4)

    return sigma_v



z=0.35
omega_m=np.arange(0.2,0.4,0.01)
sigma_8=np.arange(0.1,3.,0.1)
sigma_8_j=0.1
omega_m_i=0.2

sigma_8,omega_m=np.meshgrid(sigma_8,omega_m)
sigma_v1=np_sigma(z,sigma_8,omega_m)
sigma_approx1 = sigma_approx(z, sigma_8, omega_m)

pl.ylabel('sigma_8')
pl.xlabel('omega_m')
pl.contour(omega_m,sigma_8,sigma_v1)
pl.contour(omega_m,sigma_8,sigma_approx1)
pl.title('z=0.35')
pl.show()






