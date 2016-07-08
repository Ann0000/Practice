import cosmolopy
import numpy as np
import matplotlib.pyplot as pl
from ed_functions_memo2 import sigma_v_function
from ed_functions_cosmology2 import Ez


np_sigma=np.frompyfunc(sigma_v_function,3,1)

def  sigma_approx(z,sigma_8,omega_m):
    c_z1=0.12201-0.43413*omega_m+1.55192*omega_m**2-1.64078*omega_m**3
    c_z2=-4.4354+2.32415*np.exp(-11.335*omega_m)
    c_z3=2.68
    c_z4=1.1608-0.8429*omega_m

    sigma_v=0.006*sigma_8*np.exp(c_z1*z**-1+c_z2*z+c_z3*z**1.3+c_z4)

    return sigma_v

"""z=0.2
omega_m=20
sigma_8=20
i=0
j=0
sigma_v1=-np_sigma(z,0.2,0.1)
for i in range(omega_m):
    omega_m_i=0.2+0.01*i
    for j in range(sigma_8):
        sigma_8_j=j*0.15+0.1
        sigma_v1_i_j=-np_sigma(z, sigma_8_j,omega_m_i )
        sigma_v1.append(sigma_v1_i_j)
        print(sigma_v1)

#sigma_approx1 = sigma_approx(z, sigma_8, omega_m)
#sigma_v1=np.array([(x*y) for x in x_axis for y in y_axis ])
print(sigma_v1.size)
sigma_v1=np.reshape(sigma_v1,(20,20))"""

z=0.2
omega_m=np.arange(0.2,0.4,0.01)
sigma_8=np.arange(0.,3.,0.1)
sigma_8_j=0.1
omega_m_i=0.2
n_omega=omega_m.shape
print(n_omega)
n_sigma=sigma_8.shape
print(np_sigma(z,sigma_8_j,omega_m_i))
sigma_v=np.array([(sigma_8_j*omega_m_i)   for sigma_8_j in sigma_8 for omega_m_i in omega_m])

sigma_v=sigma_v.reshape(30,20)




pl.ylabel('sigma_8')
pl.xlabel('omega_m')
pl.contour(omega_m,sigma_8,sigma_v)
pl.title('z=0.2')
pl.show()






