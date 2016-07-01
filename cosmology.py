import cosmolopy.distance as cd
"""cosmo={'omega_M_0':0.3,'omega_lambda_0':0.7,'omega_k_0':0.0,'h':0.72}
d_co=cd.comoving_distance(6.,**cosmo)
print("Comoving distance to z=6 is %.lf Mpc"%(d_co))

from cosmolopy import*
d_a=cd.angular_diameter_distance(6,**fidcosmo)
print("Angular-diameter distance to z=6 is %.lf Mpc"%(d_a))
d_light=cd.light_travel_distance(6,**fidcosmo)
print"Light-travel distance to z=6 is %.lf Mpc"%(d_light)"""

import cosmolopy.perturbation as cp
import numpy as np
import math
log_k=np.arange(np.log10(2*math.pi/10000) ,np.log10(2*math.pi/40),0.05) #k Mpc^-1
print(log_k.shape)
k=10**log_k
print("k=",k)
cosmology={'omega_M_0' : 0.308,'omega_b_0' : 0.022,'omega_n_0' : 0.0,'N_nu' : 3,'omega_lambda_0' : 0.692,
'h' : 0.72,'n' : 0.95,'sigma_8' : 0.8}
spectrum=cp.power_spectrum(10**log_k,1.,**cosmology)
print("spectrum" ,spectrum)
import matplotlib.pyplot as mpl
mpl.plot(log_k,spectrum,color="r",lw=2, alpha=0.8)
mpl.ylabel('Power spectrum')
mpl.xlabel('wave number log10_k (Mpc^-1)')
mpl.show()

cosmology=cd.set_omega_k_0(cosmology)
print(cosmology['omega_k_0'])

from ed_functions_memo2 import sigma_v_function
from ed_functions_memo2 import sigma_v_nonlinear
