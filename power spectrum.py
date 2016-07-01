import cosmolopy
import numpy as np
import matplotlib.pyplot as pl
from ed_functions_memo2 import sigma_v_function
from ed_functions_memo2 import sigma_v_nonlinear
z=np.arange(0.01,1.,0.01)
sigma_v1=[]
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
pl.show()
# The first try to draw plot sigma_v--z.
# The next step is to


