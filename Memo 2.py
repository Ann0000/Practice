import numpy as np
import matplotlib.pyplot as plt
import emcee
import random
import time
from matplotlib.ticker import MaxNLocator
from scipy import stats

from ed_functions_cosmology2 import *
from ed_functions_memo2 import *


# defines a prior.  just sets acceptable ranges
def lnprior(theta):
    Omega_M, sigma_8, sigma_int, mu_3_int, mu_4_int, sigma_v_kms = theta  # , diag_amp

    if 0.0 < Omega_M < 0.5 and 0.0 < sigma_8 < 3.0 and 0.0 < sigma_int < 0.99 and -0.01 < mu_3_int < 0.01 and 0.0 < mu_4_int < 0.001 and 0.0 < sigma_v_kms < 20.0:
        return 0.0  # -0.5 * ((Omega_M - x_L_max)**2 / (x_sigma_range**2))
    return -np.inf


# lnprob - this just combines prior with likelihood
def lnprob(theta, z_pec, z_cos, mod, mod_error, z_bins, inv_CM_block):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, z_pec, z_cos, mod, mod_error, z_bins, inv_CM_block)


# defines likelihood.  has to be ln likelihood
def lnlike(theta, z_pec, z_cos, mod, mod_error, z_bins, inv_CM_block):
    global use_lensing
    global use_velocity
    global use_velocity_NL
    # global survey


    # ***** set parameters ******
    Omega_M, sigma_8, sigma_int, mu_3_int, mu_4_int, sigma_v_kms = theta  # , diag_amp

    use_intrinsic = 0

    if use_intrinsic == 0:
        mu_3_int = 0.0000000000000001
        mu_4_int = 0.0000000000000001
        sigma_int = 0.14

    sigma_v_kms = sigma_v_kms * 100.0

    Omega_L = 1.0 - Omega_M

    w_o = -1.0
    w_a = 0.0

    # mu_3_int=0.0
    # mu_4_int=0.0
    # sigma_v_kms=0.0

    # ***** write parameters ******
    param_file_name = 'XXX_0720_XXX'

    chain_path = 'XXX_0720_XXX'
    chain_path_file = chain_path + param_file_name
    f_handle = open(chain_path_file, 'a')
    stringOut = str(Omega_M) + ',' + str(sigma_8) + ',' + str(sigma_int) + ',' + str(mu_3_int) + ',' + str(
        mu_4_int) + ',' + str(sigma_v_kms) + '\n'
    f_handle.write(stringOut)
    f_handle.close()

    # ****** Now crunch likelihood **************

    # **************** Subtract best fit ****************

    mod_fid = np.zeros(len(z))

    # Omega_M_fit  =0.25
    # x_sigma_range = 0.2

    delta_m = np.zeros(len(z))  # 'delta_m' is the name for the residuals
    for i in range(0, len(z)):
        mod_fid = distance_modulus(z_pec[i], z_cos[i], Omega_M, Omega_L, w_o, w_a)
        delta_m[i] = mod[i] - mod_fid

    mu_vector_theory_all_bins = []
    mu_vector_data_all_bins = []

    ChSq_total = 0.0

    for j in xrange(0, len(z_bins)):
        inv_CM = inv_CM_block[:, :, j]

        # select within redshift range
        z_min = zedges[j]
        z_max = zedges[j + 1]
        select = (z > z_min) & (z < z_max)

        z_select = z[select]

        delta_m_range = delta_m[select]

        mod_error_range = mod_error[select]

        z_test = weighted_average(z_select, mod_error_range)

        # how many in bin
        N_j = len(delta_m_range) + 1

        # observed moments
        #mu_1_bin_data = weighted_average(delta_m_range, mod_error_range)
        mu_2_bin_data = kde_moment(delta_m_range, mod_error_range, 2)
        #mu_3_bin_data = kde_moment(delta_m_range, mod_error_range, 3)
        #mu_4_bin_data = kde_moment(delta_m_range, mod_error_range, 4)

        # !!!! setup for theory input mode !!!!
        # mu_1_bin_data = distance_modulus( z_test, z_test, Omega_M, Omega_L, w_o, w_a   ) -  distance_modulus( z_test, z_test, 0.25, 0.75, -1.0, 0.0   )
        # mu_2_bin_data = mu_2_total_function( z_test, 0.8, 0.25, 0.14 , 0.001 )
        # mu_3_bin_data = mu_3_total_function( z_test, 0.8, 0.25, 0.0000001 )
        # mu_4_bin_data = mu_4_total_function( z_test, 0.8, 0.25, 0.14 , 0.000000001 , 0.001 )


        mu_vector_data = np.array([mu_2_bin_data])

        # theoretical predictions for modulus in bin
        #mu_1_bin_theory = 0.0   distance_modulus( z_test, z_test, Omega_M, Omega_L, w_o, w_a   )

        # theoretical predictions for moments
        mu_2_bin_theory = mu_2_total_function(z_test, sigma_8, Omega_M, sigma_int, sigma_v_kms)
        #mu_3_bin_theory = mu_3_total_function(z_test, sigma_8, Omega_M, mu_3_int)
        #mu_4_bin_theory = mu_4_total_function(z_test, sigma_8, Omega_M, sigma_int, mu_4_int, sigma_v_kms)

        mu_vector_theory = np.array([mu_2_bin_theory])

        # calculate Chi-squared
        Delta_vector = mu_vector_data - mu_vector_theory

        ChSq = np.dot(Delta_vector, (np.dot(inv_CM, Delta_vector)))
        # print z_test , '  ChSq *** ', ChSq

        ChSq_total = ChSq_total + ChSq

    chsq__file_name = 'ChSqFile_' + 'XXX_0720_XXX'
    chsq_path = 'XXX_0720_XXX'
    path_file = chsq_path + chsq__file_name
    f_handle = open(path_file, 'a')
    stringOut = str(ChSq_total) + '\n'
    f_handle.write(stringOut)
    f_handle.close()

    plotOn = 0
    if plotOn == 1:
        import matplotlib.pyplot as plt

        DoF = (4 * 7) - 6
        plt.title(str(ChSq_total / DoF))
        """plt.subplot(4, 1, 1)
        plt.plot(z_bins, mu_1_all_bins_theory)
        plt.errorbar(z_bins, mu_1_all_bins_data, yerr=uncert_mu_1)
        plt.title(str(round(ChSq_total / DoF, 3)))"""

        plt.subplot(4, 1, 2)
        plt.plot(z_bins, mu_2_bin_theory)
        plt.errorbar(z_bins, mu_2_bin_theory, yerr=uncert_mu_2)
        plt.ylim([0.01, 0.05])

        """plt.subplot(4, 1, 3)
        plt.plot(z_bins, mu_3_all_bins_theory)
        plt.errorbar(z_bins, mu_3_all_bins_data, yerr=uncert_mu_3)
        plt.ylim([-0.001, 0.008])

        plt.subplot(4, 1, 4)
        plt.plot(z_bins, mu_4_all_bins_theory)
        plt.errorbar(z_bins, mu_4_all_bins_data, yerr=uncert_mu_4)
        plt.ylim([-0.001, 0.008])"""

        fileName = 'moments_plot_' + param_file_name + '.pdf'
        plt.savefig(fileName)
        sys.exit()

    return -0.5 * (ChSq_total)  # math.log(np.linalg.det(CM))


# settings, and send off to MeMo module
use_lensing = 0
use_velocity = 1
use_velocity_NL = 0

zero = settings(use_lensing, use_velocity, use_velocity_NL)

survey = 'a'

if survey == 'JLA_4Z':

    # ****** load JLA ****************
    Path = 'Data/'

    suffix = '.txt'
    FileName = Path + survey + suffix
    DataBlock = np.genfromtxt(FileName, skip_header=0, delimiter=',')
    z = DataBlock[:, 3]

    z_pec = DataBlock[:, 2]
    z_cos = DataBlock[:, 3]
    mod = DataBlock[:, 6]
    mod_error = DataBlock[:, 7]

    z = z_cos

    # **************** start values of parameters ****************

    Omega_M = 0.3
    test_Omega_M_min = 0.27
    test_Omega_M_max = 0.23

    sigma_8 = 0.8
    sigma_int = 0.14
    mu_3_int = 0.0001
    mu_4_int = 0.0001
    sigma_v_kms = 2.0
else:
    # **************** load catalogue ****************

    Path = 'sim_JLA_v4_'
    # Path = 'Data/sims_v5/'

    suffix = '.txt'
    FileName = Path + survey + suffix
    DataBlock = np.genfromtxt(FileName, skip_header=0, delimiter=',')
    z = DataBlock[:, 1]

    z_pec = DataBlock[:, 0]
    z_cos = DataBlock[:, 1]
    mod = DataBlock[:, 2]
    mod_error = DataBlock[:, 3]

    z = z_cos
    # delta_m_true = DataBlock[:,4]
    # delta_m_lens = DataBlock[:,5]

    # **************** start values of parameters ****************

    Omega_M = 0.25
    test_Omega_M_min = 0.23
    test_Omega_M_max = 0.27

    sigma_8 = 0.8
    sigma_int = 0.14
    mu_3_int = 0.0001
    mu_4_int = 0.0001
    sigma_v_kms = 2.0

N_bins = 9
zedge_max = 0.9
zedge_min = 0.0
zedges = np.linspace(zedge_min, zedge_max, (N_bins + 1), endpoint=True)
z_bins = (zedges[1:] + zedges[:-1]) / 2

print("z_bins:", z_bins)
print("zedges:",zedges)

# *** bootstrap covariance matrix ***

# other parameters
w_o = -1.0
w_a = 0.0

# **************** Subtract best fit ****************
nSteps = 274
test_Omega_M_all = np.linspace(test_Omega_M_min, test_Omega_M_max, nSteps)
ChSq_all = []
for test_Omega_M in test_Omega_M_all:
    mod_fid = np.zeros(len(z))
    delta_m = np.zeros(len(z))  # 'delta_m' is the name for the residuals
    for l in range(0, len(z)):
        mod_fid[l] = distance_modulus(z_pec[l], z_cos[l], test_Omega_M, (1.0 - test_Omega_M), w_o, w_a)
        delta_m[l] = mod[l] - mod_fid[l]
    ChSq = np.sum((delta_m ** 2.0 / (mod_error ** 2)))
    ChSq_all = np.append(ChSq_all, ChSq)
    print("ChSq:", ChSq)
index_min = np.argmin(ChSq_all)
Omega_M_fit = test_Omega_M_all[index_min]

mod_fid = np.zeros(len(z))

Omega_L = 1.0 - Omega_M_fit

delta_m = np.zeros(len(z))  # 'delta_m' is the name for the residuals
for k in range(0, len(z)):
    mod_fid[k] = distance_modulus(z_pec[k], z_cos[k], Omega_M_fit, Omega_L, w_o, w_a)
    delta_m[k] = mod[k] - mod_fid[k]

# make mu blocks
runs = range(0, 1000)

#mu_1_block = np.zeros((len(runs), N_bins))
mu_2_block = np.zeros((len(runs), N_bins))
#mu_3_block = np.zeros((len(runs), N_bins))
#mu_4_block = np.zeros((len(runs), N_bins))

i = 0

bootstrap_fraction = 0.7

N_total = len(delta_m)

N_boot = int(N_total * bootstrap_fraction)

for run in runs:

    boot_index = random.sample(range(0, N_total), N_boot)

    delta_m_boot = delta_m[boot_index]
    z_boot = z[boot_index]
    mod_boot = mod[boot_index]
    mod_error_boot = mod_error[boot_index]

    # loop to make CM block
    for j in xrange(0, len(z_bins)):
        # select within redshift range
        z_min = zedges[j]
        z_max = zedges[j + 1]
        select = (z_boot > z_min) & (z_boot < z_max)

        z_bin = (z_min + z_max) / 2.0
        z_select = z[select]

        delta_m_range = delta_m_boot[select]

        mod_range = mod[select]
        print("mod[select]:",mod[select])

        # print delta_m_range
        mod_error_range = mod_error[select]

        # z_bin = weighted_average(z_select, mod_error_range)

        # how many in bin
        N_j = len(delta_m_range) + 1

        #mu_1_bin_data = weighted_average(delta_m_range, mod_error_range)
        mu_2_bin_data = kde_moment(delta_m_range, mod_error_range, 2)
        #mu_3_bin_data = kde_moment(delta_m_range, mod_error_range, 3)
        #mu_4_bin_data = kde_moment(delta_m_range, mod_error_range, 4)

        # mu_2_bin_data = kde_moment( delta_m_range, mod_error_range, 2 )
        # mu_3_bin_data = kde_moment( delta_m_range, mod_error_range, 3 )
        # mu_4_bin_data = kde_moment( delta_m_range, mod_error_range, 4 )


        #mu_1_block[i, j] = mu_1_bin_data
        mu_2_block[i, j] = mu_2_bin_data
        #mu_3_block[i, j] = mu_3_bin_data
        #mu_4_block[i, j] = mu_4_bin_data

    i = i + 1

# now chew out covariance matrix

CM_block = np.zeros([1, 1, len(z_bins)])

inv_CM_block = np.zeros([1, 1, len(z_bins)])

block_assemble = 0.0
ind_good = 0

print ''
print ''

print 'CM'

# loop to make CM block
for j in xrange(0, N_bins):

    #mu_1_surveys = mu_1_block[:, j]
    mu_2_surveys = mu_2_block[:, j]
    #mu_3_surveys = mu_3_block[:, j]
    #mu_4_surveys = mu_4_block[:, j]

    # plt.scatter(mu_2_surveys,mu_3_surveys,alpha=0.3)
    # plt.show()

    mu_z_block = mu_2_surveys.T

    CM = np.cov(mu_z_block)

    # CM = CM / float(bootstrap_fraction)

    inv_CM = 1./CM

    if block_assemble > 0.0:
        CM_block = np.dstack([CM_block, CM])  # [:,:,ind_good] = inv_CM
        inv_CM_block = np.dstack([inv_CM_block, inv_CM])  # [:,:,ind_good] = inv_CM

        ind_good = ind_good + 1

    else:
        CM_block = CM
        inv_CM_block = inv_CM

        block_assemble = 1.0
        ind_good = ind_good + 1

    print ''
    print ''
    print ("j",j)
    print ("CM",CM)

# **************** MCMC ****************
print 'hi'
#


startValues = [Omega_M, sigma_8, sigma_int, mu_3_int, mu_4_int, sigma_v_kms]

# how many parameters to fit
ndim = len(startValues)

# how many walkers
nwalkers = 50
nSteps = 200

# plt.scatter(z,mod,color = 'c')
# plt.scatter(z,mod_fid,color = 'y')
# plt.plot(z_bins_good, mu_1_all_bins_theory,'b')
# plt.plot(z_bins_good, mu_1_all_bins_data,'r')
# plt.show()
# sys.exit()

pos = [startValues + 1e-3 * np.random.randn(ndim) for i in range(nwalkers)]
print("pos:",pos)

# setup the sampler
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(z_pec, z_cos, mod, mod_error, z_bins, inv_CM_block),
                                threads=1)
# run the sampler
# how many steps (will have nSteps*nwalkers of samples)
#sampler.run_mcmc(pos, nSteps)


#start time.time
tnow=time.time()
for i, result in enumerate(sampler.sample(pos,iterations=nSteps)):
    print(time.time()-tnow)
    print i
    tnow=time.time()


#print samples
plt.clf()
fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 9))
axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
axes[0].yaxis.set_major_locator(MaxNLocator(5))
axes[0].axhline(Omega_M, color="#888888", lw=2)
axes[0].set_ylabel("$Omega_M$")

axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
axes[1].yaxis.set_major_locator(MaxNLocator(5))
axes[1].axhline(sigma_8, color="#888888", lw=2)
axes[1].set_ylabel("$sigma_8$")

axes[2].plot(np.exp(sampler.chain[:, :, 2]).T, color="k", alpha=0.4)
axes[2].yaxis.set_major_locator(MaxNLocator(5))
axes[2].axhline(sigma_int, color="#888888", lw=2)
axes[2].set_ylabel("$sigma_int$")
axes[2].set_xlabel("step number")

fig.tight_layout(h_pad=0.0)
fig.savefig("line-time.png")






samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

print("samples.shape:",samples.shape)
print(samples[:,0])
print(samples[:,1])
print(samples[:,2])
print(samples[:,3])

import corner

fig = corner.corner(samples, labels=["Omega_M", "sigma_8","sigma_int", "mu_3_int", "mu_4_int", "sigma_v_kms"],truths=[Omega_M, sigma_8,sigma_int, mu_3_int, mu_4_int, sigma_v_kms])
fig.savefig("triangle.png")
plt.show()