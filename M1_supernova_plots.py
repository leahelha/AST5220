import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy import constants as const
from scipy.stats import norm

Gyr = 1/(60*60*24*365*1e9) # from s to Gyr
Mpc = 3.24*10**(-23) # from m to Mpc
Gpc = 3.24*10**(-25)

cosmo = np.loadtxt("cosmology.txt")
print(f"Shape of cosmo = {np.shape(cosmo)}")

""" All the parameters i will be using, grabbed from the txt file """
cosmo_x = cosmo[:,0]
cosmo_eta_of_x = cosmo[:,1]
cosmo_t_of_x = cosmo[:, 11]

cosmo_Hp = cosmo[:,2]
cosmo_dHpdx = cosmo[:,3]
cosmo_ddHpddx = cosmo[:, 10]

cosmo_OmegaB = cosmo[:, 4]
cosmo_OmegaCDM = cosmo[:,5]
cosmo_OmegaLambda = cosmo[:,6]
cosmo_OmegaR = cosmo[:,7]
cosmo_OmegaNu = cosmo[:,8]
cosmo_OmegaK = cosmo[:,9]

cosmo_dL = cosmo[:, 12]


""" SUPERNOVA FITTING"""
data = np.loadtxt("results_supernovafitting.txt")

converged_data = data[200:]
print(f"The shape of the converged data is {np.shape(converged_data)}\n")
best_chi = np.argmin(converged_data[:, 0])

best_fit_params = converged_data[best_chi, :]
print(f"Best params from min chi^2 {best_fit_params}\n")


selected_data = converged_data[converged_data[:, 0] < (best_fit_params[0] + 3.53)] #Data selected within 1sigma of the best fit
sig2_data = converged_data[converged_data[:, 0] < (best_fit_params[0] + 0.267)] #Data selected within 2sigma of the best fit
sig3_data = converged_data[converged_data[:, 0] < (best_fit_params[0] + 0.133)] #Data selected within 3sigma of the best fit

chi2 = selected_data[:, 0]
h_selected = selected_data[:, 1]
OmegaM_selected = selected_data[:, 2]
OmegaK_selected = selected_data[:, 3]

OmegaLambda_selected = 1 - (OmegaK_selected + OmegaM_selected)


""" dL plot with supernova best fit, fiducial cosmology and betouli observations """

betoule = np.loadtxt("Betoule_supernova.txt")


z_cosmo = np.exp(-cosmo_x)-1
z_obs = betoule[:, 0]


dL_obs = betoule[:, 1]
error_obs = betoule[:, 2]

""" Plot of luminosity distance for fiducial cosmology, observed sn data and our best fit results """
# *** THIS IS NOT RIGHT
plt.figure()
plt.plot(z_cosmo, cosmo_dL*Gpc/z_cosmo, label="Fiducial cosmology")
plt.errorbar(z_obs, dL_obs, xerr=error_obs, fmt='o', color='blue', ecolor='red', capsize=0.5, label="Observed data")
plt.plot(chi2*Gpc, label="Best fit from MCMC")
plt.title("This is not correct")
plt.legend()
plt.savefig("Figs/sn_dL_plots.pdf")


""" Confidence region 1sig, 2sig and 3sig """
plt.figure()
plt.scatter( OmegaLambda_selected, OmegaM_selected, label = "$1\sigma$")
plt.scatter(1 -(sig2_data[:,2]+sig2_data[:,3]), sig2_data[:,2], label = "$2\sigma$")
plt.scatter(1 -(sig3_data[:,2]+sig3_data[:,3]), sig3_data[:,2], label = "$3\sigma$")
# *** wanna add line for flat universe
plt.legend()
plt.xlabel('OmegaM')
plt.ylabel('OmegaLambda')
plt.title('1$\sigma$ Confidence Region')
plt.savefig("Figs/sn_Confidence_region.pdf")


std_OmegaM = np.std(OmegaM_selected)
std_OmegaLambda = np.std(OmegaLambda_selected)
std_OmegaK= np.std(OmegaK_selected)
std_h= np.std(h_selected)
print("Standard deviation of OmegaM:", std_OmegaM)
print("Standard deviation of OmegaLambda:", std_OmegaLambda)


""" Histogram of accepted H """
# *** FIX UNITS AND AXIS LABEL
plt.figure()
plt.hist(h_selected, bins=20, alpha=0.5, label='h')  #*** supposed to be H
plt.xlabel('Parameter Value')
plt.ylabel('Frequency')
plt.legend()
plt.title('Histogram of accepted h')
plt.savefig("Figs/sn_Histogram_of_H_parameters.pdf")


# Compute the mean and standard deviation of the chain values
mean_OmegaM = np.mean(OmegaM_selected)
mean_OmegaLambda = np.mean(OmegaLambda_selected)
mean_OmegaK= np.mean(OmegaK_selected)

mean_Omega_sum = mean_OmegaK + mean_OmegaLambda + mean_OmegaM

Omega_sum = OmegaK_selected + OmegaLambda_selected + OmegaM_selected
#mean_h= np.mean(h_selected)

std_Omega_sum = std_OmegaK + std_OmegaLambda + std_OmegaM # *** Check if the math is on my side here
pdf = norm.pdf(cosmo_x, loc=mean_Omega_sum, scale=std_Omega_sum)


""" Plot the histogram and Gaussian distribution """
# *** THIS IS ALSO WRONG
# *** FIX UNITS AND AXIS LABEL
plt.figure()
plt.hist(OmegaM_selected, bins=20, alpha=0.5, label='OmegaM')
plt.hist(OmegaLambda_selected, bins=20, alpha=0.5, label='OmegaLambda')
plt.hist(OmegaK_selected, bins=20, alpha=0.5, label='OmegaK')
plt.hist(OmegaM_selected+OmegaLambda_selected+OmegaK_selected, bins=20, alpha=0.5, label='$\Omega_i$')
plt.plot(cosmo_x, pdf)

plt.xlabel('Parameter Value')
plt.ylabel('Frequency')
plt.legend()
plt.title('Histogram of Parameters with Gaussian Fit')
plt.savefig("Figs/sn_Histogram_of_omega_Gaussian.pdf")




#           chi2             h      OmegaM    OmegaK    Acceptrate
# Minimum chi^2 found    29.2799   0.701711   0.255027   0.0789514


# Minimum chi^2 found 29.2799 0.701711 0.255027 0.0789514 


# Minimum chi^2 found 29.2799 0.701711 0.255027 0.0789514 

"""
// How to analyze the resulting chains:
// * Load the chains and skip the first few hundred samples (the burnin of the chains). E.g. loadtxt(file,skiprows=200) in python
// * Find the minimum chi2 and the corresponding best-fit parameters (you can use np.argmin to get index of the minvalue in python)
// * Select all samples of OmegaM and OmegaLambda (computed from OmegaM and OmegaK) that satisfy chi2 < chi2_min + 3.53 
//   (e.g. OmegaM[chi2 < chi2min + 3.53] in python)
// * Scatterplotting these gives you the 1sigma (68.4%) confidence region
// * Find the standard deviation of the samples to get the 1sigma confidence region of the parameters (assuming the posterior is a gaussian)
// * Make and plot a histogram of the samples for the different parameters (OmegaM, OmegaK, OmegaLambda, H0)
// * You can also compute the mean and standard deviation of the chain values and use this to overplot a gaussian with the same mean and variance for comparison.

"""