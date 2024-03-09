import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy import constants as const

def dL(chi2, OmegaK, x):
    A_sin = np.sin(np.sqrt(np.abs(OmegaK)) * const.Hubble * np.sqrt(chi2) / const.c)
    A_sinh = np.sinh(np.sqrt(np.abs(OmegaK)) * const.Hubble * np.sqrt(chi2) / const.c)
    B = np.sqrt(np.abs(OmegaK)) * const.H0 * np.sqrt(chi2) / const.c


    if OmegaK == 0:
        r = np.sqrt(chi2)
    elif get_OmegaK(x) < 0:
        r = np.sqrt(chi2) * A_sin / B
    elif get_OmegaK(x) > 0:
        r = np.sqrt(chi2) * A_sinh / B

    dL = r / (np.exp(x))


    return dL

Gyr = 1/(60*60*24*365*1e9) # from s to Gyr
Mpc = 3.24*10**(-23) # from m to Mpc
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


selected_data = converged_data[converged_data[:, 0] < best_fit_params[0] + 3.53] #Data selected within 1sigma of the best fit

chi2 = selected_data[:, 0]
OmegaM_selected = selected_data[:, 2]
OmegaK_selected = selected_data[:, 3]
OmegaLambda_selected = 1 - (OmegaK_selected + OmegaM_selected)


""" dL plot with supernova best fit, fiducial cosmology and betouli observations """

betoule = np.loadtxt("Betoule_supernova.txt")


z_cosmo = np.exp(cosmo_x)-1
z_obs = betoule[:, 0]


dL_obs = betoule[:, 1]
error_obs = betoule[:, 2]

dL_fit = dL(chi2, OmegaK_selected, cosmo_x)

plt.plot(z_cosmo, cosmo_dL/z_cosmo, label="Fiducial cosmology")
plt.errorbar(z_obs, dL_obs, yerr=error_obs, fmt='o', color='blue', ecolor='red', capsize=0.5, label="Observed data")
plt.plot(z_cosmo, dL_fit, label="Best fit from MCMC")
plt.show()


"""
plt.scatter(OmegaM_selected, OmegaLambda_selected)
plt.xlabel('OmegaM')
plt.ylabel('OmegaLambda')
plt.title('1$\sigma$ Confidence Region')
plt.show()


std_OmegaM = np.std(OmegaM_selected)
std_OmegaLambda = np.std(OmegaLambda_selected)
print("Standard deviation of OmegaM:", std_OmegaM)
print("Standard deviation of OmegaLambda:", std_OmegaLambda)

plt.hist(OmegaM_selected, bins=20, alpha=0.5, label='OmegaM')
plt.hist(OmegaLambda_selected, bins=20, alpha=0.5, label='OmegaLambda')
plt.xlabel('Parameter Value')
plt.ylabel('Frequency')
plt.legend()
plt.title('Histogram of Parameters')
plt.show()


# Compute the mean and standard deviation of the chain values
mean_OmegaM = np.mean(OmegaM_selected)
mean_OmegaLambda = np.mean(OmegaLambda_selected)

# Compute the standard deviation
std_OmegaM = np.std(OmegaM_selected)
std_OmegaLambda = np.std(OmegaLambda_selected)

# Generate values for the Gaussian distribution
x = np.linspace(min(OmegaM_selected), max(OmegaM_selected), 100)
gaussian_OmegaM = 1 / (std_OmegaM * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - mean_OmegaM) / std_OmegaM)**2)

y = np.linspace(min(OmegaLambda_selected), max(OmegaLambda_selected), 100)
gaussian_OmegaLambda = 1 / (std_OmegaLambda * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((y - mean_OmegaLambda) / std_OmegaLambda)**2)

# Plot the histogram and Gaussian distribution
plt.hist(OmegaM_selected, bins=20, alpha=0.5, label='OmegaM')
plt.plot(x, gaussian_OmegaM, color='blue', linestyle='--', label='Gaussian fit (OmegaM)')

plt.hist(OmegaLambda_selected, bins=20, alpha=0.5, label='OmegaLambda')
plt.plot(y, gaussian_OmegaLambda, color='orange', linestyle='--', label='Gaussian fit (OmegaLambda)')

plt.xlabel('Parameter Value')
plt.ylabel('Frequency')
plt.legend()
plt.title('Histogram of Parameters with Gaussian Fit')
plt.show()

""" 


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