import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches




""" SUPERNOVA FITTING"""
data = np.loadtxt("results_supernovafitting.txt")

converged_data = data[200:]
print(f"The shape of the converged data is {np.shape(converged_data)}\n")
best_chi = np.argmin(converged_data[:, 0])

best_fit_params = converged_data[best_chi, :]
print(f"Best params from min chi^2 {best_fit_params}\n")


selected_data = converged_data[converged_data[:, 0] < best_fit_params[0] + 3.53] #Data selected within 1sigma of the best fit


OmegaM_selected = selected_data[:, 2]
OmegaK_selected = selected_data[:, 3]
OmegaLambda_selected = 1 - (OmegaK_selected + OmegaM_selected)

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