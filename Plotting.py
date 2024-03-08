import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


""" COSMOLOGY PLOTS"""
"""

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp << ddHpddx_of_x(x)         << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

"""
cosmo = np.loadtxt("cosmology.txt")
print(f"Shape of cosmo = {np.shape(cosmo)}")
cosmo_x = cosmo[:,0]
cosmo_eta_of_x = cosmo[:,1]
cosmo_Hp = cosmo[:,2]
cosmo_dHpdx = cosmo[:,3]
cosmo_ddHpddx = cosmo[:, 10]

cosmo_OmegaB = cosmo[:, 4]
cosmo_OmegaCDM = cosmo[:,5]
cosmo_OmegaLambda = cosmo[:,6]
cosmo_OmegaR = cosmo[:,7]
cosmo_OmegaNu = cosmo[:,8]
cosmo_OmegaK = cosmo[:,9]


""" 3 Plots that demonstrate that your code works properly: """
### Regime lines

border_idx1 = np.where((np.abs((cosmo_OmegaR+cosmo_OmegaNu)-(cosmo_OmegaB+cosmo_OmegaCDM)))<0.085)[0]


border_idx2 = np.where((np.abs((cosmo_OmegaB+cosmo_OmegaCDM)-cosmo_OmegaLambda))<0.085)[0]


idx1 = border_idx1[0]
idx2 = border_idx2[-1]

print(idx1)
print(idx2)

# plt.patches(cosmo_x[:idx1], color='red', alpha=0.5)
# plt.pat(cosmo_x[idx1:idx2],color='blue', alpha=0.5)
# plt.fill_between(cosmo_x[idx2:-1], color='purple', alpha=0.5)

# Create patches for filling between lines
region1 = patches.Rectangle((cosmo_x[0], 0), cosmo_x[idx1]-cosmo_x[0], 1, color='red', alpha=0.2)
region2 = patches.Rectangle((cosmo_x[idx1], 0), cosmo_x[idx2] - cosmo_x[idx1], 1, color='blue', alpha=0.2)
region3 = patches.Rectangle((cosmo_x[idx2], 0), cosmo_x[-1] - cosmo_x[idx2], 1, color='purple', alpha=0.2)

# Add the patches to the plot
plt.gca().add_patch(region1)
plt.gca().add_patch(region2)
plt.gca().add_patch(region3)

plt.plot(cosmo_x, cosmo_OmegaR+cosmo_OmegaNu,  'red', label="Relativistic")
plt.plot(cosmo_x, cosmo_OmegaB+cosmo_OmegaCDM, 'blue', label="Matter")
plt.plot(cosmo_x, cosmo_OmegaLambda, 'purple', label="$\Omega_{\Lambda}$")
#plt.plot(cosmo_x, 1/cosmo_Hp*cosmo_dHpdx)

# plt.axvline(x=cosmo_x[idx1], color='black', linestyle='--', linewidth=1)
# plt.axvline(x=cosmo_x[idx2], color='black', linestyle='--', linewidth=1)
plt.title("Omegas")
plt.legend()
plt.show()

### 1/H * dHdx
region1 = patches.Rectangle((cosmo_x[0], -1), cosmo_x[idx1]-cosmo_x[0], 2.5, color='red', alpha=0.2)
region2 = patches.Rectangle((cosmo_x[idx1], -1), cosmo_x[idx2] - cosmo_x[idx1], 2.5, color='blue', alpha=0.2)
region3 = patches.Rectangle((cosmo_x[idx2], -1), cosmo_x[-1] - cosmo_x[idx2], 2.5, color='purple', alpha=0.2)

# Add the patches to the plot
plt.gca().add_patch(region1)
plt.gca().add_patch(region2)
plt.gca().add_patch(region3)


plt.plot(cosmo_x, 1/cosmo_Hp*cosmo_dHpdx, 'black', label=r'$\frac{1}{\math{H}} \frac{d2\math{H}}{dx2}$')
plt.plot(cosmo_x, (1/cosmo_Hp)*cosmo_ddHpddx, label=r'$\frac{1}{\math{H}} \frac{d^2\math{H}}{dx^2}$')


plt.title("H")
#plt.legend()
plt.show()





exit()
""" SUPERNOVA FITTING"""
data = np.loadtxt("results_supernovafitting.txt")

converged_data = data[200:]

best_chi = np.argmin(converged_data[:, 0])


# Minimum chi^2 found 29.2799 0.701711 0.255027 0.0789514 




#print(converged_data[0,:])


OmegaM = converged_data[:,2]
OmegaK = converged_data[:,3]

OmegaB      = 0.05;
OmegaCDM    = 0.267;

OmegaLambda = 1 - (OmegaK + OmegaM) #+ OmegaB + OmegaCDM)

best_fit_params = converged_data[best_chi, :]


selected_data = converged_data[converged_data[:, 0] < best_fit_params[0] + 3.53]
OmegaM_selected = selected_data[:, 2]
OmegaLambda_selected = 1 - (selected_data[:, 3] + OmegaM_selected)

""" 1 sigma confidence region  """
plt.scatter(OmegaM_selected, OmegaLambda_selected) #*** fix std colorbar
plt.xlabel('$\Omega_M$')
plt.ylabel('$\Omega_{\lambda}$')
plt.title('1$\sigma$ Confidence Region')

cbar = plt.colorbar()
#cbar.set_label('$\Chi^2')

#plt.show()



exit()
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


