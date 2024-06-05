import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy import constants as const
from scipy.stats import norm

Gyr = 1/(60*60*24*365*1e9) # from s to Gyr
Mpc = 3.24*10**(-23) # from m to Mpc
Gpc = 3.24*10**(-26) # from m to Gpc

cosmo = np.loadtxt("cosmology.txt")
cosmo_best = np.loadtxt("best_params_cosmology.txt")
print(f"Shape of cosmo = {np.shape(cosmo)}")

"""
                              chi2         h        OmegaM      OmegaK    Acceptrate    
Best params from min chi^2 [29.2799     0.701711   0.255027   0.0789514]
"""




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

""" For best fit cosmological parameters """
cosmo_dL_best = cosmo_best[:, 12]
cosmo_x_best = cosmo_best[:,0]


""" SUPERNOVA FITTING"""
data = np.loadtxt("results_supernovafitting.txt")

converged_data = data[200:]
print(f"The shape of the converged data is {np.shape(converged_data)}\n")
best_chi = np.argmin(converged_data[:, 0])


best_fit_params = converged_data[best_chi, :]
print(f"Best params from min chi^2 {best_fit_params}\n")


selected_data = converged_data[converged_data[:, 0] < (best_fit_params[0] + 3.53)] #Data selected within 1sigma of the best fit 1 sigma
sig2_data = converged_data[converged_data[:, 0]  < (best_fit_params[0] + 8.02)] #Data selected within 2sigma of the best fit
sig3_data = converged_data[converged_data[:, 0] < (best_fit_params[0] + 14.16)] #Data selected within 3sigma of the best fit

chi2 = selected_data[:, 0]
h_selected = selected_data[:, 1]
OmegaM_selected = selected_data[:, 2]
OmegaK_selected = selected_data[:, 3]

OmegaLambda_selected = 1 - (OmegaK_selected + OmegaM_selected) 

OmegaM_s2 = sig2_data[:, 2]
OmegaK_s2 = sig2_data[:, 3]
OmegaLambda_s2 = 1 - (OmegaK_s2 + OmegaM_s2)

OmegaM_s3 = sig3_data[:, 2]
OmegaK_s3 = sig3_data[:, 3]
OmegaLambda_s3 = 1 - (OmegaK_s3 + OmegaM_s3)

""" dL plot with supernova best fit, fiducial cosmology and betouli observations """

betoule = np.loadtxt("Betoule_supernova.txt")


z_cosmo = np.exp(-cosmo_x)-1
# print(f"z cosmo = {z_cosmo}, \n and cosmo_dL * Gpc = {cosmo_dL*Gpc}")
z_obs = betoule[:, 0]

z_best = np.exp(-cosmo_x_best)-1


dL_obs = betoule[:, 1]
error_obs = betoule[:, 2]

print(f'z_obs len = {z_obs} \n and z_cosmo = {z_cosmo[-31:]}')

"""
X    Du har Gpc = 3.24*10**(-25) men 1m er 3.24*10^(-26) Gpc så det er en faktor av 10 feil her.  DONE

Et annet problem her er at z-verdiene dine går fra rundt 0 til 10^9, mens du bare er interessert i z-verdier fra rundt 0 til rundt 1. 
Ta å lag en separat rutine der du outputter dL dataene der du velger det z-området du vil ha og så printer dette. 
Hvis du vil bruke dataene fra output() så må du sørge for å ha enormt mange punkter for å sørge for nok punkter i det intervallet 
du vil ha og så må du sette x-range til å være det du er interessert i.

"""


""" Plot of luminosity distance for fiducial cosmology, observed sn data and our best fit results """
# *** THIS IS NOT RIGHT
plt.figure()
plt.plot(z_cosmo, cosmo_dL*Gpc/z_cosmo, color="blue", label="Fiducial cosmology")
plt.plot(z_best, cosmo_dL_best*Gpc/z_best, color="orange", label="Best fit from MCMC")
plt.errorbar(z_obs, dL_obs/z_obs, yerr=error_obs/z_obs, fmt='o', color='red', ecolor='red', capsize=2, ms=2, label="Observed data")
# plt.title("$d_L$", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.xlabel('z', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.ylabel('Gpc', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.xscale('log')
plt.xlim(0.005, 1.45)
plt.ylim(3.5, 8)
plt.legend()
plt.savefig("Figs/sn_dL_plots.pdf")



""" Confidence region 1sig, 2sig and 3sig """
plt.figure()
best_lambda = 1 - (0.0789514  + 0.255027) 
# plt.scatter(OmegaM_s3,  OmegaLambda_s3, label = "$3\sigma$")
plt.scatter(OmegaM_s2, OmegaLambda_s2, color='maroon', s=3, label="$2\sigma$")
scatter = plt.scatter(OmegaM_selected, OmegaLambda_selected, c=chi2, cmap='viridis', s=3, label="$1\sigma$")
plt.colorbar(scatter, label=r"$\chi^2$")
plt.plot((0,1), (1,0), color='black', linestyle = '--', label='Flat universe')
plt.scatter(0.255027, best_lambda, color='orange', marker='*', label='Best fit')
plt.xlim(0.0, 0.6)
plt.legend()
plt.xlabel('$\Omega_{M0}$', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.ylabel('$\Omega_{\Lambda0}$', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
# plt.title('Confidence Region of $\Omega_{M0}$ and $\Omega_{\Lambda0}$', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.savefig("Figs/sn_Confidence_region.pdf")

std_OmegaM = np.std(OmegaM_selected)
std_OmegaLambda = np.std(OmegaLambda_selected)
std_OmegaK= np.std(OmegaK_selected)
std_h= np.std(h_selected)
print("Standard deviation of OmegaM:", std_OmegaM)
print("Standard deviation of OmegaLambda:", std_OmegaLambda)


""" Histogram of accepted H """
# H0 = 100*h # km /s / Mpc H0 over h
# best fit h = 0.701711

h_2sig = sig2_data[:, 1]
h_3sig = sig3_data[:, 1]

hist_sig = np.std(h_selected)
mu = np.mean(h_selected)#0.701711
best_fit = 0.701711

x_array = np.linspace(min(h_selected), max(h_selected), 1000)


plt.figure()
# plt.hist(h_3sig, bins=100, alpha=0.5, label='$H_0$', edgecolor='white', linewidth=0.5) 
# plt.hist(h_2sig, bins=100, alpha=0.5, label='$H_0$', edgecolor='white', linewidth=0.5) 
plt.hist(h_selected, color='teal', bins=100, alpha=0.5, edgecolor='white', linewidth=0.5, density=True) 
plt.axvline(x=0.701711, color='black', linestyle='--', label=f'Best $H_0$ = {best_fit:.3f}')
pdf = norm.pdf(x_array, mu, hist_sig)
plt.plot(x_array, pdf, color='navy')
plt.text(0.05, 0.95, f'$\mu$=0.702 \n $\sigma={hist_sig:.3f}$', transform=plt.gca().transAxes, 
         verticalalignment='top', bbox=dict(boxstyle="square,pad=0.3", facecolor='white', alpha=0.5, edgecolor='none'))
plt.xlabel('100 km/s/Mpc', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
# plt.ylabel('Frequency normalized')
plt.legend()
# plt.title('Posterior pdf of $H_0$', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.savefig("Figs/sn_Histogram_of_H_parameters.pdf")


""" Plot the histogram and omegas distribution """


plt.figure()
OmegaM_mean = np.mean(OmegaM_selected)
best_M = 0.255027
OmegaM_std = np.std(OmegaM_selected)

OmegaM_x = np.linspace(np.min(OmegaM_selected), np.max(OmegaM_selected), 1000)
pdf_M = norm.pdf(OmegaM_x, OmegaM_mean, OmegaM_std)


plt.hist(OmegaM_selected, color='teal', bins=100, alpha=0.5,  edgecolor='white', linewidth=0.5, density=True)
plt.axvline(x=0.255027, color='black', linestyle='--', label=r'Best $\Omega_{M0}$ = 0.255')
plt.plot(OmegaM_x, pdf_M, color='navy')
plt.text(0.05, 0.95, f'$\mu$={OmegaM_mean:.3f} \n $\sigma={OmegaM_std:.3f}$', transform=plt.gca().transAxes, 
         verticalalignment='top', bbox=dict(boxstyle="square,pad=0.3", facecolor='white', alpha=0.5, edgecolor='none') )

# plt.title('Posterior pdf of $\Omega_{M0}$', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.xlabel('$\Omega_{M0}$', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.legend()
plt.savefig("Figs/sn_Histogram_of_omegaM_Gaussian.pdf")


plt.figure()
OmegaK_mean = np.mean(OmegaK_selected)
best_K = 0.0789514
OmegaK_std = np.std(OmegaK_selected)
OmegaK_x = np.linspace(np.min(OmegaK_selected), np.max(OmegaK_selected), 1000)
pdf_K = norm.pdf(OmegaK_x, OmegaK_mean, OmegaK_std)

plt.hist(OmegaK_selected, color='teal', bins=100, alpha=0.5, edgecolor='white', linewidth=0.5, density=True)
plt.axvline(x=0.0789514, color='black', linestyle='--', label=r'Best $\Omega_{K0}$ = 0.0790')
plt.plot(OmegaK_x, pdf_K, color='navy')
plt.xlabel('$\Omega_{K0}$', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})

plt.text(0.05, 0.95, f'$\mu$={OmegaK_mean:.3f}\n$\sigma$={OmegaK_std:.3f}', transform=plt.gca().transAxes, 
         verticalalignment='top', bbox=dict(boxstyle="square,pad=0.3", facecolor='white', alpha=0.5, edgecolor='none'))

plt.legend()
# plt.title('Posterior pdf of $\Omega_{K0}$', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.savefig("Figs/sn_Histogram_of_omegaK_Gaussian.pdf")



""" Goodness of fit evaluation """
# plt.figure()
# N = 31#len(cosmo_dL)
# fit_check = chi2/N

# chi_x = np.linspace(np.min(fit_check), np.max(fit_check), 1000)

# chi_mean = np.mean(fit_check)
# chi_std = np.std(fit_check)
# pdf_chi = norm.pdf(chi_x, chi_mean, chi_std)


# plt.hist(fit_check, color='teal', bins=100, alpha=0.5,  edgecolor='white', linewidth=0.5, density=True)
# plt.plot(chi_x, pdf_chi, color='navy')
# plt.xlabel(r'$\chi ^2/N$', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
# plt.text(0.05, 0.95, f'$\mu$={chi_mean:.3f} \n $\sigma={chi_std:.3f}$', transform=plt.gca().transAxes, 
#          verticalalignment='top', bbox=dict(boxstyle="square,pad=0.3", facecolor='white', alpha=0.5, edgecolor='none'))
# plt.savefig('Figs/Goodness_of_fit.pdf')
# plt.show()




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