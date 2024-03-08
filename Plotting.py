import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy import constants as const

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
    fp << ddHpddx_of_x(x)        << " "; //***
    fp << t_of_x(x)        << " "; //***
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

"""
Gyr = 1/(60*60*24*365*1e9) # from s to Gyr
Mpc = 3.24*10**(-24) # from cm to Mpc
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




# ROUGHLY ILLUSTRATING THE DIFFERENT REGIMES IN THE PLOT
border_idx1 = np.where((np.abs((cosmo_OmegaR+cosmo_OmegaNu)-(cosmo_OmegaB+cosmo_OmegaCDM)))<0.085)[0]
border_idx2 = np.where((np.abs((cosmo_OmegaB+cosmo_OmegaCDM)-cosmo_OmegaLambda))<0.085)[0]

idx1 = border_idx1[0]
idx2 = border_idx2[-1]

plt.figure()
# Patches for filling between lines
region1 = patches.Rectangle((cosmo_x[0], 0), cosmo_x[idx1]-cosmo_x[0], 1, color='red', alpha=0.2)
region2 = patches.Rectangle((cosmo_x[idx1], 0), cosmo_x[idx2] - cosmo_x[idx1], 1, color='blue', alpha=0.2)
region3 = patches.Rectangle((cosmo_x[idx2], 0), cosmo_x[-1] - cosmo_x[idx2], 1, color='purple', alpha=0.2)

# Adding the patches to the plot
plt.gca().add_patch(region1)
plt.gca().add_patch(region2)
plt.gca().add_patch(region3)


""" Plot of Omegas """

plt.plot(cosmo_x, cosmo_OmegaR+cosmo_OmegaNu,  'red', label=r"$\Omega_{Relativistic} = \Omega_{\gamma} + \Omega_{\nu}$")
plt.plot(cosmo_x, cosmo_OmegaB+cosmo_OmegaCDM, 'blue', label=r"$\Omega_{Matter} = \Omega_{b} + \Omega_{CDM}$")
plt.plot(cosmo_x, cosmo_OmegaLambda, 'purple', label="$\Omega_{\Lambda}$")
#plt.plot(cosmo_x, 1/cosmo_Hp*cosmo_dHpdx)

# plt.axvline(x=cosmo_x[idx1], color='black', linestyle='--', linewidth=1)
# plt.axvline(x=cosmo_x[idx2], color='black', linestyle='--', linewidth=1)
plt.title("$\Omega_i(x)$")
plt.xlabel("x")
plt.legend()
plt.savefig("Figs/Omegas.pdf")



""" Hprime *its derivatives plots """

plt.figure()
### 1/H * dHdx
region1 = patches.Rectangle((cosmo_x[0], -1), cosmo_x[idx1]-cosmo_x[0], 2.5, color='red', alpha=0.2)
region2 = patches.Rectangle((cosmo_x[idx1], -1), cosmo_x[idx2] - cosmo_x[idx1], 2.5, color='blue', alpha=0.2)
region3 = patches.Rectangle((cosmo_x[idx2], -1), cosmo_x[-1] - cosmo_x[idx2], 2.5, color='purple', alpha=0.2)

# Add the patches to the plot
plt.gca().add_patch(region1)
plt.gca().add_patch(region2)
plt.gca().add_patch(region3)


plt.plot(cosmo_x, 1/cosmo_Hp*cosmo_dHpdx, 'black', label=r'$\frac{1}{\mathcal{H}(x)} \frac{d\mathcal{H}(x)}{dx}$')
plt.plot(cosmo_x, (1/cosmo_Hp)*cosmo_ddHpddx, label=r'$\frac{1}{\mathcal{H}(x)} \frac{d^2\mathcal{H}(x)}{dx^2}$')


#plt.title("$\math{H}(x)$")
plt.legend()
plt.xlabel("x")
plt.savefig("Figs/Hp_checks.pdf")
#plt.show()

""" Plot of Hprime(x)*eta(x)/c """
plt.figure()

region1 = patches.Rectangle((cosmo_x[0], 0), cosmo_x[idx1]-cosmo_x[0], 3.5, color='red', alpha=0.2)
region2 = patches.Rectangle((cosmo_x[idx1], 0), cosmo_x[idx2] - cosmo_x[idx1], 3.5, color='blue', alpha=0.2)
region3 = patches.Rectangle((cosmo_x[idx2], 0), cosmo_x[-1] - cosmo_x[idx2], 3.5, color='purple', alpha=0.2)

# Add the patches to the plot
plt.gca().add_patch(region1)
plt.gca().add_patch(region2)
plt.gca().add_patch(region3)

#print(cosmo_x[-23])  # when x = 0
plt.plot(cosmo_x[:-22], cosmo_eta_of_x[:-22]*cosmo_Hp[:-22]/const.c, label=r'$\frac{\eta(x)\mathcal{H}(x)}{c}$')
plt.legend()
#plt.xlim(np.log(1/1+1089), 0)
plt.xlabel("x")
plt.savefig("Figs/Hp_eta_checks.pdf")


""" Plot of cosmic time and conformal time """
plt.figure()
plt.plot(cosmo_x[1:], cosmo_t_of_x[1:]*Gyr, label="t") #Had to skip first index, bc it was acting up
plt.plot(cosmo_x[1:], cosmo_eta_of_x[1:]*Gyr/const.c, label="$\eta(x)/c$")
plt.axhline(y=13.8, color='black', linestyle='--', label='13.8 Gyr', alpha=0.3)

plt.yscale('log')
#plt.xlim(np.log(1/1+1089), 0)
plt.xlabel("x")
plt.ylabel("Gyr")
plt.legend()
plt.savefig("Figs/cosmic_time_and_conformal_time.pdf")


""" Hubble factor Hprime(x) """
print(f"Hubble factor {cosmo_Hp[-23]*(100/(Mpc*1000))}")

plt.figure()
plt.plot(cosmo_x, cosmo_Hp*(100/(Mpc*1000)))

plt.yscale('log')
plt.xlim(-12, 0.1)
plt.title("$\mathcal{H}(x)$")
plt.xlabel("x")
plt.ylabel("Mpc")
plt.savefig("Figs/Hubble_factor.pdf")
#plt.show()



