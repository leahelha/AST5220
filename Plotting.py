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



