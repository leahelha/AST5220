import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy import constants as const

Gyr = 1/(60*60*24*365*1e9) # from s to Gyr
Mpc = 3.24*10**(-23) # from m to Mpc
Gpc = 3.24*10**(-25)


cosmo = np.loadtxt("./cosmology.txt")
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

data_rec = np.loadtxt("./recombination.txt")

x = data_rec[:, 0]
Xe = data_rec[:, 1]
ne = data_rec[:, 2]
tau = data_rec[:, 3]
dtau = data_rec[:, 4]
ddtau = data_rec[:, 5]
gtilde = data_rec[:, 6]
dgtilde = data_rec[:, 7]
ddgtilde = data_rec[:, 8]




"""                                       PERTURBATIONS                                         """

"""
         void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);

  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";

    fp << get_delta_cdm(x,k)   << " ";
    fp << get_v_cdm(x,k)       << " ";
    fp << get_delta_b(x,k)       << " ";
    fp << get_v_b(x,k)        << " ";

    // fp << get_Source_T(x,k)  << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    
    fp << "\n";
  };

"""
# Load the data from files
pert_k_S = np.loadtxt('./perturbations_k0.001.txt')
pert_k_M = np.loadtxt('./perturbations_k0.01.txt')
pert_k_L = np.loadtxt('./perturbations_k0.1.txt')

N = len(pert_k_L)


k_list = [pert_k_S, pert_k_M, pert_k_L]
k_values = [0.001, 0.01, 0.1]

x = np.zeros((len(k_list), N))

Theta_0 = np.zeros((len(k_list), N))
Theta_1 = np.zeros((len(k_list), N))
Theta_2 = np.zeros((len(k_list), N))

Phi = np.zeros((len(k_list), N))
Psi = np.zeros((len(k_list), N))

delta_cdm = np.zeros((len(k_list), N))
delta_b = np.zeros((len(k_list), N))
v_cdm = np.zeros((len(k_list), N))
v_b = np.zeros((len(k_list), N))

Source_T = np.zeros((len(k_list), N))
Source_T_5 = np.zeros((len(k_list), N))
Source_T_50 = np.zeros((len(k_list), N))
Source_T_500 = np.zeros((len(k_list), N))

i = 0
for k in k_list:
    x[i] = k[:, 0]

    Theta_0[i] = k[:, 1]
    Theta_1[i] = k[:, 2]
    Theta_2[i] = k[:, 3]


    Phi[i] = k[:, 4]
    Psi[i] = k[:, 5]

    delta_cdm[i] = k[:, 6]
    print(f"LAST DELTA_CDM is {delta_cdm[i][-1]}")

    v_cdm[i] = k[:, 7]

    delta_b[i] = k[:, 8]
    v_b[i] = k[:, 9]
    # Pi[i] = k[:, 5] # BRUKES IKKE!!!

    # Source_T[i] = k[:, 6]
    # Source_T_5[i] = k[:, 7]
    # Source_T_50[i] = k[:, 8]
    # Source_T_500[i] = k[:, 9]

    i += 1






for k in range(len(k_list)):
    plt.plot(x[k], Theta_0[k], label=f"k = {k_values[k]}")
plt.legend()
plt.title('\Theta_0')
plt.show()

for k in range(len(k_list)):
    plt.plot(x[k], Theta_1[k], label=f"k = {k_values[k]}")
plt.legend()
plt.title('\Theta_1')
plt.show()

# for k in range(len(k_list)):
#     plt.plot(x[k], 4*Theta_0[k], label=f"k = {k_values[k]}")
# plt.legend()
# plt.title('\delta_gamma = 4\Theta_0')
# plt.show()

# for k in range(len(k_list)):
#     plt.plot(x[k], Psi[k]+Phi[k], label=f"k = {k_values[k]}")
# plt.legend()
# plt.title('Psi + Phi')
# plt.show()

for k in range(len(k_list)):
    plt.plot(x[k], Phi[k], label=f"k = {k_values[k]}")
plt.legend()
plt.title('Phi')
plt.show()

# for k in range(len(k_list)):
#     plt.plot(x[k], Psi[k], label=f"k = {k_values[k]}")
# plt.legend()
# plt.title('Psi')
# plt.show()

for k in range(len(k_list)):
    plt.plot(x[k], delta_cdm[k], label=f"delta_cdm k = {k_values[k]}")
    plt.plot(x[k], delta_b[k], '--', label=f"delta_b k = {k_values[k]}")
plt.yscale('log')
plt.legend()
plt.title('delta_cdm and delta_b')
plt.show()

for k in range(len(k_list)):
    plt.plot(x[k], v_cdm[k], label=f"v_cdm k = {k_values[k]}")
    plt.plot(x[k], v_b[k], '--', label=f"v_b k = {k_values[k]}")
plt.yscale('log')
plt.legend()
plt.title('v_cdm and v_b')
plt.show()




"""

# TODO

MAKE PLOTS OF THETAs for all three ks

make plots of the other stuff, compare with Callins


"""

""" 

# Dictionary to hold the data
perturbations = {
    'pert_k_S': pert_k_S,
    'pert_k_M': pert_k_M,
    'pert_k_L': pert_k_L
}

# Dictionary to store results
results = {}

# Process each key and data in the perturbations dictionary
for k, data in perturbations.items():
    results[f'{k}_Theta_0'] = data[:, 0]
    results[f'{k}_Theta_1'] = data[:, 1]
    results[f'{k}_Theta_2'] = data[:, 2]

    results[f'{k}_Phi'] = data[:, 3]
    results[f'{k}_Psi'] = data[:, 4]
    results[f'{k}_Pi'] = data[:, 5]  # Comment if not used

    # results[f'{k}_Source_T'] = data[:, 6]
    # results[f'{k}_Source_T_5'] = data[:, 7]
    # results[f'{k}_Source_T_50'] = data[:, 8]
    # results[f'{k}_Source_T_500'] = data[:, 9]

# Now you can access any variable like so:
# print(results['pert_k_S_Theta_0'])
"""