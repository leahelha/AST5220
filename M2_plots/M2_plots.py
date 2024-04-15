import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy import constants as const

"""
STRUCTURE OF recombination.txt
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = npts_rec_arrays;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";

    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";

    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";

    fp << sound_horizon_of_x(x) << " ";
    fp << Xe_of_x_Saha(x) << " ";

    fp << "\n";
  };
"""
data = np.loadtxt("./recombination.txt")

x = data[:, 0]
Xe = data[:, 1]
ne = data[:, 2]
tau = data[:, 3]
dtau = data[:, 4]
ddtau = data[:, 5]
gtilde = data[:, 6]
dgtilde = data[:, 7]
ddgtilde = data[:, 8]


print(f'Freeze out: Xe = {Xe[-1]}')


""" Finding the x time of last scattering """
# print(np.max(gtilde))

cosmo = np.loadtxt("cosmology.txt")
t_of_x = cosmo[ : , 11]

idx = np.where(gtilde==np.max(gtilde))
# print(x[idx])   # OUTPUT x = -6.98462

# SOLVING FOR z and t 
a_decoupled = np.exp(x[idx][0])

x_decoupled = x[idx][0]
z_decoupled = (1/a_decoupled) - 1
#t_decoupled = t_of_x[idx][0]

""" finding t from t_of_x in cosmo.cpp"""
print(f" x = {x_decoupled},  z = {z_decoupled}")

"""
Output from C++ code:
t_decoupled = 0.410948 Myr

From python code:
x = -6.98462,  z = 1078.8959792819974

"""

""" The recombination predicted by Saha"""

Xe_Saha = data[:, 10]
print(Xe_Saha)
#time of last Saha 
idx_Saha = np.where(Xe_Saha < 0.1)[0][0]
print(f"SAHA INDEX IS {idx_Saha}")
x_Saha = x[idx_Saha]  # -7.37028
a_Saha = np.exp(x[idx_Saha])
z_Saha = (1/a_Saha) - 1
print(f"x_Saha = {x_Saha}, z_Saha = {z_Saha}")

"""
Output from C++ code:
t_Saha 0.318387 Myr

From python code:
x_Saha = -7.14033, z_Saha = 1260.8447291590362

"""

""" Recombination time defining recombination as x when Xe = 0.1"""

idx_re = np.where(Xe < 0.1)
print(idx_re[0][0])

x_re = x[idx_re[0][0]]
a_re = np.exp(x_re)
z_re = (1/a_re) - 1

print(f"x_re = {x_re}, z_re = {z_re}")
"""
Output from C++ code:
t_re 0.408541 Myr

From python code:
x_re = -6.98822, z_re = 1082.790610938193

"""

""" PLOT OF TAU AND ITS DERIVATIVES """
plt.figure()
plt.xlim(-9, -4)
plt.ylim(10**(-5), 10**5, 10**1)
plt.yscale('log')
plt.plot(x, ddtau, color='firebrick', linestyle=':', label=r"$\tau''(x)$")
plt.plot(x, tau, color='sandybrown', linestyle='-', label=r"$\tau (x)$")
plt.plot(x, -dtau, color='darkorange', linestyle='--', label=r"-$\tau'(x)$")
plt.axvline(x=x_decoupled, color='magenta', linestyle='-', label="Last scattering", alpha=0.3)  # Adding vertical line at x_decoupled
plt.axvline(x=x_re, color='black', linestyle='dashdot', label="Recombination", alpha=0.4)  # Adding vertical line at x_re
plt.axhline(y=1, color='grey', linestyle='solid', alpha=0.3)  # Adding vertical line at x_Saha
plt.legend()
plt.xlabel("$x$", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
# plt.grid('True')
plt.savefig("./Figs/M2/Taus_vs_x.pdf")


""" PLOT OF G_TILDE AND DERIVATIVES """
plt.figure()  # Specifies the figure size for the third plot
plt.xlim(-7.8, -6)
plt.plot(x, ddgtilde/100, color='darkgreen', linestyle=':', label=r'$\frac{d^2\tilde{g}}{dx^2}(x) /100$')
plt.plot(x, dgtilde/10, color='olive', linestyle='--', label=r'$\frac{d\tilde{g}}{dx}(x)/ 10$')
plt.plot(x, gtilde, color='yellowgreen', linestyle='-', label=r'$\tilde{g}(x)$')
plt.axvline(x=x_decoupled, color='magenta', linestyle='-', label="Last scattering", alpha=0.3)  # Adding vertical line at x_decoupled
plt.axvline(x=x_re, color='black', linestyle='dashdot', label="Recombination", alpha=0.4)  # Adding vertical line at x_re
# plt.axvline(x=x_Saha, color='purple', linestyle='dashdot', label="Saha recombination", alpha=0.3)  # Adding vertical line at x_Saha
plt.xlabel("$x$", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.legend()
plt.savefig("./Figs/M2/Gtilde_all.pdf")


""" PLOT OF X_e """
plt.figure()
plt.xlim(-9,-3)
plt.yscale("log")
plt.plot(x, Xe, color='darkred')
plt.plot(x, Xe_Saha, color='salmon', linestyle='--', label="Saha")
plt.axvline(x=x_decoupled, color='magenta', linestyle='-', label="Last scattering", alpha=0.3)  # Adding vertical line at x_decoupled
plt.axvline(x=x_re, color='black', linestyle='dashdot', label="Recombination", alpha=0.4)  # Adding vertical line at x_re
plt.axvline(x=x_Saha, color='purple', linestyle='dashdot', label="Saha recombination", alpha=0.3)  # Adding vertical line at x_Saha
# plt.axhline(y=, color='dodgerblue', linestyle='dashdot', label="Freeze out")  # Adding horizontal line at Xe=0
plt.xlabel("$x$", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.ylabel("$X_e$", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.grid('True')
plt.legend()
plt.savefig("./Figs/M2/Xe_vs_x.pdf")
plt.show()


""" 
plt.figure(figsize=(10, 6))  # Specifies the figure size for the first plot
plt.xlim(-12, 0)
plt.plot(x, gtilde, label=r'$\tilde{g}(x)$')
plt.xlabel('x')
plt.ylabel(r'$\tilde{g}(x)$')
plt.legend()
# plt.title(r'Plot of $\tilde{g}$ vs. x')
plt.savefig("./Figs/M2/Gtilde_vs_x.pdf")
# plt.show()  # Show the first plot

# Plot dg_tilde/dx vs x
plt.figure(figsize=(10, 6))  # Specifies the figure size for the second plot
plt.xlim(-12, 0)
plt.plot(x, dgtilde, label=r'$\frac{d\tilde{g}}{dx}(x)$')
plt.xlabel('x')
plt.ylabel(r'$\frac{d\tilde{g}}{dx}(x)$')
# plt.title(r'Plot of $\frac{d\tilde{g}}{dx}$ vs. x')
plt.legend()
plt.savefig("./Figs/M2/Dgtilde_vs_x.pdf")
# plt.show()  # Show the second plot

# Plot d^2g_tilde/dx^2 vs x
plt.figure(figsize=(10, 6))  # Specifies the figure size for the third plot
plt.xlim(-12, 0)
plt.plot(x, ddgtilde, label=r'$\frac{d^2\tilde{g}}{dx^2}(x)$')
plt.xlabel('x')
plt.ylabel(r'$\frac{d^2\tilde{g}}{dx^2}(x)$')
# plt.title(r'Plot of $\frac{d^2\tilde{g}}{dx^2}$ vs. x')
plt.legend()
plt.savefig("./Figs/M2/DdGtilde_vs_x.pdf")
"""