import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy import constants as const

"""
STRUCTURE OF recombination.txt
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
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
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);

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

plt.figure()
plt.yscale('log')
plt.plot(x, tau, label=r"$\tau (x)$")
plt.plot(x, dtau, label=r"$\tau'(x)$")
plt.plot(x, dtau, label=r"$\tau''(x)$")
plt.savefig("Taus_vs_x.pdf")
plt.legend()
plt.show()


# Plot g_tilde vs x
plt.figure(figsize=(10, 6))  # Specifies the figure size for the first plot
plt.plot(x, gtilde, label=r'$\tilde{g}(x)$')
plt.xlabel('x')
plt.ylabel(r'$\tilde{g}(x)$')
# plt.title(r'Plot of $\tilde{g}$ vs. x')
plt.legend()
plt.savefig("Gtilde_vs_x.pdf")
plt.show()  # Show the first plot

# Plot dg_tilde/dx vs x
plt.figure(figsize=(10, 6))  # Specifies the figure size for the second plot
plt.plot(x, dgtilde, label=r'$\frac{d\tilde{g}}{dx}(x)$')
plt.xlabel('x')
plt.ylabel(r'$\frac{d\tilde{g}}{dx}(x)$')
# plt.title(r'Plot of $\frac{d\tilde{g}}{dx}$ vs. x')
plt.legend()
plt.savefig("Dgtilde_vs_x.pdf")
plt.show()  # Show the second plot

# Plot d^2g_tilde/dx^2 vs x
plt.figure(figsize=(10, 6))  # Specifies the figure size for the third plot
plt.plot(x, ddgtilde, label=r'$\frac{d^2\tilde{g}}{dx^2}(x)$')
plt.xlabel('x')
plt.ylabel(r'$\frac{d^2\tilde{g}}{dx^2}(x)$')
# plt.title(r'Plot of $\frac{d^2\tilde{g}}{dx^2}$ vs. x')
plt.legend()
plt.savefig("DdGtilde_vs_x.pdf")
plt.show()  # Show the third plot
