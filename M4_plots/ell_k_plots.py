import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy import constants as const

h = 0.67
Gyr = 1/(60*60*24*365*1e9) # from s to Gyr
Mpc = 3.24*10**(-23) # 1 m in Mpc
Gpc = 3.24*10**(-25) 
c = const.c

cosmo = np.loadtxt("cosmology.txt")
cosmo_x = cosmo[:,0]
cosmo_eta_of_x = cosmo[:,1]
cosmo_Hp_of_x = cosmo[:,2]
TempCMB = cosmo[:, 14]

TempCMB0 = 2.7255#TempCMB[np.argmin(abs(cosmo_x))]

eta_0 = cosmo_eta_of_x[np.argmin(abs(cosmo_x))]
print(f'Eta_0 = {eta_0}, x = {cosmo_x[np.argmin(abs(cosmo_x))]}')

data_cell = np.loadtxt("./cells.txt")
data_SW = np.loadtxt("./cells_SW.txt")
data_ISW = np.loadtxt("./cells_ISW.txt")
data_Doppler = np.loadtxt("./cells_Doppler.txt")
data_Polarization = np.loadtxt("./cells_Polarization.txt")

l = data_cell[:, 0]
Cell = data_cell[:, 1] # Cell is normalized in text
Cell_SW = data_SW[:, 1]
Cell_ISW = data_ISW[:, 1]
Cell_Doppler = data_Doppler[:, 1]
Cell_Polarization = data_Polarization[:, 1]

planck_data = np.loadtxt('./M4_plots/planck_cell_low.txt')

ell_planck = planck_data[:,0]
C_ell_planck = planck_data[:,1]
upper_error = planck_data[:,2]
lower_error = planck_data[:,3]

error = [lower_error, upper_error]

# Angular powerspectrum
plt.figure()
# plt.plot(l, Cell*(l*(l+1))/(2*np.pi)*(10**6*TempCMB0)**2 )
lnwidth = 0.5
plt.plot(l, Cell, color='blue')
plt.plot(l, Cell_SW, linestyle='-', linewidth=lnwidth, label='SW')
plt.plot(l, Cell_ISW, linestyle='-',linewidth=lnwidth, label='ISW')
plt.plot(l, Cell_Doppler,linestyle='-', linewidth=lnwidth, label='Doppler')
plt.plot(l, Cell_Polarization, linestyle='-',linewidth=lnwidth, label='Polarization')
plt.errorbar(ell_planck, C_ell_planck, elinewidth=0.5, yerr=error, fmt='o', color='red', ecolor='red', capsize=1, ms=1, label="Planck data")
# plt.xlim(2, l[-1])
plt.xscale("log")
plt.ylabel(r"$C_{\ell}\cdot\frac{\ell\cdot(\ell + 1)}{2\pi} 10^6 T_{CMB}^2$", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.xlabel("$\ell$", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.legend()
plt.savefig('./Figs/M4/C_ell.pdf')



""" Make plot of Theta_100(k) from l=100"""
"""
 std::ofstream fp(filename.c_str());

  auto print_data = [&] (const double k) {

    fp << k  << " ";
    fp << get_matter_power_spectrum(0, k) << " ";
    fp << thetaT_ell_of_k_spline[0](k)    << " ";
    fp << thetaT_ell_of_k_spline[5](k)    << " ";
    fp << thetaT_ell_of_k_spline[15](k)    << " ";
    fp << thetaT_ell_of_k_spline[19](k)    << " ";
    fp << thetaT_ell_of_k_spline[25](k)    << " ";
    fp << thetaT_ell_of_k_spline[29](k)    << " ";

    fp << "\n";
  };

  Vector ells{ 
        2*,    3,    4,    5,    6,    7,    8,    10,   12,   15,   
        20,   25,   30,   40,   50,   60*,   70,   80,   90,   100*,  
        120,  140,  160,  180,  200,  225*,  250,  275,  300,  350*,  
        400,  450,  500,  550,  600,  650,  700,  750,  800,  850,  
        900,  950,  1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 
        1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 
        1900, 1950, 2000};

"""


# Matter-power spectrum

data_2 = np.loadtxt("./matter_transfer.txt")
 
k = data_2[:, 0]
pofk = data_2[:, 1]

theta_ell_2 = data_2[:, 2]
theta_ell_7 = data_2[:, 8]
theta_ell_60 = data_2[:, 3] # idx 15
theta_ell_100 = data_2[:, 4] # 19
theta_ell_500 = data_2[:, 5] # 32
theta_ell_1000 = data_2[:, 6] #42
theta_ell_2000 = data_2[:, 7]
# [2, 60, 100, 500, 1000]

# from M1: rm_time: x = -8.65778, z = 5753.744937399699, t = 0.0232026 Myr
rm_time = -8.65778
x_eq = cosmo_x[np.argmin(abs(cosmo_x-rm_time))]

print(x_eq)
k_eq = cosmo_Hp_of_x[np.where(cosmo_x==x_eq)]/c.value

plt.figure()
plt.plot(k/h /Mpc, pofk*(Mpc**3)*(h**3))
plt.axvline(k_eq/h/Mpc, color='black', linestyle='--', label=r'$k_{eq}$')
plt.xlabel('k [1/h]', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.ylabel(r'$P(k)(Mpc/h)^3$', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig('./Figs/M4/Pk.pdf')


# # Theta_l

plt.figure()

# ells = [ 2, 7, 60, 100, 225, 350 ]


linewidth = 1
# plt.plot(k*eta_0, theta_ell_2, linewidth=linewidth, label = r'$\ell$ = 2')
plt.plot(k*eta_0, theta_ell_7, linewidth=linewidth, label = r'$\ell$ = 7')
plt.plot(k*eta_0, theta_ell_60, linewidth=linewidth, label = r'$\ell$ = 60')
plt.plot(k*eta_0, theta_ell_100,linewidth=linewidth, label = r'$\ell$ = 100')

plt.plot(k*eta_0, theta_ell_500, linewidth=linewidth, label = r'$\ell$ = 500')
plt.plot(k*eta_0, theta_ell_1000, linewidth=linewidth, label = r'$\ell$ = 1000')
# plt.plot(k*eta_0, theta_ell_2000, linewidth=linewidth, label = r'$\ell$ = 2000')
plt.xlim(0, 1500)
plt.ylim(-0.005, 0.02)
plt.ylabel(r'$\Theta_{\ell}(k)$', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.xlabel("$k\eta_0$", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.legend()
plt.savefig('./Figs/M4/Theta_ell_k.pdf')



# the C_l integrand, abs(theta_ell)**2/k
plt.figure()

"""
Hvis du bruker log-y så er det lurt å sette en minimums verdi for y aksen så den ikke løper helt ned til 10^-300 og du ikke ser noe i plottet.
For x så er det en mulighet å plotte k*eta0 siden Theta_ell(k) ~ const * j_ell(k eta0) i den enkleste approximajonen for SW leddet og vi vet 
at j_ell(x) peaker rundt ell ~ x og starter å svinge etter det så er kanskje lettere å se om alt se om ting ser riktig ut med denne x aksen.

"""
ells = [ 2, 7, 60, 100, 225, 350 ]

scale_2 = 2*(2+1)
scale_7 =  7*(7+1)
scale_60 = 60*(60+1)
scale_100 = 100*(100+1)
scale_500 = 500*(500+1)
scale_1000 = 1000*(1000+1)
scale_2000 =  2000*(2000+1)

print(np.where(abs(theta_ell_500)**2/k == np.max(abs(theta_ell_500)**2/k))[0][0])

linewidth = 0.5
# plt.plot(k*eta_0, abs(theta_ell_2)**2/k, linewidth=linewidth, label = r'$\ell$ = 2')
plt.plot(k[1:]*eta_0, scale_7*(abs(theta_ell_7[1:])**2)/k[1:], linewidth=linewidth, color='blue', label = r'$\ell$ = 7')
plt.axvline(k[ np.where((abs(theta_ell_7)**2)/k == np.max((abs(theta_ell_7)**2)/k))[0][0] ]*eta_0, linestyle='--', linewidth = 1, color='blue', label='')

plt.plot(k[1:]*eta_0, scale_60*(abs(theta_ell_60[1:])**2)/k[1:], linewidth=linewidth, label = r'$\ell$ = 60', color='orange',)
plt.axvline(k[ np.where((abs(theta_ell_60)**2)/k == np.max((abs(theta_ell_60)**2)/k))[0][0] ]*eta_0, linestyle='--', linewidth = 1, color='orange', label='')

plt.plot(k[1:]*eta_0, scale_100* (abs(theta_ell_100[1:])**2)/k[1:],linewidth=linewidth, label = r'$\ell$ = 100', color='gold')
plt.axvline(k[    np.where((abs(theta_ell_100)**2)/k == np.max((abs(theta_ell_100)**2)/k))[0][0]    ]*eta_0, linestyle='--', linewidth = 1, color='gold', label='')

plt.plot(k[1:]*eta_0, scale_500*(abs(theta_ell_500[1:])**2)/k[1:], linewidth=linewidth, label = r'$\ell$ = 500', color='olivedrab',)
plt.axvline(k[    np.where((abs(theta_ell_500)**2)/k == np.max((abs(theta_ell_500)**2)/k))[0][0]    ]*eta_0, linestyle='--', linewidth = 1, color='olivedrab', label='')

plt.plot(k[1:]*eta_0, scale_1000*(abs(theta_ell_1000[1:])**2)/k[1:], linewidth=linewidth, label = r'$\ell$ = 1000', color='magenta',)
plt.axvline(k[ np.where((abs(theta_ell_1000)**2)/k == np.max((abs(theta_ell_1000)**2)/k))[0][0]]*eta_0, linestyle='--', linewidth = 1, color='magenta', label='')

# plt.plot(k[1:]*eta_0, scale_2000*(abs(theta_ell_2000[1:])**2)/k[1:], linewidth=linewidth, label = r'$\ell$ = 2000', color='cyan')
# plt.axvline(k[ np.where((abs(theta_ell_2000)**2)/k == np.max((abs(theta_ell_2000)**2)/k))[0][0]]*eta_0, linestyle='--', linewidth = 1, color='cyan', label='')

plt.ylim(0, 1.5e24)
plt.xlim(-0.01,1500)
# plt.yscale('log')
plt.ylabel(r'$\ell(\ell + 1)|\Theta_{\ell}(k)|^2$/k', fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.xlabel("$k\eta_0$", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.legend()
plt.savefig('./Figs/M4/Integrand.pdf')

plt.show()



