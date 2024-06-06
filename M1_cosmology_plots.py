import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy import constants as const

""" COSMOLOGY PLOTS"""
"""
        chi2                h             OmegaM    OmegaK    Acceptrate    
Best params from min chi^2 [29.2799     0.701711   0.255027   0.0789514]


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
    fp << ddHpddx_of_x(x)        << " "; 
    fp << t_of_x(x)        << " "; 
    fp << get_luminosity_distance_of_x(x) << " ";
    fp << deta_of_x_dx(x) << " ";

    fp << get_TCMB(x) << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

"""
Gyr = 1/(60*60*24*365*1e9) # from s to Gyr
Mpc = 3.24*10**(-23) # from m to Mpc
Gpc = 3.24*10**(-25)
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
cosmo_deta = cosmo[:, 13]


""" Acceleration of the universe 
Acceleration is 0 when x = = -0.486721, z = 0.6269726206889359 and t = 7756.11 Myr
"""
cosmo_H = cosmo_Hp/np.exp(cosmo_x)
x_index = np.where((cosmo_dHpdx*cosmo_H)>0)[0][0] 

acc_z = 1/np.exp(cosmo_x[x_index-1]) -1
acc_x = cosmo_x[x_index-1] # This is the moment cosmo_dHpdx*cosmo_H is closest to 0
print(f'Acceleration is 0 when x = {acc_x} and z = {acc_z}') 




""" Calculating the times """
tol = 0.00005

rm_time = np.where( (abs( (cosmo_OmegaR+cosmo_OmegaNu) - (cosmo_OmegaB+cosmo_OmegaCDM)<tol )))[0]
md_time = np.where(abs((cosmo_OmegaB+cosmo_OmegaCDM) - (cosmo_OmegaLambda) )<tol)[0]


z_rm_time = 1/np.exp(cosmo_x[rm_time[0]]) - 1
z_md_time = 1/np.exp(cosmo_x[md_time[-1]]) - 1


print(f'Radiation-matter equality time  x = {cosmo_x[rm_time[0]]}, z = {z_rm_time}')
print(f'Matter-dark energy equality time  x = {cosmo_x[md_time[0]]}, z = {z_md_time}') #insert these in c code to find cosmic time
""" 
Radiation-matter equality time  = [-8.65778 -8.65772 -8.65766 ...  4.99988  4.99994  5.     ]
Matter-dark energy equality time  = [-0.255858]

c code:
Time of radiation-matter equality 0.0232026 Myr
Time of matter-dark energy equality 10379.9 Myr

rm_time: x = -8.65778, z = 5753.744937399699, t = 0.0232026 Myr
md_time: x = -0.255858, z = 0.2915693120751712, t = 10379.9 Myr 
"""




# ROUGHLY ILLUSTRATING THE DIFFERENT REGIMES IN THE PLOT
border_idx1 = rm_time#np.where((np.abs((cosmo_OmegaR+cosmo_OmegaNu)-(cosmo_OmegaB+cosmo_OmegaCDM)))<tol)[0] #Radiation-matter domination border
border_idx2 = md_time#np.where((np.abs((cosmo_OmegaB+cosmo_OmegaCDM)-cosmo_OmegaLambda))<tol)[0]  #Matter-dark energy domination border


idx1 = border_idx1[0]
idx2 = border_idx2[-1]

print(f'idx 1 = {idx1} and idx2 = {idx2}')

# plt.figure()
# # Patches for filling between lines
# region1 = patches.Rectangle((cosmo_x[0], 0), cosmo_x[idx1]-cosmo_x[0], 1, color='orange', alpha=0.2)
# region2 = patches.Rectangle((cosmo_x[idx1], 0), cosmo_x[idx2] - cosmo_x[idx1], 1, color='blue', alpha=0.2)
# region3 = patches.Rectangle((cosmo_x[idx2], 0), cosmo_x[-1] - cosmo_x[idx2], 1, color='magenta', alpha=0.2)

# # Adding the patches to the plot
# plt.gca().add_patch(region1)
# plt.gca().add_patch(region2)
# plt.gca().add_patch(region3)


# """ Plot of Omegas """

# plt.plot(cosmo_x, cosmo_OmegaR+cosmo_OmegaNu,  'orange', label=r"$\Omega_{R} = \Omega_{\gamma} + \Omega_{\nu}$")
# plt.plot(cosmo_x, cosmo_OmegaB+cosmo_OmegaCDM, 'blue', label=r"$\Omega_{M} = \Omega_{b} + \Omega_{CDM}$")
# plt.plot(cosmo_x, cosmo_OmegaLambda, 'purple', label="$\Omega_{\Lambda}$")
# plt.plot(cosmo_x, (cosmo_OmegaR+cosmo_OmegaNu+cosmo_OmegaB + cosmo_OmegaCDM + cosmo_OmegaLambda), color='black', linestyle='--', label='Sum')
# #plt.plot(cosmo_x, 1/cosmo_Hp*cosmo_dHpdx)

# # plt.title("$\Omega_i(x)$", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
# plt.xlabel("x", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
# plt.legend()
# plt.savefig("Figs/Omegas.pdf")



""" ANALYTICAL DERIVATIVES """



""" Hprime *its derivatives plots """


# plt.figure()
# ### 1/H * dHdx
# region1 = patches.Rectangle((cosmo_x[0], -1), cosmo_x[idx1]-cosmo_x[0], 2.5, color='orange', alpha=0.2)
# region2 = patches.Rectangle((cosmo_x[idx1], -1), cosmo_x[idx2] - cosmo_x[idx1], 2.5, color='blue', alpha=0.2)
# region3 = patches.Rectangle((cosmo_x[idx2], -1), cosmo_x[-1] - cosmo_x[idx2], 2.5, color='magenta', alpha=0.2)

# # Add the patches to the plot
# plt.gca().add_patch(region1)
# plt.gca().add_patch(region2)
# plt.gca().add_patch(region3)

# plt.plot(cosmo_x, 1/cosmo_Hp*cosmo_dHpdx, 'black', label=r'$\frac{1}{\mathcal{H}(x)} \frac{d\mathcal{H}(x)}{dx}$')
# plt.plot(cosmo_x, (1/cosmo_Hp)*cosmo_ddHpddx, color ='blue', label=r'$\frac{1}{\mathcal{H}(x)} \frac{d^2\mathcal{H}(x)}{dx^2}$')

# plt.hlines(y=-1, xmin=cosmo_x[0], xmax=cosmo_x[rm_time[0]], colors='orange', linestyles='--') # Radiation dominated Hp'/Hp
# plt.hlines(y=1, xmin=cosmo_x[0], xmax=cosmo_x[rm_time[0]], colors='orange', linestyles='--') # Radiation dominated Hp''/Hp

# plt.hlines(y=-0.5, xmin=cosmo_x[rm_time[0]], xmax=cosmo_x[md_time[0]], colors='blue', linestyles='--') # Matter dominated Hp'/Hp
# plt.hlines(y=0.25, xmin=cosmo_x[rm_time[0]], xmax=cosmo_x[md_time[0]], colors='blue', linestyles='--') # Matter dominated Hp''/Hp

# plt.hlines(y=1, xmin=cosmo_x[md_time[0]], xmax=cosmo_x[-1], colors='magenta', linestyles='--') # Dark energy dominated Hp'/Hp
# plt.hlines(y=1, xmin=cosmo_x[md_time[0]], xmax=cosmo_x[-1], colors='magenta', linestyles='--') # Dark energy dominated Hp''/Hp

# plt.legend()
# plt.xlabel("x", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
# plt.savefig("Figs/Hp_checks.pdf")



# """ Plot of Hprime(x)*eta(x)/c """
# plt.figure()

# region1 = patches.Rectangle((cosmo_x[0], 0), cosmo_x[idx1]-cosmo_x[0], 3.5, color='orange', alpha=0.2)
# region2 = patches.Rectangle((cosmo_x[idx1], 0), cosmo_x[idx2] - cosmo_x[idx1], 3.5, color='blue', alpha=0.2)
# region3 = patches.Rectangle((cosmo_x[idx2], 0), cosmo_x[-1] - cosmo_x[idx2], 3.5, color='magenta', alpha=0.2)

# # Add the patches to the plot
# plt.gca().add_patch(region1)
# plt.gca().add_patch(region2)
# plt.gca().add_patch(region3)

# eta_end = np.where(abs(cosmo_x)<0.001)[0][0]  # when x = 0, bc theres no point continuing eta checks passed that
# print(eta_end, 'ETA END')
# plt.plot(cosmo_x[:eta_end], cosmo_eta_of_x[:eta_end]*cosmo_Hp[:eta_end]/const.c, color ='blue', label=r'$\frac{\eta(x)\mathcal{H}(x)}{c}$')
# plt.plot(cosmo_x[:eta_end], cosmo_deta[:eta_end]*cosmo_Hp[:eta_end]/const.c, color='black', label=r'$\frac{d\eta(x)}{dx}\frac{\mathcal{H}(x)}{c}$')
# # plt.plot(cosmo_x, cosmo_eta_of_x*cosmo_Hp/const.c)
# plt.legend()
# plt.ylim(0, 3.5)
# plt.xlim(-18, 5)
# # plt.xlim(np.log(1/1+1089), 0)
# plt.xlabel("x", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
# plt.savefig("Figs/Hp_eta_checks.pdf")

# plt.show()

# """ Plot of cosmic time and conformal time """

# plt.figure()

# h = np.max(cosmo_t_of_x)*Gyr + 10*3

# region1 = patches.Rectangle((cosmo_x[0], 0), cosmo_x[idx1]-cosmo_x[0], h, color='orange', alpha=0.2)
# region2 = patches.Rectangle((cosmo_x[idx1], 0), cosmo_x[idx2] - cosmo_x[idx1], h, color='blue', alpha=0.2)
# region3 = patches.Rectangle((cosmo_x[idx2], 0), cosmo_x[-1] - cosmo_x[idx2], h, color='magenta', alpha=0.2)

# # Add the patches to the plot
# plt.gca().add_patch(region1)
# plt.gca().add_patch(region2)
# plt.gca().add_patch(region3)

# plt.plot(cosmo_x, cosmo_t_of_x*Gyr, color ='blue', label="t") 
# plt.plot(cosmo_x, cosmo_eta_of_x*Gyr/const.c, color ='orange', label="$\eta(x)/c$")
# plt.axhline(y=13.8, color='black', linestyle='--', label='13.8 Gyr', alpha=0.3)

# plt.yscale('log')
# #plt.xlim(np.log(1/1+1089), 0)
# plt.xlabel("x", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
# plt.ylabel("Gyr", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
# plt.legend()
# plt.savefig("Figs/cosmic_time_and_conformal_time.pdf")


""" Hubble factor Hprime(x) """


print(f"Hubble factor {cosmo_Hp[-1]*(1/(Mpc*100000))} fo x = {cosmo_x[-1]}")

plt.figure()

h = np.max(np.exp(-cosmo_x[:rm_time[0]]))#np.max(cosmo_Hp)*(1/(Mpc*100000)) 
region1 = patches.Rectangle((cosmo_x[0], 0), cosmo_x[idx1]-cosmo_x[0], h, color='orange', alpha=0.2)
region2 = patches.Rectangle((cosmo_x[idx1], 0), cosmo_x[idx2] - cosmo_x[idx1], h, color='blue', alpha=0.2)
region3 = patches.Rectangle((cosmo_x[idx2], 0), cosmo_x[-1] - cosmo_x[idx2], h, color='magenta', alpha=0.2)

# Add the patches to the plot
plt.gca().add_patch(region1)
plt.gca().add_patch(region2)
plt.gca().add_patch(region3)

plt.plot(cosmo_x, cosmo_Hp*(1/(Mpc*100000)), color ='blue',) #(100/(Mpc*1000)))
plt.plot(cosmo_x[:rm_time[0]], np.exp(-cosmo_x[:rm_time[0]]), color='blue', alpha = 0.5,  linestyle='--', label='1/a')
plt.axvline(acc_x, color='black', linestyle='--', label='Start of acceleration')
plt.yscale('log')
plt.legend()
#plt.xlim(-12, 0.1)
# plt.title("$\mathcal{H}(x)$", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.xlabel("x", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.ylabel("100 km/s / Mpc", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
plt.savefig("Figs/Hubble_factor.pdf")

plt.show()



