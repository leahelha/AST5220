import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy import constants as const


file = np.loadtxt("./Xe_test.txt")


Xe = file[:, 0]
ne = file[:, 1]
x = file[:, 2]

# plt.xlim(-12,0)
# plt.yscale("log")
# plt.plot(x, Xe)
# plt.axvline(x=x_decoupled, color='magenta', linestyle='-', label="Decoupling")  # Adding vertical line at x_decoupled
# plt.axvline(x=x_re, color='black', linestyle='dashdot', label="Recombination")  # Adding vertical line at x_re
# plt.axvline(x=x_Saha, color='purple', linestyle='dashdot', label="Saha")  # Adding vertical line at x_Saha
# # plt.axhline(y=, color='blue', linestyle='dashdot', label="Freeze out")  
# plt.xlabel("$x$", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
# plt.ylabel("$X_e$", fontdict={'fontsize': 14, 'fontname': 'Times New Roman'})
# plt.savefig("./Figs/M2/Xe_vs_x.pdf")