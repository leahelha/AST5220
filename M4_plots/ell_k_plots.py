import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy import constants as const


data_rec = np.loadtxt("./cells.txt")

l = data_rec[:, 0]
Cell = data_rec[:, 1] # Cell is normalized in text

plt.plot(l, Cell)
plt.title("Cell")
# plt.xlim(2, l[-1])
plt.yscale("log")
plt.xscale("log")
plt.xlabel("ell")
plt.show()


""" Make plot of Theta_100(k) from l=100"""