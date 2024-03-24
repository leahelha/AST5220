import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy import constants as const


file = np.loadtxt("Xe_test.txt")


Xe = file[:, 0]
ne = file[:, 1]
x = file[:, 2]

plt.xlim(-12,0)
plt.yscale("log")
plt.plot(x, Xe)
plt.xlabel("x")
plt.ylabel("Xe")
plt.show()