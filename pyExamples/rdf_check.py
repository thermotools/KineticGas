import numpy as np
import matplotlib.pyplot as plt
from pykingas.LJSpline import LJSpline
from matplotlib import colormaps
from matplotlib.colors import Normalize

cmap = colormaps['jet']
norm = Normalize(vmin=1,vmax=3,clip=False)
norm2 = Normalize(vmin=0.1,vmax=0.9,clip=False)

ljs = LJSpline(1.0,1.0,40.0)

Ts = np.linspace(0.5,20,10)
g_lst = []
for T in Ts:
    gs = ljs.get_rdf(0.00000001, T, [0.5,0.5])
    print("T =", T, "RDF = ", gs)
    g_lst.append(gs[0])

plt.plot(Ts, g_lst)
plt.show()
plt.close()

temps = [1.0,1.25,1.5,2.0,2.5,3]
rhos = np.linspace(0.001, 0.9, 30)

for T in temps:
    g = []
    for rho in rhos:
        g.append(ljs.get_rdf(rho, T, [0.5,0.5])[0][0])
    plt.plot(rhos,g, color = cmap(norm(T)), label = r"$T^* =$ " + str(T))


plt.ylabel(r"$g^{LJ/s}(\sigma)$", fontsize = 18)
plt.xlabel(r"$\rho^*$", fontsize = 18)
plt.tick_params(axis='both', direction = "in", length = 7, labelsize = 14)
plt.legend(fontsize = 14)
plt.ylim([0.78,2.4])
plt.xlim([0,0.93])
plt.show()
plt.tight_layout()
plt.close()
