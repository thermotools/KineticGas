import numpy as np
import matplotlib.pyplot as plt
from pykingas.LJSpline import LJSpline
from pykingas.MieKinGas import MieKinGas
from pykingas.LJTS import LJTS
from matplotlib import colormaps
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy.constants import Boltzmann as kB
from scipy.constants import Avogadro as Na 

temps_sc = np.array([1.2, 1.3, 1.5, 2.0, 3.0])
rho_sc = np.linspace(0.001, 0.9, 30)

temps_sub = np.array([1.2, 1.3])
rho_sub =   np.linspace(0.7,0.9, 10)

temps_sc_2 = np.array([ 0.75, 0.8, 0.9, 1.5, 2.0])
rho_sc_2 = np.linspace(0.001, 0.9, 30)

temps_sub_2 = np.array([0.7, 0.8])
rho_sub_2 =   np.linspace(0.7,0.9, 10)


# Potential parameters
sig = 3.42e-10
eps_div_k = 124.0
eps = eps_div_k*kB
mw = 40.0
m = mw*1e-3/Na

# Plots
fig, ax = plt.subplots(1,2, figsize = (2*6.4, 4.8), layout = "tight")
norm = Normalize(vmin=np.min(temps_sc_2),vmax=np.max(temps_sc_2),clip=False)
cmap = colormaps['cividis']
plt.rc('font', family='serif')

#Ititalizing KineticGas objects
ljs = LJSpline(sig,eps_div_k,mw)
lj = MieKinGas('LJF', mole_weights=[mw, mw], sigma=[sig, sig], eps_div_k=[eps_div_k, eps_div_k], la=[6, 6], lr=[12, 12])

### LJ ###
# for t in (temps_sub):
#     rdf_c = []
#     t_si = t * eps_div_k
#     for r in rho_sub:
#         n = r / sig**3
#         rdf_c += [lj.get_rdf(n, t_si, [0.5, 0.5])[0][0]]
#     ax[0].plot(rho_sub, rdf_c, color = cmap(norm(t)))

for t in (temps_sc):
    rdf_c = []
    t_si = t * eps_div_k
    for r in rho_sc:
        n = r / sig**3
        rdf_c += [lj.get_rdf(n, t_si, [0.5, 0.5])[0][0]]
    ax[0].plot(rho_sc, rdf_c, color = cmap(norm(t)))

### LJS ###

# for t in (temps_sub_2):
#     rdf_c = []
#     t_si = t * eps_div_k
#     for r in rho_sub_2:
#         n = r / sig**3
#         rdf_c += [ljs.get_rdf(n, t_si, [0.5, 0.5])[0][0]]
#     ax[1].plot(rho_sub, rdf_c, color = cmap(norm(t)))

for t in (temps_sc_2):
    rdf_c = []
    t_si = t * eps_div_k
    for r in rho_sc_2:
        n = r / sig**3
        rdf_c += [ljs.get_rdf(n, t_si, [0.5, 0.5])[0][0]]
    ax[1].plot(rho_sc, rdf_c, color = cmap(norm(t)))
    
    

cbar = fig.colorbar(ScalarMappable(norm = norm, cmap = cmap), ax = ax[1])
cbar.set_label(r"$T^*$", fontsize = 22, rotation = 0,labelpad = 20)
fig.figure.axes[2].tick_params(axis="y", labelsize=18) 

ax[0].set_xlabel(r"$\rho^*$", fontsize = 22)
ax[1].set_xlabel(r"$\rho^*$", fontsize = 22)

ax[0].set_ylabel(r"$g(\sigma)$", fontsize = 22)
ax[1].set_ylabel(r"$g(\sigma)$", fontsize = 22)
ax[0].set_title("LJ", fontsize = 24)
ax[1].set_title("LJ-spline", fontsize = 24)

ax[0].tick_params(axis='both', direction = "in", length = 7, labelsize = 20)
ax[1].tick_params(axis='both', direction = "in", length = 7, labelsize = 20)

plt.show()
plt.close()