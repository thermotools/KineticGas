import numpy as np
import matplotlib.pyplot as plt
from pykingas.LJSpline import LJSpline
from scipy.constants import Boltzmann, Avogadro
from matplotlib import colormaps
from matplotlib.colors import Normalize

cmap = colormaps['jet']
norm = Normalize(vmin=1,vmax=3,clip=False)
norm2 = Normalize(vmin=0.1,vmax=0.9,clip=False)

sig = 3.4e-10
epsilon_div_k = 117.0
MW = 40.0
m = MW*1e-3/Avogadro
eps = epsilon_div_k*Boltzmann
prefac = [(m*eps)**0.5/sig**2,
           Boltzmann/sig**2*(eps/m)**0.5,
           sig*(eps/m)**0.5]

ljs = LJSpline(sig,epsilon_div_k,MW)

def KristiansenLowDensity(T, transport_id): #id = 0 => viscosity, id = 1 => conductivity, id = 2 => selfdiffusion
    i = transport_id
    A = [[-2.2319,0.5703,0.3485,0.0811,-0.0879,0.9005],
         [-0.9044,0.5699,0.3470,0.0815,-0.0723,0.9017],
         [-1.9791,0.5419,0.3771,0.1146,-0.2484,1.0316]]
    return np.exp(A[i][0] + A[i][1]*np.log(T) + (A[i][2] + A[i][3]*np.log(T))*np.tanh( (np.log(T)-A[i][4]) / A[i][5] ))

temps = [1.0,2.0,5.0] #,10.0]
rhos = np.linspace(0.0001, 0.3, 10)

fig, ax = plt.subplots(1,3, figsize = (15,4))
i = 1
for T in temps:
    v = []
    l = []
    d = []
    for rho in rhos:
        print(i)
        i += 1
        vm = sig**3*Avogadro/rho
        t = T*epsilon_div_k
        v.append(ljs.viscosity(t,vm,[0.5,0.5]) /prefac[0])
        l.append(ljs.thermal_conductivity(t,vm,[0.5,0.5])/prefac[1])
        d.append(ljs.selfdiff_ljs(t,vm,[0.5,0.5])/prefac[2] * rho)
    v_e = KristiansenLowDensity(T, 0)
    l_e = KristiansenLowDensity(T,1)
    d_e = KristiansenLowDensity(T,2)
    ax[0].plot(rhos,v, color = cmap(norm(T)), label = "$T^*$ = " +str(T))
    ax[0].plot([0,0.3], [v_e,v_e], "--", color = cmap(norm(T)))
    ax[1].plot(rhos,l, color = cmap(norm(T)), label = "$T^*$ = " +str(T))
    ax[1].plot([0,0.3], [l_e,l_e], "--", color = cmap(norm(T)))
    ax[2].plot(rhos,d, color = cmap(norm(T)), label = "$T^*$ = " +str(T))
    ax[2].plot([0,0.3], [d_e,d_e], "--", color = cmap(norm(T)))

ax[0].set_ylabel(r"$\eta^*$", fontsize = 18)
ax[0].set_xlabel(r"$\rho^*$", fontsize = 18)
ax[0].tick_params(axis='both', direction = "in", length = 7, labelsize = 14)
ax[0].legend(fontsize = 12, loc = 2)

ax[1].set_ylabel(r"$\lambda^*$", fontsize = 18)
ax[1].set_xlabel(r"$\rho^*$", fontsize = 18)
ax[1].tick_params(axis='both', direction = "in", length = 7, labelsize = 14)
ax[1].legend(fontsize = 12)

ax[2].set_ylabel(r"$D \rho^*$", fontsize = 18)
ax[2].set_xlabel(r"$\rho^*$", fontsize = 18)
ax[2].tick_params(axis='both', direction = "in", length = 7, labelsize = 14)
ax[2].legend(fontsize = 12)

plt.tight_layout()
plt.show()
plt.close()