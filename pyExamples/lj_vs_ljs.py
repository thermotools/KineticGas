import numpy as np
import matplotlib.pyplot as plt
from pykingas.LJSpline import LJSpline
from pykingas.MieKinGas import MieKinGas
from matplotlib import colormaps
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# Potential parameters
sig = 3.42e-10
eps_div_k = 124.0
mw = 40.0


#Ititalizing KineticGas objects
ljs = LJSpline(N = 2)
lj = MieKinGas('LJF', mole_weights=[mw, mw], sigma=[sig, sig], eps_div_k=[eps_div_k, eps_div_k], la=[6, 6], lr=[12, 12])

units_ljs = ljs.get_reducing_units()
units_lj = lj.get_reducing_units()

# Some densities and temperatures in LJ units
T_stars = np.array([1.5,2.0,3.0])

# Plots
fig, ax = plt.subplots(1,3, figsize = (3*6.4, 4.8), layout = "tight")
norm = Normalize(vmin=np.min(T_stars),vmax=np.max(T_stars),clip=False)
cmap = colormaps['winter']
plt.rc('font', family='serif')

# Looping through isotherms:
Npoints = 10

count = 0
for T in T_stars:
    print("Isotherms completed:", str(count) + "/" + str(len(T_stars)))
    count += 1
    rho_stars = np.linspace(0.01,0.81,Npoints)
    cond_current = [[],[]] # Lists of thermal conductivities for current isotherm
    visc_current = [[],[]] # Lists of shear viscosities for current isotherm
    Drho_current = [[],[]] # Lists of selfdiffusivity*density for current isotherm
    for R in rho_stars:
        cond_current[0].append(lj.thermal_conductivity(T * units_lj.T,1/(R*units_lj.rho), [0.5,0.5]) / units_lj.tcond)
        cond_current[1].append(ljs.thermal_conductivity(T*units_ljs.T,1/(R*units_ljs.rho), [0.5,0.5]) / units_ljs.tcond)
        visc_current[0].append(lj.viscosity(T*units_lj.T,1/(R*units_lj.rho), [0.5,0.5]) / units_lj.visc)
        visc_current[1].append(ljs.viscosity(T*units_ljs.T,1/(R*units_ljs.rho), [0.5,0.5]) / units_ljs.visc)
        Drho_current[0].append(lj.selfdiffusion(T*units_lj.T,1/(R*units_lj.rho)) * R / units_lj.D)
        Drho_current[1].append(ljs.selfdiffusion(T*units_ljs.T,1/(R*units_ljs.rho)) * R / units_lj.D)
    ax[0].plot(rho_stars, cond_current[0], "--", color = cmap(norm(float(T))))
    ax[0].plot(rho_stars, cond_current[1], color = cmap(norm(float(T))))
    ax[1].plot(rho_stars, visc_current[0], "--", color = cmap(norm(float(T))))
    ax[1].plot(rho_stars, visc_current[1], color = cmap(norm(float(T))))
    ax[2].plot(rho_stars, Drho_current[0], "--", color = cmap(norm(float(T))))
    ax[2].plot(rho_stars, Drho_current[1], color = cmap(norm(float(T))))


plt.suptitle("Solid lines: Lennard-Jones/spline" + " "*15 + "Dashed lines: Lennard-Jones", fontsize = 22)

ax[0].set_xlabel(r"$\rho^*$", fontsize = 22)
ax[1].set_xlabel(r"$\rho^*$", fontsize = 22)
ax[2].set_xlabel(r"$\rho^*$", fontsize = 22)

ax[0].set_ylabel(r"$\lambda^*$", fontsize = 22)
ax[1].set_ylabel(r"$\eta^*$", fontsize = 22)
ax[2].set_ylabel(r"$D^* \rho^*$", fontsize = 22)

ax[0].tick_params(axis='both', direction = "in", length = 7, labelsize = 20)
ax[1].tick_params(axis='both', direction = "in", length = 7, labelsize = 20)
ax[2].tick_params(axis='both', direction = "in", length = 7, labelsize = 20)

cbar = fig.colorbar(ScalarMappable(norm = norm, cmap = cmap), ax = ax[2])
cbar.set_label(r"$T^*$", fontsize = 22, rotation = 0,labelpad = 20)
fig.figure.axes[3].tick_params(axis="y", labelsize=18) 

plt.show()
plt.close()