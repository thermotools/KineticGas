import numpy as np
import matplotlib.pyplot as plt
from pykingas.LJSpline import LJSpline
from pykingas.LJTS import LJTS
from pykingas.MieKinGas import MieKinGas
from scipy.constants import Boltzmann, Avogadro
from scipy.optimize import curve_fit

sig = 3.4e-10
epsilon_div_k = 117.0
MW = 40.0
m = MW*1e-3/Avogadro
eps = epsilon_div_k*Boltzmann
prefac = [(m*eps)**0.5/sig**2,
           Boltzmann/sig**2*(eps/m)**0.5,
           sig*(eps/m)**0.5]

print(prefac)

#ljs = LJSpline(sig, epsilon_div_k, MW, True)
#ljts = LJTS(sig, epsilon_div_k, MW, True)
#lj = MieKinGas("AR", [MW,MW], [sig,sig], [epsilon_div_k,epsilon_div_k], la=[6,6], lr=[12,12], is_idealgas= True, N = 2)

# rho = 1.0
# vm = sig**3*Avogadro/rho
# T_vals = np.linspace(0.2,50.2, 25)
# x = [0.5,0.5]

# l = []
# v = []
# d = []
# i = 0
# for T in T_vals:
#     print(i)
#     i += 1
#     t = epsilon_div_k * T
    #l.append(ljts.thermal_conductivity(t,vm,x) / prefac[1])
    #v.append(ljts.viscosity(t, vm, x) / prefac[0])
    #d.append(ljts.selfdiff_ljs(t, vm, x) / prefac[2])
#d = [0.03168408370495478, 0.19241189498038946, 0.32877633380041843, 0.43948530812681924, 0.5362431997806039, 0.6244061652774723, 0.7058190829354826, 0.7820678240085702, 0.8541555852555965, 0.9228856809624]

# def func_v(T, a0, a1, a2, a3, a4, a5):
#     return np.exp(a0 + a1*np.log(T) + (a2 + a3*np.log(T))*np.tanh( (np.log(T)-a4) / a5 ))

# def func_l(T, a0, a1, a2, a3, a4, a5):
#     return np.exp(a0 + a1*np.log(T) + (a2 + a3*np.log(T))*np.tanh( (np.log(T)-a4) / a5 ))

# def func_d(T, a0, a1, a2, a3, a4, a5):
#     return np.exp(a0 + a1*np.log(T) + (a2 + a3*np.log(T))*np.tanh( (np.log(T)-a4) / a5 ))
    
# popt_d, pcov_d = curve_fit(func_d, T_vals, np.array(d))
# popt_l, pcov_l = curve_fit(func_l, T_vals, np.array(l))
# popt_v, pcov_v = curve_fit(func_v, T_vals, np.array(v))

# print("POPT V:", popt_v)
# print("POPT L:", popt_l)
# print("POPT D:", popt_d)

# fig, ax = plt.subplots(1,3, figsize = (14, 4))

# ax[0].plot(T_vals, v, "ro", label = r"Enskog values viscosity")
# ax[0].plot(T_vals, func_v(np.array(T_vals), *popt_v), "r-")
# ax[0].legend()

# ax[1].plot(T_vals, l, "bo", label = r"Enskog values thermal conductivity")
# ax[1].plot(T_vals, func_l(np.array(T_vals), *popt_l), "b-")
# ax[1].legend()

# ax[2].plot(T_vals, d, "go", label = r"Enskog values self-diffusion")
# ax[2].plot(T_vals, func_d(np.array(T_vals), *popt_d), "g-")
# ax[2].legend()

# ax[0].set_xlabel(r"$\rho^*$")
# ax[0].set_ylabel(r"$\lambda^*$")

# ax[1].set_xlabel(r"$\rho^*$")
# ax[1].set_ylabel(r"$\eta^*$")

# ax[2].set_xlabel(r"$\rho^*$")
# ax[2].set_ylabel(r"$D^* \rho^*$")

# fig.suptitle("Parametrization of Enskog values for LJTS")

# plt.show()
# plt.close()


