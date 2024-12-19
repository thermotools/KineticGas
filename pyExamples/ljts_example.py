import numpy as np
import matplotlib.pyplot as plt
from pykingas.LJSpline import LJSpline
from pykingas.LJTS import LJTS
from pykingas.MieKinGas import MieKinGas
from scipy.constants import Boltzmann, Avogadro

sig = 3.4e-10
epsilon_div_k = 117.0
MW = 40.0
m = MW*1e-3/Avogadro
eps = epsilon_div_k*Boltzmann
prefac = [(m*eps)**0.5/sig**2,
           Boltzmann/sig**2*(eps/m)**0.5,
           sig*(eps/m)**0.5]

ljs = LJSpline(sig, epsilon_div_k, MW, True)
ljts = LJTS(sig, epsilon_div_k, MW, True)
lj = MieKinGas("AR", [MW,MW], [sig,sig], [epsilon_div_k,epsilon_div_k], la=[6,6], lr=[12,12], is_idealgas= True, N = 2)

rho = 1.0
dummy_v = sig**3*Avogadro/rho
# rho_red =  1. / (dummy_v*sig**3*Avogadro)

for T in [5.0]:
    t = t = T*epsilon_div_k
    print(r"T =", T)
    #print("LJ", lj.(t,dummy_v, [0.5,0.5]) / prefac[0], "LJS: ", ljs.viscosity(t,dummy_v, [0.5,0.5]) / prefac[0], " | ", "LJTS: ", ljts.viscosity(t,dummy_v, [0.5,0.5]) / prefac[0], "\n")
    #print("LJTS cond: ", ljts.thermal_conductivity(t,dummy_v, [0.5,0.5]) / prefac[1], "\n")
    print("LJTS selfdiff: ", ljts.selfdiff_ljs(t,dummy_v, [0.5,0.5]) / prefac[2], "\n")

