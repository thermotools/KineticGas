import numpy as np
import matplotlib.pyplot as plt
from pykingas.LJSpline import LJSpline
from pykingas.MieKinGas import MieKinGas
from pykingas.LJTS import LJTS
from matplotlib import colormaps
from matplotlib.colors import Normalize
from scipy.constants import Boltzmann as kB
from scipy.constants import Avogadro as Na 
import pandas as pd

cmap = colormaps['winter']
plt.rc('font', family='serif')

# Potential parameters
sig = 3.42e-10
eps_div_k = 124.0
eps = eps_div_k*kB
mw = 40.0
m = mw*1e-3/Na

# Conversion factor to get transport coefficients to LJ units
fac = [(m*eps)**0.5/sig**2,
         kB/sig**2*(eps/m)**0.5,
         sig*(eps/m)**0.5]

#Ititalizing KineticGas objects
ljs = LJSpline(sig,eps_div_k,mw)
lj = MieKinGas('LJF', mole_weights=[mw, mw], sigma=[sig, sig], eps_div_k=[eps_div_k, eps_div_k], la=[6, 6], lr=[12, 12])

# Some densities and temperatures in LJ units
T_stars = np.array([1.0,2.0,3.0, 5.0])

norm = Normalize(vmin=np.min(T_stars),vmax=np.max(T_stars),clip=False)

# Plots
fig, ax = plt.subplots(1,3, figsize = (3*6.4, 4.8), layout = "tight")

# Looping through isotherms:
Npoints = 10

count = 1
T_out = []
R_out = []
cond_current = [] # Lists of thermal conductivities for current isotherm
visc_current = [] # Lists of shear viscosities for current isotherm
Drho_current = [] # Lists of selfdiffusivity*density for current isotherm

for T in T_stars:
    print(count)
    count += 1
    rho_stars = np.array([0.01,0.11,0.21,0.31,0.41,0.51,0.61,0.71,0.81])
    for R in rho_stars:
        v_si = sig**3 * Na / R # Molar volume in SI units
        t_si = T * eps_div_k # Temperature in SI units
        cond_current.append(ljs.thermal_conductivity(t_si,v_si, [0.5,0.5]))
        visc_current.append(ljs.viscosity(t_si,v_si, [0.5,0.5]))
        Drho_current.append(ljs.selfdiffusion(t_si,v_si) / v_si)
        T_out.append(T)
        R_out.append(R)

df = pd.DataFrame({"T" : T_out, "rho" : R_out, "cond" : cond_current, "visc" : visc_current, "Drho" : Drho_current})
df.to_csv("transport_properties_ljs.csv")