from pykingas.MieKinGas import MieKinGas
from thermopack.cubic import PengRobinson

eos = PengRobinson('H2,CO') # ThermoPack saftvrmie does not have parameters for CO (see: pyExamples/custom_eos.py)
mie = MieKinGas('H2,CO', use_eos=eos) # Modelling a hydrogen/carbon monoxide mixture (syngas)

T = 1300 # Kelvin (typical steam reforming temperature)
p = 100e5 # pascal (= 100 bar, approximate pressure of syngas stream)
x = [0.7, 0.3] # Molar composition

visc = mie.viscosity_tp(T, p, x, N=2) # Computing viscosity as a function of pressure and temperature

print(f'Shear Viscosity at T = {T} K, p = {p / 1e5:.1f} bar is : {visc * 1e6:.2f} micro Pa s')

# Internally, this just computes the molar volume using
vm, = mie.eos.specific_volume(T, p, x, mie.eos.VAPPH)
print(f'Molar density at T = {T} K, p = {p / 1e5:.1f} bar is : {vm} m3 / mol')
# And feeds the call to
visc = mie.viscosity(T, vm, x, N=2) # Computing viscosity as a function of temperature and volume

print(f'Shear Viscosity at T = {T} K, Vm = {vm:.1f} m3 / mol is : {visc * 1e6:.2f} micro Pa s')

try:
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib import colormaps
    import numpy as np
    import os
except ImportError:
    print('You need pandas and matplotlib if you want to see the figure examples.')
    exit(0)

# data = pd.read_excel(f'{os.path.dirname(__file__)}/data/co2_viscosity.xlsx')
# data = data[(data['p'] < 50) & (data['p'] > 0)]
# data.to_csv(f'{os.path.dirname(__file__)}/data/co2_viscosity.csv')
data = pd.read_csv(f'{os.path.dirname(__file__)}/data/co2_viscosity.csv')

mie = MieKinGas('CO2')

isotherms = sorted(set(data['T']))
norm = Normalize(min(isotherms), max(isotherms))
cmap = colormaps['cool']
for T in isotherms:
    isodata = data[data['T'] == T].sort_values('p')
    visc = np.array([mie.viscosity_tp(T, p * 1e6, [0.5, 0.5], N=1) for p in isodata['p']])

    plt.plot(isodata['p'], np.array(visc) * 1e6, color=cmap(norm(T)), label=T)
    plt.scatter(isodata['p'], isodata['eta'] * 1e3, color=cmap(norm(T)))

plt.legend(title=r'$T$ (K)', ncols=2)
plt.title(r'Shear viscosity of CO$_2$' + '\ncompared to experimental data')
plt.xlabel(r'$p$ (MPa)')
plt.ylabel(r'$\eta$ ($\mu$Pa s)')
plt.show()
