import numpy as np
import matplotlib.pyplot as plt
from pykingas.HardSphere import HardSphere
from pykingas.MieKinGas import MieKinGas

# Initialize model for a fluid mixture Argon/Krypton
hs = HardSphere('AR,KR')

T, Vm = 300, 1e-3 # State point
x = [0.1, 0.9] # Mole fractions
N = 2 # Enskog approximation order (must be >= 2 for thermal conductivity)
cond = hs.thermal_conductivity(T, Vm, x, N=N) # Compute thermal conductivity
print('Thermal conductivity of Argon/Krypton mixture with mole fractions', x)
print('At T =', T, 'K', 'Vm =', Vm, 'm3 / mol')
print('is', cond)
print()


# Initialize model for pure Neon, specify a default value for the Enskog approximation order
hs = HardSphere('NE', N=2)

# State point
T, Vm = 300, 1e-3 # State point
# Mole fractions must still be supplied for two components because
# pure fluids are treated as binary mixtures of identical particles
# x = [0.5, 0.5] is good for numerical stability
x = [0.5, 0.5]

cond = hs.thermal_conductivity(T, Vm, x) # Compute thermal conductivity (uses N=2, as specified upon initialization)
print('Thermal conductivity of Argon (using HS potential)')
print('At T =', T, 'K', 'Vm =', Vm, 'm3 / mol')
print('is', cond)
print()

mie = MieKinGas('C1') # Methane model using Mie interaction potential
cond = mie.thermal_conductivity(T, Vm, x, N=N) # All other models than HardSphere are significantly slower ...
print('Thermal conductivity of Argon (using Mie potential)')
print('At T =', T, 'K', 'Vm =', Vm, 'm3 / mol')
print('is', cond)
print()

# ... but many computations at the same temperature are fast
Vm_list = np.linspace(1e-4, 1e-3)
cond_list = np.empty_like(Vm_list)
for i, Vm in enumerate(Vm_list):
    # Only needs to compute collision integrals at 300 K in the first iteration
    # Subsequent iterations are fast.
    cond_list[i] = mie.thermal_conductivity(T, Vm, x, N=N)

plt.plot((1 / Vm_list) * 1e3, cond_list)
plt.xlabel(r'$\rho$ [kmol m$^{-3}$]')
plt.ylabel(r'$\lambda$ [W m$^{-1}$ K$^{-1}$]')
plt.show()

T = 310 # New temperature
cond = mie.thermal_conductivity(T, Vm, x, N=N) # First computation at new temperature is slow
T = 300 # Resetting to old temperature
cond = mie.thermal_conductivity(T, Vm, x, N=N) # Still remembers tabulated collision integrals from before (is fast)

