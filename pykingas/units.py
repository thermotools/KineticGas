"""
The dumb struct Units, holder for the reducing scale of common variables, for easy construction
from a given model.

Usage:

    model = MieKinGas('AR')
    unt = model.get_reducing_units()
    T_red = np.linspace(1, 3)
    T = T_red * unt.T # T in kelvin
    rho = 0.3 * unt.rho # density in mol / m3
    D = model.interdiffusion(T, 1 / rho, [1]) # D in m2 / s
    D_red = D / unt.D

    etc.
"""
from .libpykingas import cppUnits

class Units:
    # Dummy class used for type hinting
    def __init__(self, m_unit, L_unit, T_unit):
        self.T, # Temperature (K)
        self.E, # Energy (J)
        self.L, # Length (m)
        self.m, # Mass (kg)
        self.V, # Volume (m3)
        self.t, # Time (s)
        self.F, # Force (N)
        self.speed, # Speed (m / s)
        self.rho, # Density (mol / m3)
        self.D, # Diffusion (m2 / s)
        self.p, # Pressure (Pa)
        self.visc, # Shear viscosity (Pa s)
        self.kvisc, # Kinematic viscosity (m2 / s)
        self.tdiff, # Thermal diffusivity (m2 / s)
        self.tcond = [0 for _ in range(15)] # Thermal conductivity (W / K m)

Units = cppUnits