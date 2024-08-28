from pykingas.MieKinGas import MieKinGas
from thermopack.saftvrmie import saftvrmie
from thermopack.pcsaft import pcsaft
from thermopack.cpa import SRK_CPA

"""
KineticGas uses an equation of state internally to compute thermodynamic factors (chemical potential derivatives),
and for the TP-property methods, to compute the molar volume of a mixture at a given temperature and pressure

In some situations, we may want to use a different eos than the default, for example if we are modelling a mixture
that ThermoPack (which supplies the default EoS) does not have parameters for.

The default EoS for the `MieKinGas` class is the ThermoPack `saftvrmie` eos
"""

comps = 'H2O,C1' # Modelling a water/methane mixture
eos = SRK_CPA(comps)
mie = MieKinGas(comps, use_eos=eos) # Will now use CPA instead of SAFT-VR Mie

comps = 'CO2,H2O'
eos = pcsaft(comps) # CO2/Hydrogen mixture
mie = MieKinGas(comps, use_eos=eos) # Will use PC-SAFT instead of SAFT-VR Mie

comps = 'H2'
eos = saftvrmie(comps, parameter_reference='HYVA-FH0')
mie = MieKinGas(comps, use_eos=eos) # EoS now uses parameter set HYVA-FH0 from ThermoPack, rather than the default

"""
Note: If we specify parameters for MieKinGas, these are forwarded to the EoS, this is nice when we 
for example wish to model a "pure" Mie fluid
"""

comps = 'LJF'
mie = MieKinGas(comps, lr=[16, 16]) # exponent lr=16 will be forwarded to saftvrmie

"""
However, we may not always want this behaviour, so we can use the `use_default_eos_param` kwarg
"""

comps = 'H2'
eos = saftvrmie(comps, parameter_reference='HYVA-FH0')
mie1 = MieKinGas(comps, lr=[10, 10], sigma=[3e-10, 3e-10], use_eos=eos)
# mie1.eos uses the specified parameters, rather than the HYVA-FH0 values

mie2 = MieKinGas(comps, lr=[10, 10], sigma=[3e-10, 3e-10], use_eos=eos, use_default_eos_param=True) 
# mie2.eos uses parameter set HYVA-FH0 from ThermoPack, rather than the specified parameters

"""
Writing your own EoS class is relatively simple, it must simply contain methods with signatures equivalent
to the below dummy class.
"""

import numpy as np

class DummyEOS:

    def __init__(self, comps):
        self.ncomps = len(comps.split(',')) # Not required, only used for this example
        self.VAPPH = None # Object must have this attribute, it can be whatever you like.

    def chemical_potential_tv(self, T, V, n, dmudt=None, dmudv=None, dmudn=None):
        """
        The chemical potential of the mixture. Note: Only dmudn is actually used by the py_KineticGas class.
        Also note: We need to have V in the signature, even though it is not used, in order to be compatible with
        the signature of chemical_potential_tv in thermopack.
        
        Args:
            T (float) : Temperature [K]
            V (float) : Volume [m3]
            n (list[float]) : Mole numbers [mol]
            dmudn (bool) : Flag to activate derivative calculation
        Returns
            tuple : chemical potential [J / mol] and derivatives
        """
        print('DummyEOS is computing chemical potential')
        if dmudn is True:
            dmudn = np.identity(self.ncomps) + 3 # Compute dmudn_tv
            return 0, dmudn
        
        raise NotImplementedError('Only dmudn is implemented, because that is all we need for the pykingas package.')

    def specific_volume(self, T, p, n, phase, dvdt=None, dvdp=None, dvdn=None):
        """
        Compute molar volume. Note that we must have `n` and `phase` in the signature in
        order to be compatible with the signature of equation of state objects from ThermoPack.
        
        Args:
            T (float) : Temperature [K]
            p (float) : Pressure [Pa]
            n (list[float]) : mole numbers [mol]
            phase (int) : phase flag (see: ThermoPack)

        Returns:
            (float,) : molar volume [m3 / mol]
        """
        print('DummyEOS is computing specific volume!')
        vm = 5 # Compute specific volume
        return_tuple = (vm, )
        if dvdn is not None:
            dvdn = 2. # We only need to implement dvdn for KineticGas to be happy
            return_tuple += (dvdn, )
        return return_tuple # NOTE: ThermoPack 2.1 returns everything as tuples

    def pressure_tv(T, V, n):
        """
        Compute pressure for an ideal gas, we must take `x` as an argument to be compatible with the generic ThermoPack
        signature. This method is required for the solvent `frame_of_reference` diffusion and thermal diffusion coefficients.
        &&
        Args:
            T (float) : Temperature [K]
            Vm (float) : molar volume [m3 / mol]
            x (list[float]) : mole fractions [-]

        Returns:
            (tuple) : (Pressure,)
        """
        return_tuple = (n * T / V, ) # NOTE: ThermoPack 2.1 returns everything as tuples
        return return_tuple
    
    def idealenthalpysingle(T, ci, dhdt=None):
        h_id = 0. # We don't need this value for KineticGas, just to mimic the ThermoPack signature
        if dhdt is not None:
            dhdt = 24. # Compute ideal gas isobaric heat capacity
            return_tuple = (h_id, dhdt)
            return return_tuple

        return (h_id, )

comps = 'H2O'
eos = DummyEOS(comps)
mie = MieKinGas(comps, use_eos=eos)

T = 300
p = 1e6
Vm, = eos.specific_volume(T, p, None, eos.VAPPH) # Note the comma! (Because specific_volume returns a tuple)
print('Computing with custom eos at specified (T, V) ... ')
print(f'Diffusion coefficient with dummy eos : {mie.interdiffusion(T, Vm, [0.5, 0.5], N=1)}')
print(f'Shear viscosity with dummy eos : {mie.viscosity(T, Vm, [0.5, 0.5], N=2)}')
print(f'Thermal conductivity with dummy eos : {mie.thermal_conductivity(T, Vm, [0.5, 0.5], N=2)}')
print()
print('Computing with custom EoS at specified (T, P) ... ')
print(f'Diffusion coefficient with dummy eos : {mie.interdiffusion_tp(T, p, [0.5, 0.5], N=1)}')
print(f'Shear viscosity with dummy eos : {mie.viscosity_tp(T, p, [0.5, 0.5], N=2)}')
print(f'Thermal conductivity with dummy eos : {mie.thermal_conductivity_tp(T, p, [0.5, 0.5], N=2)}')