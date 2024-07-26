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

"""
Note that as a default, parameters specified for MieKinGas will not be forwarded to the EoS, so we can use
different parameter sets for the EoS and RET models if we wish.
"""

comps = 'H2'
eos = saftvrmie(comps, parameter_reference='HYVA-FH0')
mie = MieKinGas(comps, use_eos=eos)

"""
This comes with the caveat that if we wish to specify parameters for MieKinGas, and want these parameters
to be used also for the EoS (for example if we are comparing to simulations, rather than trying to predict
the behaviour of real fluids) we need to specify this.
"""

comps = 'LJF'
mie = MieKinGas(comps, lr=[16, 16], use_default_eos_param=False) # exponent lr=16 will be forwarded to saftvrmie

"""
Writing your own EoS class is relatively simple, it must simply contain methods with signatures equivalent
to the below dummy class.
"""

import numpy as np

class DummyEOS:

    def __init__(self, comps):
        self.ncomps = len(comps.split(',')) # Not required, only used for this example
        self.VAPPH = None # Object must have this attribute, it can be whatever you like.

    def chemical_potential_tv(self, T, V, n, dmudn=None):
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

    def specific_volume(self, T, p, n, phase):
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
        return return_tuple # NOTE: ThermoPack 2.1 returns everything as tuples


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