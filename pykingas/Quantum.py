from .libpykingas import cpp_Quantum
from .py_KineticGas import py_KineticGas

class Quantum(py_KineticGas):

    def __init__(self, comps, is_idealgas=True):
        """Constructor
        Interface to quantum mechanical stuff.
        """
        super().__init__(comps, is_idealgas=True)
        self.cpp_kingas = None
        self.eos = None
    
    def potential(self, i, j, r):
        """Utility
        Potential
        """
        return self.cpp_kingas.potential(i, j, r)
    
    def potential_r(self, i, j, r):
        """Utility
        Potential derivative
        """
        return self.cpp_kingas.potential_r(i, j, r)
    
    def potential_rr(self, i, j, r):
        """Utility
        Potential second derivative
        """
        return self.cpp_kingas.potential_rr(i, j, r)

    def cross_section(self, i, j, l, E, reduced=False):
        """Utility
        Calculate the collision cross section. If `reduced=True`, return the cross section divided by the hard-sphere cross section.
        """
        return self.cpp_kingas.reduced_cross_section(i, j, l, E) if (reduced is True) else self.cpp_kingas.cross_section(i, j, l, E)

    def omega(self, i, j, n, s, T):
        r"""Utility
        Calculate the collision integral $\Omega^{(n, s)}$ as defined in The Limits of the Feynman-Hibbs corrections ... paper (see cite page).

        This method uses quantum mechanical or classical calculation based on whether `self.get_quantum_active()` is `True`
        """
        return self.cpp_kingas.omega(i, j, n, s, T)

    def quantum_omega(self, i, j, n, s, T):
        r"""Utility
        Calculate the quantal collision integral $\Omega^{(n, s)}$ as defined in The Limits of the Feynman-Hibbs corrections ... paper (see cite page).
        """
        return self.cpp_kingas.quantum_omega(i, j, n, s, T)

    def de_broglie_wavelength(self, i, T):
        """Utility
        Get the de Broglie wavelength of species `i` at temperature `T`.
        """
        return self.cpp_kingas.de_broglie_wavelength(i, T)

    def wave_function(self, i, j, l, E, r_end, dr=0.1):
        """Utility
        Solve the Schrödinger equation for the two-particle wave function at energy E, out to the distance `r_end`
        &&
        Args:
            i, j (int): Species indices
            l (int): Angular momentum quantum number
            E (float): Dimensionless energy (E / epsilon[i][j])
            r_end (float): Maximum particle separation (m)
            dr (float): Step size
        Returns:
            list[list[float]] : [r, psi, delta], the positions, wave function value, and local phase shift out to r_end
        """
        return self.cpp_kingas.wave_function(i, j, l, E, r_end, dr)
    
    def phase_shift(self, i, j, l, E):
        r"""Utility
        Compute the phase shift for a collision with angular momentum quantum number `l` and energy `E`
        Args:
            i, j (int): Species indices
            l (int): Angular momentum quantum number
            E (float): Dimensionless energy (E / epsilon[i][j])
        Returns:
            float: The relative phase shift $(- \pi / 2, \pi / 2)$
        """
        return self.cpp_kingas.phase_shift(i, j, l, E)
    
    def JKWB_phase_shift(self, i, j, l, E):
        r"""Utility
        Compute the phase shift for a collision with angular momentum quantum number `l` and energy `E`, using the JKWB approximation
        Args:
            i, j (int): Species indices
            l (int): Angular momentum quantum number
            E (float): Dimensionless energy (E / epsilon[i][j])
        Returns:
            float: The relative phase shift $(- \pi / 2, \pi / 2)$
        """
        return self.cpp_kingas.JKWB_phase_shift(i, j, l, E)
    
    def get_reducing_units(self, i=0, j=None):
        """Utility
        See `py_KineticGas`.
        """
        j = j if (j is not None) else i
        return self.cpp_kingas.get_reducing_units(i, j)
    
    def get_de_boer(self, i=None, j=None):
        """Utility
        Get the de Boer parameter
        """
        if (i is None) and (j is None) : return self.cpp_kingas.get_de_boer()
        elif (j is None): j = i
        return self.cpp_kingas.get_de_boer(i, j)
    
    def set_de_boer_mass(self, i, de_boer):
        """Utility
        Set the particle mass to get the specified de Boer parameter 
        """
        self.cpp_kingas.set_de_boer_mass(i, de_boer)

    def JKWB_upper_E_limit(self, i=0, j=None):
        """Utility
        Get the upper energy limit for when the JKWB approximation is automatically applied.
        """
        j = i if (j is None) else j
        return self.cpp_kingas.JKWB_upper_E_limit(i, j)
    
    def set_quantum_active(self, active):
        """Utility
        Activate/deactivate quantum mechanical calculation of things.
        """
        self.cpp_kingas.set_quantum_active(active)

    def get_quantum_active(self):
        """Utility
        Get the current quantum_active state.
        """
        return self.cpp_kingas.get_quantum_active()

class FH_Corrected:

    def __init__(self): pass
    
    def set_FH_order(self, FH_order):
        """Utility
        Set the Feynman-Hibbs order
        """
        self.cpp_kingas.set_FH_order(FH_order)
    
    def potential(self, r, T):
        return self.cpp_kingas.potential(0, 0, r, T)
    
    def potential_r(self, r, T):
        return self.cpp_kingas.potential_r(0, 0, r, T)
    
    def potential_rr(self, r, T):
        return self.cpp_kingas.potential_rr(0, 0, r, T)
    
    def potential_dn(self, r, T, n):
        """Utility
        Calculate the `n`'th derivative of the potential wrt. distance.
        """
        return self.cpp_kingas.potential_dn(0, 0, r, T, n)

    def get_r_min_eff(self, T):
        """Utility
        Get the position of the potential minimum at a specified temperature.
        """
        return self.cpp_kingas.get_r_min(0, 0, T)
    
    def get_sigma_eff(self, T):
        """Utility
        Get the position of the potential root at a specified temperature.
        """
        return self.cpp_kingas.get_sigma_eff(0, 0, T)
    
    def get_eps_eff(self, T):
        """Utility
        Get the potential well depth at a specified temperature.
        """
        return self.cpp_kingas.get_eps_eff(0, 0, T)
    
    def get_alpha_eff(self, T):
        """Utility
        Get the dimensionless van der Waals energy, (reduced using sigma_eff and eps_eff).  
        """
        return self.cpp_kingas.get_alpha_eff(0, 0, T)
