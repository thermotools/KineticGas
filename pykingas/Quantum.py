from .libpykingas import cpp_Quantum
from .py_KineticGas import py_KineticGas

class Quantum(py_KineticGas):

    def __init__(self, comps, is_idealgas=True):
        super().__init__(comps, is_idealgas=True)
        self.cpp_kingas = None
        self.eos = None
    
    def potential(self, i, j, r):
        return self.cpp_kingas.potential(i, j, r)
    
    def potential_r(self, i, j, r):
        return self.cpp_kingas.potential_r(i, j, r)
    
    def potential_rr(self, i, j, r):
        return self.cpp_kingas.potential_rr(i, j, r)

    def cross_section(self, i, j, l, E, reduced=False):
        return self.cpp_kingas.reduced_cross_section(i, j, l, E) if (reduced is True) else self.cpp_kingas.cross_section(i, j, l, E)

    def quantum_omega(self, i, j, n, s, T):
        return self.cpp_kingas.quantum_omega(i, j, n, s, T)

    def de_broglie_wavelength(self, i, T):
        return self.cpp_kingas.de_broglie_wavelength(i, T)

    def wave_function(self, i, j, l, E, r_end, dr=0.1):
        return self.cpp_kingas.wave_function(i, j, l, E, r_end, dr)
    
    def phase_shift(self, i, j, l, E):
        return self.cpp_kingas.phase_shift(i, j, l, E)
    
    def JKWB_phase_shift(self, i, j, l, E):
        return self.cpp_kingas.JKWB_phase_shift(i, j, l, E)
    
    def get_reducing_units(self, i=0, j=None):
        j = j if (j is not None) else i
        return self.cpp_kingas.get_reducing_units(i, j)
    
    def get_de_boer(self, i=None, j=None):
        if (i is None) and (j is None) : return self.cpp_kingas.get_de_boer()
        elif (j is None): j = i
        return self.cpp_kingas.get_de_boer(i, j)
    
    def set_de_boer_mass(self, i, de_boer):
        self.cpp_kingas.set_de_boer_mass(i, de_boer)

    def JKWB_upper_E_limit(self, i=0, j=None):
        j = i if (j is None) else j
        return self.cpp_kingas.JKWB_upper_E_limit(i, j)
    
    def set_quantum_active(self, active):
        self.cpp_kingas.set_quantum_active(active)

    def get_quantum_active(self):
        return self.cpp_kingas.get_quantum_active()
    