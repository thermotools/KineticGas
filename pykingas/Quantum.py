from .libpykingas import cpp_Quantum
from .py_KineticGas import py_KineticGas

class Quantum(py_KineticGas):

    def __init__(self, comps):
        super().__init__(comps)
        self.cpp_kingas = cpp_Quantum(comps)
        self.eos = None
    
    def cross_section(self, i, j, l, E):
        return self.cpp_kingas.cross_section(i, j, l, E)

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
    