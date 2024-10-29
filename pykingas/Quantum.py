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
    