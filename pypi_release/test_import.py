from pykingas.MieKinGas import MieKinGas

print('\033[96mImported MieKinGas\033[0m')

mie = MieKinGas('H2,C1')
mie.viscosity(300, 0.02, [0.5, 0.5], N=1)

print('\033[92mComputed Something!\033[0m')