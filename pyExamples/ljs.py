import numpy as np
import matplotlib.pyplot as plt
from pykingas.LJSpline import LJSpline


ljs = LJSpline(3.4e-10,120.0,40.0)

rdf = ljs.get_rdf((3.4e-10)**-(3),300,[0.5,0.5])

print(r"$\rho^* = 1$")
print(r"$T^* = 1$")
print("RDF at contact: ", rdf)

rho_array = np.array([1/(20e-3)])#np.linspace(0.1,1.4,8)*1e6
Vm_array = 1/rho_array
eta_list = []
lambda_list = []
d_list = []

for vm in Vm_array:
    eta_list.append(ljs.viscosity(300,vm,[0.5,0.5]))
    lambda_list.append(ljs.thermal_conductivity(300,vm,[0.5,0.5]))
    d_list.append(ljs.selfdiff_ljs(300,vm,[0.5,0.5]))

print(eta_list)
print(lambda_list)
print(d_list)
plt.plot(rho_array, eta_list, "-o")
plt.show()