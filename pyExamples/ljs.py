import numpy as np
import matplotlib.pyplot as plt
from pykingas.LJSpline import LJSpline


ljs = LJSpline(1e-10,100.0,40.0)

rdf = ljs.get_rdf(1e30,100,[0.5,0.5])

print(r"$\rho^* = 1$")
print(r"$T^* = 1$")
print("RDF at contact: ", rdf)

rho_array = np.linspace(0.1,1.4,8)*1e6
Vm_array = 1/rho_array
mu_list = []

for vm in rho_array:
    mu_list.append(ljs.viscosity(100,vm,[0.5,0.5]))

plt.plot(rho_array, mu_list, "-o")
plt.show()