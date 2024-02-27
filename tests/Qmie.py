import matplotlib.pyplot as plt
import numpy as np
from plottools.cyclers import NormedCmap
from pykingas.QuantumMie import QuantumMie
from scipy.constants import Avogadro

qmie = QuantumMie('H2', FH_orders=1)
print(qmie.get_fh_order())
T_lst = [35, 60, 90]
p_lst = np.linspace(1e3, 1e7, 50)
# p_lst = np.linspace(1e3, 1.1e3, 2)

cmap = NormedCmap('cool', T_lst)
for T in T_lst:
    Vm_lst = np.array([qmie.eos.specific_volume(T, p, [1], qmie.eos.MINGIBBSPH)[0] for p in p_lst])
    rho_lst = Avogadro / Vm_lst

    mie = qmie.get_effective_mie_model(T)
    s = mie.sigma[0][0]
    plt.plot(rho_lst * s**3, [qmie.get_rdf(rho, T, [0.5, 0.5])[1][0] for rho in rho_lst], color=cmap(T))
    plt.plot(rho_lst * s**3, [mie.get_rdf(rho, T, [0.5, 0.5])[1][0] for rho in rho_lst], linestyle='', color=cmap(T), marker='|')

plt.show()