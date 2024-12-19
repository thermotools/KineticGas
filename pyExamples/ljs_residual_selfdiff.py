import numpy as np
import matplotlib.pyplot as plt
from pykingas.LJSpline import LJSpline
from scipy.constants import Boltzmann, Avogadro

sig = 3.42e-10
epsilon_div_k = 124.0

MW = 40.0
m = MW*1e-3/Avogadro
eps = epsilon_div_k*Boltzmann
prefac = [(m*eps)**0.5/sig**2,
           Boltzmann/sig**2*(eps/m)**0.5,
           sig*(eps/m)**0.5]

ljs = LJSpline(3.42e-10,124.0,40.0)
ljs_ig = LJSpline(3.42e-10,124.0,40.0, is_ideal=True)

rhos = np.linspace(0.2,1,100)
Ts = [1,2,3,5]

for t in Ts:
    d_res = []
    for r in rhos:
        vm = sig**3*Avogadro/r
        T = epsilon_div_k*t
        d_res.append((ljs.selfdiff_ljs(T,vm, [0.5,0.5]) - ljs_ig.selfdiff_ljs(T,vm, [0.5,0.5])) / prefac[2] * r)
    plt.plot(rhos,d_res, label = r"$T^* =$" + str(t))

plt.legend()
plt.show()
plt.close()
