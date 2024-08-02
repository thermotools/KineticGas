import warnings

from pykingas.MieKinGas import MieKinGas
from thermopack.saftvrmie import saftvrmie
import pytest
from scipy.constants import Avogadro
from tools import check_eq

@pytest.mark.parametrize('T', [0.5, 0.8, 1.0, 1.3, 2.0, 3.0])
@pytest.mark.parametrize('rho', [0.01, 0.1, 0.3, 0.7])
@pytest.mark.parametrize('lr', [12, 15, 30])
def test_vs_thermopack(T, rho, lr):
    sigma = 3e-10
    eps_div_k = 100
    la = 6
    kgt = MieKinGas('LJF', sigma=[sigma, sigma], eps_div_k=[eps_div_k, eps_div_k], la=[la, la], lr=[lr, lr])
    tp = saftvrmie('LJF')
    tp.set_pure_fluid_param(1, 1, sigma, eps_div_k, la, lr)

    if not hasattr(tp, 'rdf_at_contact'):
        warnings.warn('This test requires a more recent version of ThermoPack!', FutureWarning)
        return

    print(T, rho)
    rho = rho / (Avogadro * sigma**3)
    particle_density = rho * Avogadro
    T = T * eps_div_k

    rdf_tp = tp.rdf_at_contact(T, 1, [rho], 2)[0][0]
    rdf_kgt = kgt.get_rdf(particle_density, T, [0.5, 0.5])[0][0]

    tc, _, _ = tp.critical([1])
    if T < tc:
        pd, _ = tp.bubble_pressure(T, [1])
        vl, = tp.specific_volume(T, pd, [1], tp.LIQPH)
        vg, = tp.specific_volume(T, pd, [1], tp.VAPPH)
        if (vl < 1 / rho) and (1 / rho < vg):
            return # Two-phase region

    assert check_eq(rdf_tp, rdf_kgt, 1e-3)