"""
    BEWARE!

    These tests are used to catch changes that cause a massive, unexpected performance hit. Therefore, it is set up to run various calculations
    and checks the calculation time against some values that I've hardcoded in at one point. I've given all the calculations times some margin,
    so they shouldn't fail unless they use significantly more time than what they did when I wrote this file. 

    With that said: If you run these tests on a potato, they will definitely fail.

    Generated on an M1 MacBook Pro using 10 cores. 
"""
from pykingas.MieKinGas import MieKinGas
from pykingas.HardSphere import HardSphere
from pykingas.multiparam import ModTangToennies, Patowski, PatowskiFH, FH_ModTangToennies
import time

def test_mie():
    kin = MieKinGas('H2', is_idealgas=False)

    t0 = time.perf_counter()
    for T in range(300, 500, 30):
        kin.viscosity(T, 1, [1], N=2)
    t1 = time.perf_counter()
    non_ideal_time = t1 - t0

    t0 = time.perf_counter()
    for T in range(300, 500, 30):
        kin.viscosity(T, 1, [1], N=2)
    t1 = time.perf_counter()
    cache_time = t1 - t0

    f = 50 # Roughly the expected speedup factor from caching
    print(cache_time * f, non_ideal_time)
    assert cache_time < non_ideal_time / f

    kin = MieKinGas('H2', is_idealgas=True)

    t0 = time.perf_counter()
    for T in range(300, 500, 30):
        kin.viscosity(T, 1, [1], N=2)
    t1 = time.perf_counter()

    ideal_time = t1 - t0

    # print(ideal_time, non_ideal_time)
    assert ideal_time < non_ideal_time
    assert ideal_time < 3.0
    assert non_ideal_time < 3.0


def test_multicomp():
    kin = MieKinGas('AR,KR,XE', la=[6.1, 6.1, 6.1]) # To prevent using correlations

    T = 700
    Vm = 1 / 20
    x = [0.1, 0.6, 0.3]

    t0 = time.perf_counter()
    for x2 in (0.3, 0.6, 0.8):
        x = [0.1, x2, 0.9 - x2]
        kin.viscosity(T, Vm, x, N=2)
        kin.interdiffusion(T, Vm, x, N=2)
        kin.thermal_conductivity(T, Vm, x, N=2)
    t1 = time.perf_counter()
    elapsed = t1 - t0
    assert elapsed < 3.5

def test_hardsphere():
    kin = HardSphere('KR', is_idealgas=False)

    t0 = time.perf_counter()
    for T in range(300, 500, 30):
        kin.viscosity(T, 1, [1], N=2)
    t1 = time.perf_counter()

    non_ideal_time = t1 - t0

    kin = HardSphere('KR', is_idealgas=True)

    t0 = time.perf_counter()
    for T in range(300, 500, 30):
        kin.viscosity(T, 1, [1], N=2)
    t1 = time.perf_counter()

    ideal_time = t1 - t0

    # HardSphere rdf and transfer lengths are so fast that we can't reliably compare the non-ideal to the ideal time
    assert ideal_time < 0.1
    assert non_ideal_time < 0.1

def test_patowski():
    kin = Patowski('H2', quantum_active=False)

    t0 = time.perf_counter()
    for T in range(300, 500, 60):
        kin.viscosity(T, 1, [1], N=2)
    t1 = time.perf_counter()

    elapsed = t1 - t0
    # print(elapsed)
    assert elapsed < 7.0

    t0 = time.perf_counter()
    for T in range(300, 500, 60):
        kin.viscosity(T, 1, [1], N=2)
    t1 = time.perf_counter()

    elapsed_with_cache = t1 - t0
    # print(elapsed_with_cache)
    assert elapsed_with_cache < elapsed

    t0 = time.perf_counter()
    for T in range(100, 600, 5):
        kin.second_virial(0, 0, T)
    t1 = time.perf_counter()
    elapsed = t1 - t0
    assert elapsed < 0.02

def test_modtangtoennies():
    kin = ModTangToennies('AR', quantum_active=False)

    t0 = time.perf_counter()
    for T in range(300, 500, 60):
        kin.viscosity(T, 1, [1], N=2)
    t1 = time.perf_counter()

    elapsed = t1 - t0
    # print(elapsed)
    assert elapsed < 5.0

    t0 = time.perf_counter()
    for T in range(300, 500, 60):
        kin.viscosity(T, 1, [1], N=2)
    t1 = time.perf_counter()

    elapsed_with_cache = t1 - t0
    # print(elapsed_with_cache)
    assert elapsed_with_cache < elapsed

    t0 = time.perf_counter()
    for T in range(100, 600, 1):
        kin.second_virial(0, 0, T)
    t1 = time.perf_counter()
    elapsed = t1 - t0
    # print(elapsed)
    assert elapsed < 0.05

def test_FH():
    kin = FH_ModTangToennies('AR', FH_order=1)

    t0 = time.perf_counter()
    kin.selfdiffusion(200, 1, N=2)
    t1 = time.perf_counter()

    elapsed = t1 - t0
    # print(elapsed)
    assert elapsed < 13.0

test_mie()