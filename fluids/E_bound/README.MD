# Bound state energies

Some substances (like Argon) have a very rich roto-vibrational spectrum. In order to avoid cluttering the fluid files with massive tables of bound-state energies, 
the parameter `Quantum::has_E_bound_file` can be set to `true` in the fluid file. In that case, bound state energies will be retrieved from the `fluid.dat` file 
in this directory, instead of the `Quantum::E_bound_div_k` element of the fluid file.

## Units and organisation of energy files

All bound state energies are supplied in Kelvin. The conversion to Joules is

$$
E (\mathrm{J}) = E (\mathrm{K}) \cdot \mathrm{Boltzmann}
$$

and the conversion to $\mathrm{cm}^{-1}$ is

$$
E (\mathrm{cm}^{-1})  = E (\mathrm{K}) \cdot 10^{-2} \cdot \mathrm{Boltzmann} / (\mathrm{Planck} \cdot \mathrm{Speed}\ \mathrm{of}\ \mathrm{light}).
$$

The files are organised such that each line indicates a vibrational state, with the first line being the ground state. 
Each column indicates a rotational state, with the first column being the ground state.

The file

```
DummyFluid.dat
-1.5,-1.2,-0.7,-0.1
-0.9,-0.4,
-0.3
```

Thus indicates that `DummyFluid` has three vibrational states, with the ground state (both rotational and vibrational) being at -1.5 K. The zeroth vibrational state has three rotational states, (-1.5, -1.2, -0.7, and -0.1). The highest vibrational state has only one rotational mode: The state (-0.3) corresponds to the second excited vibrational state, in the rotational ground state.

The bound-state energy files are parsed by the method `Quantum::get_E_bound_from_file`, in `Quantum.[h/cpp]`