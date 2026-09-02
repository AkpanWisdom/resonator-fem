# resonator-fem

Finite element modal analysis of prestressed beams, coupled to an
environmental history model. Written in C++ with Eigen.

A clamped beam cannot change length, so temperature and humidity changes
produce axial stress rather than free expansion. That stress enters the
geometric stiffness and shifts the natural frequency. Because the material
creeps, the stress depends on the *history* of the environment, not only
its current state.

## Formulation

Two-node Euler-Bernoulli beam elements with Hermitian cubic shape functions.
Element matrices:

Free vibration gives the generalized symmetric eigenproblem

with `N > 0` in tension. Tension stiffens the beam and raises frequency;
compression softens it, and at the buckling load the first frequency
vanishes.

Axial stress follows from full restraint, `A₀σ + Σγᵢ + ε_env = 0`, where each
Kelvin unit integrates exactly over a step and the retardation times scale
with a microprestress state variable driven by the *rate* of environmental
change.

## Verification

Natural frequencies against the analytical clamped-clamped solution
`f_n = (βₙL)²/(2πL²)·√(EI/ρA)`, with `cos(βL)cosh(βL) = 1`:

| elements | f₁ error | f₂ error | f₃ error |
|---|---|---|---|
| 4  | 0.1327 % | 0.9250 % | 2.1357 % |
| 8  | 0.0084 % | 0.0632 % | 0.2370 % |
| 16 | 0.0005 % | 0.0040 % | 0.0153 % |
| 32 | 0.0000 % | 0.0003 % | 0.0010 % |
| 64 | 0.0000 % | 0.0000 % | 0.0001 % |

Errors fall by roughly 16× per refinement, the expected fourth-order rate
for cubic elements.

Further checks, all automated in `tests/`:

- static midspan deflection reproduces `PL³/192EI` to machine precision
- assembled mass sums to `ρAL`
- `f₁ → 0` as compression approaches `N_cr = 4π²EI/L²`
- high tension approaches the string limit `f = √(N/ρA)/2L`
- Kelvin relaxation matches the closed-form single-unit solution

## Path dependence

Three temperature histories reaching the same final temperature produce
different final frequencies:

| case | final σ (Pa) | drift (ppm) |
|---|---|---|
| fast ramp, then hold | −6.12e4 | −79,300 |
| slow ramp | −8.77e4 | −113,900 |
| up, down, up | −1.06e5 | −138,300 |

![path dependence](docs/figures/path_dependence.png)

The eigenproblem is re-solved only when the axial stress moves beyond a
tolerance, which reduces 2000 environmental steps to fewer than 200 solves.

## Build

Examples:

Eigen is fetched automatically by CMake.

## Limitations

Euler-Bernoulli kinematics, so shear deformation and rotary inertia are
neglected — adequate for slender members only. Microprestress and Kelvin
compliance parameters are illustrative rather than calibrated to a specific
material, so the drift magnitudes demonstrate the mechanism rather than
predict a device. The coupling is one-way: stress affects frequency, but
vibration does not feed back on the environmental state.

## License

MIT
