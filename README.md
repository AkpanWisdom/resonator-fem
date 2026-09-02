# resonator-fem

Finite element modal analysis of prestressed beams. C++ with Eigen.

Computes natural frequencies of a beam under axial load, and how those
frequencies drift when temperature and humidity histories generate stress
in a restrained member.

![mode shapes](docs/figures/mode_shapes.png)

First three mode shapes, clamped-clamped.

![path dependence](docs/figures/path_dependence.png)

Three temperature histories reaching the same final temperature produce
different final stress and frequency.
