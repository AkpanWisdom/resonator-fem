# Beam element formulation

2-node Euler-Bernoulli, length L, EI, rho*A.
DOFs {v1, th1, v2, th2}.

## Shape functions

v = a0 + a1 x + a2 x^2 + a3 x^3, four nodal conditions.
xi = x/L.

N1 = 1 - 3 xi^2 + 2 xi^3
N2 = L (xi - 2 xi^2 + xi^3)
N3 = 3 xi^2 - 2 xi^3
N4 = L (-xi^2 + xi^3)

## Matrices

K_ij   = int EI  N_i'' N_j'' dx
M_ij   = int rhoA N_i   N_j  dx
Kg_ij  = int N    N_i'  N_j' dx

K = (EI/L^3) *
[  12    6L   -12    6L  ]
[  6L  4L^2   -6L  2L^2  ]
[ -12   -6L    12   -6L  ]
[  6L  2L^2   -6L  4L^2  ]

M = (rhoA*L/420) *
[ 156    22L    54   -13L  ]
[ 22L   4L^2   13L  -3L^2  ]
[  54    13L   156   -22L  ]
[-13L  -3L^2  -22L   4L^2  ]

Kg = (N/30L) *
[  36    3L   -36    3L  ]
[  3L  4L^2   -3L   -L^2 ]
[ -36   -3L    36   -3L  ]
[  3L   -L^2  -3L   4L^2 ]

N > 0 tension. K_tot = K + Kg.

## Eigenproblem

K_tot phi = omega^2 M phi

## Targets

Clamped-clamped: f_n = (bn L)^2 / (2 pi L^2) sqrt(EI/rhoA)
bn L = 4.730041, 7.853205, 10.995608   [cos(bL)cosh(bL) = 1]

N_cr = 4 pi^2 EI / L^2

Consistent mass, not lumped.