# resonator-fem

Finite element modal analysis of prestressed beams. C++ with Eigen.

Computes natural frequencies of a beam under axial load, and how those
frequencies drift when temperature and humidity histories generate stress
in a restrained member.

## Build

    cmake -S . -B build
    cmake --build build
    ctest --test-dir build

## Run

    ./build/modal_convergence
    ./build/prestress_sweep
    ./build/coupled_history

## License

MIT
