# C++ solver at the critical point

Faster C++ implementation of the critical-point solver (Algorithm 5 in the Appendix of the paper), used for the smallest fields, ln H = −50 and −54.

- `parametri.hpp`: parameters (J = J_c, grid, quadrature order, tolerances, values of ln H in `LM_LIST`).
- `gruzberg.hpp`: the solvers, at finite H and at H = 0 (reference solution).
- `main.cpp`: for each H = e^(−Lm) in `LM_LIST`, solves with Δt = 0.01 and writes `psi_Lm<Lm>.dat` (columns t, P_H(t), 1 − P_H(t), Δt); also writes the H = 0 solution to `psi_H0.dat`.

Compile and run:

```
g++ -std=c++17 -O3 main.cpp -o gruzberg
./gruzberg
```

The Bessel function comes from the C++17 special functions (GCC). With compilers that lack them (e.g. Apple clang), Boost is used instead: add `-I/path/to/boost`.
