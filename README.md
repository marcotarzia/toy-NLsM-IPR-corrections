# Toy non-linear σ-model on the Bethe lattice: cavity-equation solvers

Codes accompanying

> T. Rizzo, M. Tarzia, and A. Vecchiarelli,
> *Anomalous corrections to the inverse participation ratio in the localized phase on the Bethe lattice: insights from a solvable "toy" non-linear σ-model* (2026),

and the companion Letter *Corrections to the IPR in the localized phase of the Anderson model on the Bethe lattice*.

The codes solve the cavity equation of the toy non-linear σ-model on the hyperbolic superplane (Zirnbauer; Gruzberg and Mirlin, 1996) on a Bethe lattice with branching ratio k = 2, in the presence of a symmetry-breaking field H. They compute the field distribution P_H(t) and the functions P_0(t), f_typ(t) and f_min(Δt) used in the scaling analysis, in the paramagnetic phase J < J_c and at the critical point J = J_c. The algorithms are described in the Appendix "Numerical solution of the cavity equations of the toy model" of the paper.

## Contents

| Path | Description |
|---|---|
| `mathematica/toy_NLsM.nb` | Mathematica notebook for the paramagnetic (localized) phase J < J_c (Algorithms 1–4 of the Appendix): P_H(t) for several H, P_0(t) and the constant c, the bulk zero mode f_typ(t), the boundary-layer function f_min(Δt), and the IPR proxy I_2(H). It also contains a first exploration of the critical case. |
| `mathematica/toy_NLsM_critical.nb` | Mathematica notebook for the critical point J = J_c (Algorithm 5): P_H(t) at J_c and the logarithmic corrections. |
| `cpp/` | Faster C++ implementation of the critical-point solver (Algorithm 5), used for the smallest fields (ln H = −50 and −54). |

## Parameters used in the paper

- k = 2.
- J = 0.01 (J < J_c), for which β = 0.7530159; J_c = 0.0265402, from k K_0(J_c) = K_{−1/2}(J_c).
- Grid on [t_min, t_max] with t_max = 10 and Δt ≈ 0.08. An end-corrected trapezoidal rule gives an error of order Δt⁴.

## How to cite

If you use these codes, please cite the paper above.

## License

MIT (see `LICENSE`).
