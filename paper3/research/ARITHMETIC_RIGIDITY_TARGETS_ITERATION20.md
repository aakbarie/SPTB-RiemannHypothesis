# Paper III — Iteration 20: arithmetic rigidity targets and exact prime-ramp transfer

**Status:** Mathematical research note and independent finite-cutoff check. The complete derivation, ranked research targets, figures and executable verification code are in the accompanying Iteration 20 research package. **Not a new proof of RH.**

## Exact prime-derived filter transfer
For von Mangoldt weights, define
- A(x)=sum_{n<=x} Lambda(n)(x/n-1) (ramp used in the 2026 Alpöge proof account);
- B(x)=sum_{n<=x} Lambda(n)(1-n/x) (Paper III iteration-19 first Riesz mean).

Their Mellin transforms on Re(s)>1 are F(s)/[s(s-1)] and F(s)/[s(s+1)], respectively, with F(s)=-zeta'(s)/zeta(s). Consequently, with t=log x,

  B(exp t)=A(exp t)-2 int_0^t exp(-(t-u)) A(exp u) du.

The two causal filters are the Green functions for (D^2-D) and (D^2+D) with the same arithmetic source, hence formally adjoint differential expressions.

## Critical-scale dissipation theorem
Set gamma_E=Euler's constant, kappa=3/2+gamma_E,

  a(t)=exp(-t/2)[A(exp t)-exp(t)(t-1-gamma_E)],
  b(t)=exp(-t/2)[B(exp t)-exp(t)/2],
  v(t)=int_0^t exp(-3(t-u)/2) a(u)du.

Then b+kappa exp(-3t/2)=a-2v and v'=a-(3/2)v, v(0)=0, so

  int_0^T a^2 = int_0^T (b+kappa exp(-3t/2))^2
                + 2v(T)^2 + 2 int_0^T v^2.

This is an **unconditional positive dissipation identity**, not a proof that the input energy is bounded. Thus sqrt(C_B(X)) <= ||a||_{L²(0,log X)} + kappa/sqrt(3).

On a zeta-zero mode rho, the transfer factor is (rho-1)/(rho+1), nonzero for 0<Re rho<1. Therefore it suppresses neither hypothetical off-critical zeros nor actual critical zeros. Generic one-sided window upper bounds do **not** transfer directly to L²: a(t)=-exp(t/4) is bounded above and has unbounded energy.

The pointwise polynomial identity underlying the passivity calculation has been added to the Lean pilot with an axiom audit. This is NOT a Lean formalization of the convolution, source equation, integration or RH-related claims.

## Target ranking and solution prospects
1. **Sparse one-sided arithmetic window theorem:** investigate/compress the newly independent-kernel-verified Alpöge proof mechanism. Its arithmetic supply uses iterated theta lifts, multiplicity/divisibility, Bost slopes and Galois windows; a Landau/Mellin argument then excludes off-line zeros. Independent Lean kernels accepted the mathematical RH statements as documented by Sutherland (MIT), 2026-11-21. Broad expository review continues.
2. **Suzuki--Weil signed kernel positivity:** identify a relative scale-extension inequality for the prime-defined screw kernel. It is RH-equivalent but unlike the automatic positive observation Gram genuinely detects zero location; local positivity alone cannot be extended by an ordinary spectral-gap argument because compact spectra accumulate at zero.
3. **Direct prime-pair critical Gram:** bound centered C_B(X)=O(log X) using exact Lambda correlations. Existing positivity and PNT estimates are insufficient (explicit positive-weight counterexample from iteration 19).
4. **Rank-trace/moment enhancement:** extend Alpöge--Furman's prior 2/3 line-zero result beyond the two-moment/bandwidth obstruction via new higher-order arithmetic information.

## Numerical validation
Actual von Mangoldt weights to X=100000, log-time Gauss--Legendre integration split at each prime power:
- A energy: 3.120796648654101;
- B energy: 1.2397516992391653;
- positive dissipation term: 2.1253416734232666;
- equality residual: -3.686038019951e-14;
- pointwise transfer residual max: 8.174e-13.
No extrapolation to arbitrary X or RH is justified.

Primary research links:
- https://math.mit.edu/~drew/grh/
- https://alpo.ge/remark.pdf
- https://nasqret.github.io/riemann-rev/atlas/
- https://doi.org/10.1112/jlms.12785
- https://arxiv.org/abs/2608.13637
