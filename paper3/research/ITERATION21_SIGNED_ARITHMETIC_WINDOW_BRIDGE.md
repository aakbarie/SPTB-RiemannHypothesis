# Paper III — Iteration 21: prime windows and signed Suzuki kernel

Status: revised 51-page manuscript and 5-page addendum compiled locally. This is research, not an independent RH proof. External Suzuki and Alpöge results are distinguished from our algebraic identities.

## Track A: independently checked arithmetic window
G(y) = sum_{n<=y} Lambda(n)(y/n - 1) - y(log y - 1 - gamma_E).
For each theta in (1/2,1), there is C_theta such that for every Y0 some Y>=Y0 satisfies G(y)<=C_theta y^theta for ALL y in [Y,Y^(3/2)]. The checked WindowBridge theorem implies RH. Conversely RH and the Riesz explicit formula imply G(y)=O(sqrt y). The arithmetic theta-lift, divisibility, Bost slopes, and cohomological estimate were not independently proved here.
Source: https://github.com/nasqret/riemann-rev/blob/main/references/archive-bridges/WindowBridge.lean
MIT verification: https://math.mit.edu/~drew/grh/

## Track B: exact Riesz-Suzuki intertwiner
Let R(t)=e^(-t/2) sum_{n<=e^t}Lambda(n)(1-n e^(-t)), and P(t)=sum_{n<=e^t}Lambda(n)/sqrt n (t-log n).
As causal distributions P''=(D+1/2)(D+3/2)R. Pointwise:
P(t)=R(t)+2 int_0^t R(u)du + (3/4) int_0^t (t-u)R(u)du.
Mellin multiplier quotient: (s+1/2)(s+3/2)/s^2.
Suzuki's exact prime-defined function is Psi(t)=Q(t)-P(t), with archimedean Q from equation (1.1) of Suzuki JLMS (2023). Suzuki Theorem 1.7 already proves RH iff Psi(t)>=0 for every t. No new RH proof is claimed.
Source: https://doi.org/10.1112/jlms.12785

## Indefinite prime-power refinement
For n=p^k, L=log n, w=Lambda(n)/sqrt n, the signed-kernel contribution is:
K_n(t,u)=w[(|t-u|-L)_+ - (|t|-L)_+ - (|u|-L)_+].
On [-a,a] only prime powers n<=exp(2a) enter; their activation threshold is a=L/2.
For L/2<s<L, on sample {s,-s}, the signed increment is [[0,c],[c,0]], with c=w(2s-L)>0 and eigenvalues +c,-c.
At n=2,s=.45, c=0.10138458037777. The full 2-point Suzuki kernel eigenvalues are ~0.03591411059424 and 0.12301928881685 (finite illustrative computation).
Psi has negative derivative jumps -Lambda(n)/sqrt n at log prime powers; it cannot be convex across them.

## Obstruction to purely local positivity extension
Smooth even bumps supported beyond the old aperture leave the entire old kernel exactly unchanged but can produce negative far diagonal entries. Suzuki's compact-window integral operators are trace class and lack a uniform coercive lower bound on infinite-dimensional L2. For a positive finite old block A, the signed extension requires C-B*inv(A)B>=0. Semidefinite A additionally needs a range condition and pseudoinverse. Arithmetic cross-block control, not monotone positive Gram refinement, is the missing target.

## Numerical and formal status
Actual Lambda(n), n<=10000, 60-digit precision: the prime-sum reconstruction residual is below 2.4e-58; finite transform residual below 2e-62. A smooth bump gives a negative far diagonal without changing the local kernel. Analytic identities have written proofs and finite tests; the Suzuki and Alpöge theorems are external. Four elementary theorems — prime_ramp_threshold, signed_prime_kick_indefinite, signed_schur_cleared, signed_schur_block_psd — passed Lean 4.33.1 compilation and the axiom audit (only propext, Classical.choice, Quot.sound) in https://github.com/aakbarie/SPTB-RiemannHypothesis/actions/runs/44918857512 . The full analytic transfer and Suzuki equivalences were not formalized.

## Open objective
Find an independent arithmetic source for the scalar inequality P(t)<=Q(t) for all t, or a relative scale-extension Schur inequality. The window theorem implies these statements through RH, but that logical deduction is not an independent structural conservation law.
