# Paper 3: metric rigidity formalization

Theorems about the real two-dimensional generator
A = [[eta, -gamma], [gamma, eta]] and positive metric
G = [[a,c],[c,b]]. The central claim is that
A-transformation cannot preserve G when eta != 0.
The quantitative inequality bounds |eta| by a normalized
diagonal symmetry defect.

These claims do not establish existence of any arithmetic metric
for zeta, and cannot be used as a proof of the Riemann Hypothesis.

Build from this directory:
lake exe cache get
lake build StructuralSpectralRigidity
lake env lean Check.lean

Until GitHub Actions succeeds and reports the axioms,
the source is a verification candidate, not a kernel-checked proof.
