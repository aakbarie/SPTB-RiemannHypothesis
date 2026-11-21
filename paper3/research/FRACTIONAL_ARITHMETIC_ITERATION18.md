# Paper III — Iteration 18: fractional arithmetic measures

**Status:** New Paper III section written, compiled locally in LaTeX, and numerically checked. The full revised manuscript is a 38-page PDF; the standalone technical addendum is 5 pages. The new results are written proofs, **not Lean verified**. The current draft research branch is not a proof of RH.

## Recovered historical constructions
The author's earlier *Structural Derivation of the Zeta Function and the Riemann Hypothesis* and *Appendix A: Structural Derivation of the Riemann Hypothesis via Harmonic Recursion* are preserved verbatim in the separately packaged research archive; neither claimed RH proof is adopted.

For σ>1 define μσ=Σ(n≥1)n^-σ δ_{log n}. For 0<d<1 and Re(s)>0:

- Lμσ(s)=ζ(σ+s).
- L(I^d μσ)(s)=s^-d ζ(σ+s).
- The normalized atomic measure has D1=0 and entropy log ζ(σ)-σ ζ'(σ)/ζ(σ).
- On any fixed bounded positive interval, the normalized fractional smoothing has D1=1, independently of d.

For νσ=Σ(n≥2)Λ(n)n^-σ δ_{log n}, its transforms are -ζ'/ζ at σ+s and -s^-d ζ'/ζ, respectively. Upon meromorphic continuation, zeta zeros give poles of residue **minus** their multiplicity, and multiplying by a nonvanishing fractional factor cannot move these poles.

## Counterexample to old entropy claim
With w(u)=1/2+sin(u)/4 and π_u=(1-w(u))δ0+w(u) Unif(0,e^-u), the Shannon entropy at every fixed quantizer resolution eventually becomes identically zero, but D1(π_u)=w(u) keeps oscillating.

Thus the fractional integration order d, information dimension D1, and ARFIMA long-memory parameter are distinct concepts.

## Empirical diagnostics
Exact integrated bin masses were used to avoid sampling singular fractional densities. At σ=1.6, T=3, fitted entropy slopes at the largest grids were approximately atomic 0.000, d=.2 0.909, d=.5 0.998, d=.8 1.000. Prime-power truncated transforms converge numerically to -ζ'/ζ. The analytic derivations do not rely on these tests.

## Research implication
Fractional integration is a useful positive arithmetic observation transform. The von Mangoldt measure links prime powers to zeta-zero poles. Neither entropy scaling nor smoothing independently supplies a positive invariant metric compatible with a genuine arithmetic spectral operator; that remains the substantive RH-related open problem.

**Verification:** Both LaTeX documents compiled with no unresolved references or overfull boxes. New section is not part of the Lean CI theorem declarations, and no kernel verification for these analytic claims is asserted.
