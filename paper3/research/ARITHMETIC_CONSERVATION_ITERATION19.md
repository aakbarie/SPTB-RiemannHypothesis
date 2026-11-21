# Paper III / Iteration 19: Arithmetic conservation and prime-pair energy

Status: mathematical research draft, **NOT an RH proof**, and new analytic assertions **NOT Lean-verified**. Complete 46-page manuscript and 9-page standalone addendum compiled locally.

## Rigorous arithmetic source
For the locally finite positive prime-power measure ν₀=Σ_{n≥2}Λ(n)δ_{log n} and causal first-Riesz kernel k(t)=(1−e^{−t})1_{t≥0}, S=k*ν₀ satisfies (∂²+∂)S=ν₀. In Re(s)>1, L S(s)=−ζ'(s)/[s(s+1)ζ(s)].

Set r(t)=e^(−t/2)[S(t)−e^t/2]. From the classical Riesz explicit formula:
- RH iff r(t)=O(1), and iff sup_{T≥1} (1/T)∫₀ᵀ|r|²<∞.
- For B=sup_ρ Re ρ and η*=(B−1/2)_+, the **unconditional** exact limsup exponents are
  limsup_{t→∞} log(1+|r(t)|)/t = η* and
  limsup_{T→∞} log(1+∫₀ᵀ|r|²)/(2T) = η*.
These are detection equivalences, not an independent bound that forces η*=0.

## Explicit finite prime-pair correlation kernel
For w_n(x)=1_{x≥n}(1/x−n/x²), define K_X(m,n)=∫₁ˣw_m(x)w_n(x)dx and J_X(n)=∫₁ˣw_n(x)dx. Then
K_X(m,n)=[−1/x+(m+n)/(2x²)−mn/(3x³)]_{max(m,n)}^X,
J_X(n)=log(X/n)−1+n/X, and
C(X)=Σ_{m,n≤X}Λ(m)Λ(n)K_X(m,n)−Σ_{n≤X}Λ(n)J_X(n)+(X−1)/4.
Exact: C(X)=∫₁ˣ[S(log x)/x−1/2]²dx, so RH iff C(X)=O(log X). K_X is positive semidefinite unconditionally; its positivity alone does **not** prove the upper bound.

## Counterexample
For a>0, β in (1/2,1), perturbed weights Λ(n)(1+a n^(β−1)) retain nonnegativity, prime-power support, and Σ_{n≤x}Λ_modified(n)=x+o(x), but their meromorphic prime transform develops a pole at s=β. The critical-scale bound fails. The perturbation is intentionally not the same Euler product or functional equation as ζ.

## Independent checks and remaining goal
Mellin and xi-reflection checks agree to ≈10⁻⁶⁶ (65-digit arithmetic); a finite-X prime-pair Gram identity matched piecewise quadrature within 3×10⁻⁶⁰. Positive weights checked to cutoff X=1.5 million. Both PDF builds have no unresolved references or overfull boxes. None of this proves a bound for all X.

**Next open problem:** derive from genuinely arithmetic identities the cancellation estimate C(X)≤C log X for all large X. This is equivalent to RH and is NOT proven. Neither source balance nor zeta reflection alone supplies it. The new identities are not yet Lean-checked; do not label them formally verified.
