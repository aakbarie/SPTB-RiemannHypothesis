# Manuscript integration audit

Read-only audit of the current `aakbarie/SPTB-RiemannHypothesis` manuscript and local complementary-observation draft. No repository change made.

## Conditional reverse implication strengthens the proposed criterion

Assume every eta_rho=beta-sigma is nonpositive, and alpha>=1. Include signed frequencies from expanding cosines; let A_k be total absolute amplitude coefficient in [k,k+1). The local zero count gives A_k <= C log(2+|k|)/(1+|k|)^alpha, hence sum A_k^2 finite.

Fix chi in C_c^infty(R), equal to one on [0,1], supported in [-1,2]. Translate t=n+s, n>=0. The coefficient of each signed frequency becomes c_rho exp(eta_rho n) exp(i gamma_rho n), whose magnitude is at most |c_rho|. The remaining factor exp(eta_rho s) and its first two derivatives are uniformly bounded on support chi, since eta lies in a fixed bounded interval. Two integrations by parts therefore bound each cross kernel by C(1+|gamma-gamma'|)^(-2), with C independent of n and spectral truncation. Unit-bin convolution yields

    integral_n^{n+1} |H_sigma(t)|^2 dt <= C sum_k A_k^2.

Pass to the full field by local L2 convergence. Summing unit intervals proves integral_0^T |H_sigma|^2 <= C(T+1). Therefore E_H(T)<=C(T+1), regardless of the mesh.

Subject to the complementary-mesh converse's independent review, this proves equivalence of (i) beta<=sigma for every zero; (ii) E_H(T)=O(T); and (iii) E_H(T)=O(T^p) for some fixed p, for all sufficiently large real T with exact mesh h=kappa/logT. At sigma=1/2, functional-equation symmetry gives an RH criterion. This is an equivalence theorem, not an unconditional RH proof.

## A consistent derivative-penalized version

Define separately the full-field fit P_I H and cutoff-field fit P_I H_X, with X=1/h. Set

    F_sep(T)=E_H(T)+lambda(T) sum_I integral_I |(H_X-P_I H_X)'|^2,
    0<=lambda(T)<=C h(T)^2.

This definition deliberately differs from using H_X' minus the slope of P_I H. For F_sep the standard H1 stability of affine L2 projection applies to the same differentiable function H_X. On any interval of length l, its slope has the representation

    b = (6/l^3) integral_0^l s(l-s) g'(s) ds,
    l |b|^2 <= (6/5) integral_0^l |g'|^2.

Thus ||(g-Pg)'||_2^2 <= (22/5)||g'||_2^2 by the elementary squared-triangle bound. This constant is independent of l and works for the shortened last block.

Under eta<=0 the truncated derivative's bin coefficient masses satisfy B_k<=C log(2+|k|)(1+|k|)^(1-alpha). The same translated-window argument gives

    ||H_X'||_L2(0,T)^2 <= C(T+1) sum_{|k|<=X+1} B_k^2
                        <= C(T+1)(X+1)log^2(2+X).

For X=1/h tending to infinity and lambda<=C/X^2, F_sep(T)=O(T). Its amplitude term is unchanged, so the proposed converse applies. No minimum spacing or Hypothesis (S) is needed for this nonsharp upper bound.

## Exact manuscript locations requiring action

- `paper_clean/parts/part1.tex`, definition `def:harmonic-field`, lines108–134: alpha=1 does not ensure absolute convergence; square-summability alone does not establish L2 convergence for arbitrary clustered frequencies. Replace with unit-bin convergence proof and specify positive ordinates/multiplicity convention.
- Same file, `eq:SPTB`, lines15–25 and `rmk:truncation`, lines125–134: currently combines full amplitude/fits and cutoff derivative with ambiguous differentiation. Define E alone for main criterion, or the explicit separate-fit F_sep above.
- Same file, `eq:canon`, lines168–177: arbitrary comparable mesh must not be silently substituted in the converse. State exact h=kappa/logT and quantifier over all large real T.
- Same file, `lem:affine-lb`, lines213–296, and `lem:derivative-penalty`, lines298–318: false as stated. For eta=omega=1, f=1+t-t^3/3+O(t^4). Choosing 1+t yields amplitude residual <=Delta^7/63+O(Delta^8), contradicting claimed C Delta^3. Choosing derivative constant 1 gives <=Delta^5/5+O(Delta^6), contradicting c Delta. Retire claims and dependent sharp-detection results pending separate repair.
- Same file, `lem:high-freq`, lines371–384: displayed proof does not justify diagonal-only bound for clustered frequencies; replace if needed with unit-bin estimate.
- Same file, `thm:variance`, lines394–416: theorem statement omits beta<=sigma, and its unconditionality paragraph is incorrect. Conditional O(T) proof above can replace the entire argument for E or F_sep.
- Same file lines434–435: full field norm is bounded using a truncated coefficient sum without a tail argument. Replace with full unit-bin bound.
- `paper_clean/parts/part2.tex`, `thm:bias`, line199, and `prop:barrier`, lines260–271: overstate universal sharp detection and conflict with subsequent caveat at lines300–306. New converse supplies non-polynomial growth, not all-T sharp rate 2 eta.
- Same file `conj:horocycle`, lines293–309: retain as historical scope, state the new precisely delimited full-field criterion separately until reviewed.
- Same file `prop:attained-edge`, lines319ff: termwise negligibility does not justify infinite-sum negligibility; new pole proof avoids this.
- `paper_clean/appendices/appA.tex`, lines43–60: integration by parts for slope of full H cannot replace full derivative by cutoff derivative. F_sep repairs this by projecting H_X separately.
- Same appendix lines82–90: eta<=0 is essential; diagonal contributions from zeros exactly at beta=sigma cannot be replaced by contributions from all zeros as an asymptotic equality. Upper bounds suffice.
- Same appendix lines96–107: derivative-square summability threshold is alpha>3/2, not alpha>1; alpha=3/2 borderline omitted. Constants also depend on exact cutoff normalization.

## Integration sequence

1. General complementary-mesh reconstruction theorem for L2_loc functions.
2. Full zero-field construction, including alpha=1, using local zero counts.
3. Conditional uniform local L2 bound under beta<=sigma.
4. Primitive meromorphic transform and pole contradiction.
5. Exact amplitude criterion and optional F_sep corollary.
6. Explicitly separate moving-cutoff numerical diagnostic from full-field theorem, and retain numerical evidence as finite-window exploration.

No claim of novelty is established by this audit. No independent upper-bound premise for RH is proved.
