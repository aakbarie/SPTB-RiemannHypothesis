# Review of the standalone manuscript

Reviewed manuscript: `core.md`, 10 September 2026.

This is an AI mathematical review of the analytic manuscript. It is not human peer review, a novelty assessment, or formal verification. The numerical section was not yet appended when this review was performed.

## Verdict

No fatal mathematical gap was found in the stated criterion for the fixed full zero field, alpha at least one, the exact origin-anchored mesh, and residual bounds at every sufficiently large real horizon. The argument proves an equivalence criterion, not an unconditional bound or the Riemann hypothesis.

## Corrections and clarifications

1. **Abstract: residues are nonzero, not positive.** The expression in equation (27) is a positive real coefficient sum divided by `2 z`, and is generally complex. Replace “positive residues” with “nonzero residues” or “noncancelling residues arising from positive coefficients.” The proof itself uses the correct nonzero-residue fact.
2. **Nullspace scope after equation (14).** The two-mesh argument shows that zero residual on both meshes forces a single affine function on the interior reconstruction interval `J_T`. Qualify “on the interval” explicitly as `J_T`; a statement about common knots elsewhere in the entire observation domain is unnecessary.
3. **Generic reconstruction hypothesis.** The paragraph starting “Suppose f ... satisfies (3)” should state `E_f(T)=O(T^p)` explicitly. Equation (3) is written for the zeta field, whereas Theorem 2 concerns arbitrary locally square-integrable functions.
4. **Standalone wording.** “Summing over both conjugates, as in the manuscript” and “the manuscript's broader admissible mesh regime” should identify the earlier SPTB manuscript [4], to distinguish it from the present manuscript.
5. **References pending.** Reference [1, Corollary 1] must actually supply the Riemann-von Mangoldt estimate with an `O(log T)` remainder, which implies (20). The numbered bibliography must be reconciled with the final document; unused references need no forced citations. This review cannot certify placeholders not yet populated.

## Mathematical checks

- The step-hinge matrix has positive first principal minor and determinant `q^4(1-q)^4/12`. On the indicated compact range its trace bound gives the stated explicit coercivity constant. Scaling correctly distinguishes the value-jump factor `h'` from the slope-jump factor `(h')^3`.
- Equation (9) follows exactly from the two logarithmic mesh widths. The range `[1/8,3/8]` follows from the displayed logarithm bounds. A crossing block cannot contain a second first-mesh knot, and distinct selected crossing blocks are disjoint.
- The jump traces belong to the finite spline projection, so no unsupported pointwise traces of an arbitrary L2 function are used. The positive-length convention for the first intersecting piece resolves a knot at the left endpoint of `J_T`.
- Summing jump increments yields the fourth-power reconstruction loss in (14). The subsequent logarithmic factor is absorbed by one further power of the horizon. The initial affine component is correctly retained in the implied constant.
- The overlapping geometric intervals have the asserted overlap. The rescaled affine Gram estimate controls both coefficient increments. Since `r=p+5>3`, both geometric sums converge in the required direction and yield (18), after adjoining an initial compact interval.
- The local series construction uses square-summable bin masses, not merely square-summable individual coefficients. The zero count controls multiplicity and permits arbitrarily close ordinates. Integration by parts against a fixed compact cutoff supplies a summable bin interaction kernel. This establishes the required local L2 convergence at alpha equal to one.
- The primitive has an absolutely and locally uniformly convergent series. Its identification with the integral of the L2 limit follows on every compact interval by Cauchy-Schwarz.
- Termwise Laplace integration is valid in a sufficiently far-right half-plane. The meromorphic expression converges normally away from its locally finite pole set. At any individual right-half-plane pole, multiplicities add positive coefficients to a nonzero residue; no attained spectral edge is required.
- Polynomial growth of the primitive gives a holomorphic transform on the whole open right half-plane. Analytic continuation on that half-plane with the locally finite poles removed then contradicts any nonremovable pole there.
- Under confinement, translation multiplies the coefficients by factors of magnitude at most one. The same fixed-window estimate is uniform in the translated integer, proving the linear full energy bound and completing all four implications.
- At sigma equal to one half, the functional equation and conjugation symmetry convert the one-sided confinement statement to RH. This does not establish the residual bound independently.
- The integer-knot piecewise-affine example correctly excludes arbitrary sparse observations and arbitrary asymptotically comparable mesh choices from the general reconstruction theorem. It does not refute a zeta-specific extension, and the manuscript states that distinction.

## Limits that must remain visible

- The analytic observations use one fixed full field. A numerical cutoff growing with the horizon does not constitute a proof of the full-field asymptotic hypothesis.
- The exact canonical mesh is part of the theorem. A bound uniform over a wider admissible class can imply the hypothesis if that class contains the exact mesh; a bound for one arbitrary comparable mesh need not do so.
- Finite numerical zero lists and synthetic growing modes can test implementation and finite-field behavior. They do not verify RH or exclude off-line zeros outside a cutoff.
- No result here supplies the older sharp exponential residual lower bound at every horizon. The pole contradiction proves impossibility of a uniform polynomial residual bound in the presence of a right-half-plane pole.

The analytic manuscript is suitable for assembly after the listed wording corrections and bibliography reconciliation. Independent human mathematical review remains appropriate before treating the result as established scholarship.

## Resolution record

All five presentation findings were addressed in the assembled `manuscript.md`: residues are called nonzero; the affine ambiguity is localized to J_T; the generic residual hypothesis is explicit; the earlier manuscript is identified as reference [4]; and reference [1] supplies the zero-counting estimate. Trudgian's Corollary 1 was checked directly in arXiv:1208.5846v2. The numerical section received the separate review recorded in `review_numerics.md`.
