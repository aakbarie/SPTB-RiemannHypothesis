# Developing the full-field confinement criterion

The revised [proof draft](complementary_spline_observations.md) states a precisely delimited equivalence for the **fixed full harmonic field**, alpha >= 1, and the exact canonical mesh `h(T) = kappa / log(T)`:

1. Every zero has beta <= sigma.
2. Full-field integrated squared amplitude is O(T).
3. The blockwise affine amplitude residual E(T) is O(T).
4. E(T) is O(T^p) for some fixed finite p.

The bounds are for every sufficiently large real T. The forward linear estimate is conditional on confinement. At sigma = 1/2, the theorem is an RH equivalence criterion; it does not independently establish any of its equivalent boundedness conditions.

## Review record

Three separate AI agents developed/audited the mesh, spectral, and manuscript arguments. A fourth agent, without the preceding conversation, audited the assembled revision. This is AI-agent review, not independent human peer review or a formally verified proof.

- [Mesh reconstruction audit](reviews/mesh_audit.md): coercivity, jump summation, affine drift, and counterexamples delimiting the assumptions.
- [Spectral audit](reviews/spectral_audit.md): local L2 convergence at alpha = 1, primitive and pole argument, conditional O(T) estimate. Its optional full-energy limsup result is separate from the main criterion.
- [Manuscript integration audit](reviews/manuscript_audit.md): exact source defects and a proposed consistent derivative-penalty variant. That optional variant is not part of the main theorem.
- [Final referee report](reviews/final_referee.md): assessment of the assembled proof and its relation to the original conjecture.

The reports found no fatal gap in the precisely stated full-field, exact-mesh argument. The original manuscript's broader canonical regime is **not automatically covered**: comparability to 1/log(T) alone does not supply the complementary observations. A premise uniform over all canonical choices would include the exact choice, but a premise for one arbitrary comparable mesh is not addressed.

## Reproducible local checks

[check_complementary_mesh.py](check_complementary_mesh.py) uses only the Python standard library. Run from the repository root:

```sh
python3 docs/research/check_complementary_mesh.py
```

It verifies the local Gram matrix and determinant as exact rational polynomial identities, a hidden-step example, and finite mesh configurations. These computations supplement the analytic proof and do not certify the infinite-dimensional theorem.

## Paper development sequence

1. State the general complementary-mesh reconstruction theorem for locally square-integrable functions, including the common fixed function and exact mesh.
2. Establish the full zero field by the unit-bin convergence argument; remove the incorrect absolute-convergence assertion at alpha = 1.
3. Prove the conditional uniform unit-window bound and hence full-field O(T) energy under confinement.
4. Prove the primitive/pole implication and state the amplitude criterion.
5. Repair the definition of the derivative penalty before claiming a result for the existing SPTB functional. Full-amplitude projection and truncated-derivative projection must be distinguished.
6. Remove or repair the original false affine/derivative lower-bound lemmas and their dependent sharp detection claims. The new converse does not use them.
7. Keep moving-cutoff simulations separate from the fixed-full-field theorem; a changing amplitude cutoff requires a transfer estimate.
8. Seek independent human mathematical review and assess the related approximation-theory and analytic-number-theory literature before claiming publication readiness or novelty.

The main manuscript and compiled PDFs have not been revised in this research update. Source-location references in the reports describe the manuscript inspected during this review.
