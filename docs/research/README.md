# Building out the Horocycle argument

This folder contains a proposed proof for the **fixed full harmonic field** with the exact canonical mesh `Delta(T) = kappa / log(T)`. It is a research draft awaiting independent mathematical review, not an established theorem or an RH proof.

## Files

- [complementary_spline_observations.md](complementary_spline_observations.md): full argument, assumptions, source alignment, and review targets.
- [check_complementary_mesh.py](check_complementary_mesh.py): exact rational polynomial checks of the local projection matrix and determinant, plus finite mesh-placement checks.

Run the checks from the repository root:

```sh
python3 docs/research/check_complementary_mesh.py
```

The script uses only the Python standard library. Its finite checks do not certify the complete analytic argument.

## Proposed paper development

1. **Fix the statement and quantifiers in Part 1.** Specify the same full field at every observation horizon, its local L2 interpretation, the exact canonical mesh, and the bound for every sufficiently large real T. Distinguish truncation of the derivative penalty from truncation of the amplitude itself.
2. **Review the complementary-mesh estimate.** Check the local step/hinge projection matrix, disjoint crossing blocks, and polynomial reconstruction of affine coefficients across overlapping intervals.
3. **Review the infinite-field argument.** Check unit-bin zero-counting control, local L2 convergence at alpha = 1, convergence of the primitive, and the isolated-pole contradiction.
4. **Integrate only after review.** Candidate material for Part 2: the complementary-observation lemma, polynomial-growth transfer, and full-field converse. Move detailed projection algebra and convergence arguments into appendices.
5. **Reconcile existing statements.** Reassess the current finite-configuration and attained-edge claims against the reviewed result. A qualitative converse does not establish the manuscript's sharper exponential lower bounds.
6. **Keep the RH premise explicit.** The argument assumes polynomial boundedness; it does not establish that bound independently. The variance-regime assumption cannot supply an unconditional RH proof.
7. **Handle numerical cutoffs separately.** Existing finite-zero simulations do not test the full-field theorem. A moving amplitude cutoff needs an additional transfer estimate.

The manuscript and compiled PDFs have not been edited as part of adding these research materials. No novelty claim is made.
