# Witness Functors revision record

Date: 21 November 2026.

## Purpose

Restore the original question about representational limits while repairing the
mathematical claims. A representation can hide a target distinction even when a
model fits everything the representation displays. The revised paper separates
those losses and gives explicit conditions for exposing a diagnostic through a
refinement. It treats the framework's own observation, model, loss, and transport
choices as possible objects of examination.

## Claim dispositions

| Earlier claim or omission | Revised treatment |
| --- | --- |
| Approximation cost alone describes representational adequacy | Separate observation loss from model error; prove the exact orthogonal decomposition in a finite weighted setting. |
| Refinement preserves information without specifying its meaning | Define fiber relations and factorization; preservation means recovery of prior observation values. |
| Higher-order witnessing follows from recursion | Prove the least diagnostic refinement property and show why deterministic postprocessing cannot recover collapsed distinctions. Diagnostic access remains an explicit obligation. |
| Squared loss is a Lawvere-enriched profunctor under addition | Give a counterexample to the additive triangle law; construct a distributor for unsquared distance using directed excess and state both actions. |
| Arbitrary system morphisms induce witness functoriality | Specify linear contractions preserving the signal and admissible models. |
| An infimum alone proves a right Kan-extension claim | Remove that assertion; supply the explicit universal property of diagnostic augmentation instead. |
| Information-preserving coordinate changes preserve cost | Require isometry or transport of the metric and model; show invertible scaling can change a fixed-metric cost. |
| A sheaf assumption bounds global error by the largest local error | Prove direct-sum additivity and the cost of compatibility constraints; retain the constant-fit counterexample. |
| WF equals the original SPTB functional | Prove the exact projection-gap identity; give the cubic example with gap 3/325. |
| A single-mode bound validates all-spectrum detection | Prove a fixed-width complex-mode identity, track dimensionless mesh scaling, and state the independent projected-remainder bound needed for sums. |
| SPTB's scale-uniform affine bound can be imported | Exhibit its rescaling counterexample and state the missing restrictions. |
| A new RH proof validates the separate SPTB estimates | Cite external verification as source context only; the damped-field result explicitly states zero confinement and alpha > 2. |
| A scalar phase penalty distinguishes every holonomy and is additive | Retain phase in a lifted path category; separate additive edge cost from phase mismatch in the appendix. |
| Receptive-field covariance proves WF as a strict supertheory | Replace with an explicit requirement to establish cost isometry for the chosen observation space; remove unsupported empirical validation claims. |
| Monotonicity applies to held-out prediction error | Restrict the theorem to exact nested-model optima under the same loss. |
| Physical, developmental, or mystical equivalence follows | Keep Kegan, Attar, and the tea reflection as motivation; state testable domain-specific obligations separately. |

## Scope of retained and new material

- Retain Hilbert projection, contraction monotonicity, dual certificates, and
  Hardy-space pole visibility from the reconstruction, with proofs.
- Add a finite observation theory, a diagnostic augmentation property, the
  impossibility of strict refinement by postprocessing alone, and a concrete
  example of making model metadata observable.
- Supply a genuinely specified enriched construction, rather than merely
  removing the categorical ambition.
- Move bridge history and phase to an appendix so the central argument remains
  about representations and their limitations.
- Omit the reconstruction's entropy--Fourier margin theorem from this manuscript.
  It is a separate result and is not needed for the representational argument;
  its earlier source has not been modified or repudiated here.
- Include a compiler contract and a finite certificate implementation. No
  arbitrary-paper parsing or automatic proof certification is claimed.

The projection identities, metric constructions, and factorization arguments are
standard mathematical ingredients. Their organization here is not evidence of a
new universal law. Novelty, usefulness, and cross-domain applicability require
further literature review and independent evaluation.

## Source basis and limits

SPTB inspected at main commit `14ba2b963ad4e2a81cdf8288507aae9da08e1f98`.
The original witness PDF, reconstructed TeX, and bridge and SPTB audit documents
were read before revision. The current main branch contains no witness manuscript;
`paper_witness/` therefore adds its editable source and compiled result. The older
SPTB source is unchanged, and the repository README links this revision and states
that its analytic obligations remain.

This is a substantive research revision. It has written proofs, finite checks,
symbolic checks, and numerical sanity checks; it has not received an independent
referee report or a formal proof-assistant verification.
