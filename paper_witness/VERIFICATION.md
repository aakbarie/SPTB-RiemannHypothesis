# Verification record

Run date: 21 November 2026.

## Mathematical checks

Command: `python3 verify_witness.py`.

| Check | Result |
| --- | --- |
| Finite observation-loss, expressibility, pair-bound, and decomposition cases | 1,215 exact rational cases passed |
| Nested-observation energy identities | 4,860 exact rational cases passed |
| Diagnostic-refinement universal-property instances | 538 qualifying instances passed |
| Distance-distributor inequalities | 24,025 finite cases passed |
| Directed-excess triangle inequalities | 29,791 finite cases passed |
| Cubic example, independently solved normal equations | Slopes 9/10 and 129/130; exact gap 3/325 |
| Complex exponential translation identity | Maximum relative discrepancy 4.936289043619882e-15 in the stated test grid |
| Dimensionless width/penalty rescaling | All six width/penalty combinations passed |
| False affine-bound scaling regression | c_star = 0.08343491543648048; linear-in-width scaling passed |
| Different Z/3 phases at equal scalar penalty | Counterexample reproduced |
| Invalid finite certificate inputs | Zero weights, mismatched lengths, and absent provenance rejected |

Exact enumeration used all 15 partitions of four states and all 81 targets with
values in {-1, 0, 1}, with unequal positive rational weights. Diagnostic augmentation
used all binary diagnostics. These are finite-instance checks, not formal proofs
of general statements.

Command: `python3 verify_witness.py --certificate example.json`.
The supplied example preserves the coarse observation, is a strict refinement,
and changes observation loss from 1 to 0. Each of the two obstruction pairs gives
the exact lower bound 1/2; projection gain is 1. The diagnostic's external
provenance remains visible in the output.

Dependencies used: NumPy 2.3.5, SciPy 1.17.0, SymPy 1.14.0. The finite certificate
mode and rational enumerations use the Python standard library.

## Manuscript build and inspection

- Built `main.tex` to a 17-page `main.pdf` with pdfTeX and latexmk.
- Final LaTeX log has no undefined references, undefined citations, overfull or
  underfull boxes, or LaTeX warnings.
- All 17 rendered pages were inspected for clipping, missing mathematics, table
  layout, page numbering, and reference readability.
- The execution image had missing generated TeX configuration, format, and font-map
  files. A local format and font map were generated from its installed TeX Live
  assets, and filesystem lookup was enabled for the build. No manuscript content
  was substituted to bypass the LaTeX build. A normal configured TeX Live
  installation uses the provided `make all` command.

## Limits

The manuscript contains written proofs checked during this revision. It has not
been independently refereed or formalized in a proof assistant. The numerical
tests do not establish infinite-spectrum cancellation bounds, the original
SPTB alpha=1 regime, a universal conservation law, or psychological/spiritual
equivalence. The compiler contract is a specification; the finite certificate
checker is the implemented example.
