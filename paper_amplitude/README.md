# Amplitude confinement paper

[Read the paper (PDF)](manuscript.pdf) or [editable manuscript](manuscript.md).

This standalone manuscript develops the reviewed exact-mesh amplitude criterion and reports actual finite-field numerical experiments. It proves a conditional equivalence under the stated full-field assumptions. It supplies no independent RH upper bound. The older `paper_clean` manuscript remains a separate historical source and is not superseded silently.

## Contents

- `manuscript.md`: authoritative editable manuscript, including methods, results, proofs, and references.
- `manuscript.tex` and `manuscript.pdf`: generated typeset versions.
- `build.sh`: rebuild LaTeX and PDF with Pandoc and pdfLaTeX.
- `numerics/run_experiments.py`: reproducible numerical study; run from any directory.
- `numerics/*.csv`: zero ordinates, all quadrature orders, mesh sensitivity, and matrix scan.
- `numerics/numerical_results.md`: full numerical report.
- `numerics/environment.json` and `requirements.txt`: recorded environment and Python dependencies.
- `check_complementary_mesh.py`: exact rational algebra and finite placement checks.
- `review_manuscript.md` and `review_numerics.md`: separate AI reviews and resolution record.

## Reproduce

Use Python 3 with a virtual environment, then run from this directory:

```sh
python -m pip install -r numerics/requirements.txt
python check_complementary_mesh.py
python numerics/run_experiments.py
bash build.sh
```

The numerical program regenerates the zero table, results, figures, and numerical report deterministically apart from environment metadata and PDF timestamps. The numerical tables in the manuscript are a reviewed snapshot; if parameters or computations change, update those tables and statements before rebuilding the paper. Pandoc and a LaTeX installation are additional system requirements for the PDF build. The generated `.tex` can also be compiled directly with pdfLaTeX.

The main study compares quadrature orders 32/64/128 and uses independent analytic energy integration. A separate reviewer checked one case using adaptive integration. These are ordinary floating-point accuracy checks, not interval-certified bounds. Known critical-line finite sums and synthetic growing signals cannot establish the full-field asymptotic hypothesis.
