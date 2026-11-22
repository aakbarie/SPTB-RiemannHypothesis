# Witness Functors: Representational Limits, Approximation, and Certified Refinement

Revised research manuscript by Akbar Akbari Esfahani, 21 November 2026.

Read [the PDF](main.pdf), [the LaTeX source](main.tex), or the
[revision record](REVISION_NOTES.md).

The paper separates a target distinction lost by an observation from error
introduced by a restricted model. It connects finite observation refinements,
Hilbert witnesses, an explicit enriched distance distributor, and the proof
obligations needed for SPTB and refinement bridges. Its interpretation concerns
how the choices governing a representation can themselves become observable.

## Euler product example

Crossley’s constructed Euler products give an obstruction pair for the paper’s fibre criterion. Assuming RH for $\zeta$ and the critical-line property for the quadratic characters of conductors 5, 13 and 65, the model permits either a double zero at $1/2$ or two simple zeros at $1/2\pm i\delta$. Choose a circle $s=1/2+Re^{it}$ enclosing just those zeros, with $R>\delta$. In both cases the coefficient of $e^{-it}$ in $-L'/L$ is $-2/R$. Expanding the pole terms gives $0$ and $2\delta^2/R^3$, respectively, for the coefficient of $e^{-3it}$. That third coefficient distinguishes the two inputs. [Crossley, *A prime disagreement model with two prescribed zero patterns*](https://doi.org/10.5281/zenodo.24417708).

## Build

With a working TeX Live installation containing latexmk, Latin Modern, AMS,
mathtools, geometry, microtype, booktabs, hyperref, xcolor, and enumitem:

```sh
cd paper_witness
make all
```

The build does not require shell escape. The checked PDF is committed for
readers without LaTeX. The mathematical content is a research draft with written
proofs, not a Lean formalization or an independent referee certification.

## Verification

```sh
python3 -m pip install -r requirements.txt
make check
python3 verify_witness.py --certificate example.json
```

The finite certificate mode uses only the Python standard library. It requires
explicit state labels, positive rational weights summing to one, observation
labels, target values, diagnostic values, and diagnostic provenance. It computes
fiber projections, obstruction-pair bounds, strictness, and the exact refinement
identity. It cannot recover a target or diagnostic from observations that omit
them. Labels must be hashable JSON scalars.

The checks enumerate finite instances with rational arithmetic, independently
solve the polynomial projection example symbolically, and evaluate the spectral
scaling identities by quadrature. See [the verification record](VERIFICATION.md).
They test examples and regression cases; they do not prove general theorems.

## Scope and source lineage

This directory adds the revised witness manuscript to the SPTB repository.
Its sources are the author's original ten-page Witness Functors PDF and the
21 November reconstruction and bridge/SPTB audits. The revision reworks those
sources rather than presenting the original categorical claims as established.
The earlier manuscripts are identified in the bibliography; the SPTB source
snapshot is pinned there by commit SHA.

The SPTB source in `paper_clean/` retains its own analytic proof obligations.
The Paper 3 Lean draft in PR #5 is a separate project. This revision does not
modify either manuscript or import the Paper 3 branch. A general scientific-paper
compiler is specified as future work; the finite checker is the implemented
example, not an automated compiler of arbitrary scientific theories.
