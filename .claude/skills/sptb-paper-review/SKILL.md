---
name: sptb-paper-review
description: Drives a 4-dimension adversarial review (rigor, literature, numerics, exposition) of the Spline–Penalized Tail Bound (SPTB) Riemann-Hypothesis paper in this repo, synthesizes prioritized findings, then applies mechanical/Tier-1 fixes, re-runs the R numerics, and rebuilds the PDF in a loop. Use when the user asks to review, critique, audit, referee, or revise the SPTB / Riemann Hypothesis / "L-functions tail bound" paper here.
---

# SPTB Paper Review & Revise

## Overview
This skill runs a full-pipeline, multi-agent adversarial review-and-revise loop on the
**Spline–Penalized Tail Bound (SPTB)** manuscript — an RH-adjacent paper claiming a
"detection theorem" for off-critical zeros of automorphic L-functions. It (1) dispatches
four skeptical review subagents in parallel, (2) synthesizes a prioritized findings report,
(3) applies agreed mechanical/Tier-1 fixes, (4) re-runs the R numerics + rebuilds the PDF,
and (5) loops until paper tables match the committed CSVs and all refs resolve.

This is the user's own thesis work. Be rigorous but constructive.

## When to use
Trigger when the user asks to review, referee, audit, critique, check, or revise the SPTB /
Riemann Hypothesis paper in this repo — e.g. "review my RH paper", "audit the SPTB numerics",
"is the detection theorem sound?", "get this ready for submission".

## Repository map (verified paths — use these exactly)
- **CANONICAL manuscript root**: `paper_clean/main.tex` (amsart, 12pt; March 2026 version).
  `\input`s `parts/part1.tex`..`part4.tex`, `parts/part2b_explicit_formula.tex` (the
  explicit-formula/Turán–Pintz reframe), then `appendices/appA.tex`..`appD.tex`
  (appD IS included — the "for L-Functions" title is supported).
  - Bibliography is a manual `\begin{thebibliography}{99}` with ~21 entries, **cited in the body**.
  - Clean page-break policy (`\partbreak` between major parts only — no per-section `\clearpage`).
    `\numberwithin{equation}{section}` is active.
  - `c_der` is `4\pi^2` (corrected from the erroneous `1/12`); the variance table in part3
    uses the regenerated monotone `variance_table.csv` numbers.
- **Numerics**: `numerical_analysis/sptb_analysis.R`. Outputs land in
  `numerical_analysis/sptb_out/*.csv` and `*.png`. Figures consumed by the paper are in
  `paper_clean/figs/`.
- **Data**: `data/zeros_1.txt` (~1.8 MB, ~100k Odlyzko ordinates). `data/readme.md` OVERSTATES
  this as ~50 MB / 300k zeros — flag the mismatch.
- **Archived (do NOT edit)**: the older Oct-2025 copy was moved to
  `archive/paper-legacy-2025-10/` (formerly `paper/`). It is superseded; never build or cite from it.
  The failure-mode checklist below records the ORIGINAL audit findings against that legacy copy —
  many (conjecture* error, orphaned appD, 6-entry uncited bib, blank-page hack, c_der=1/12,
  fabricated variance table) are ALREADY RESOLVED in `paper_clean/`. Re-verify before re-reporting.

### How the R script actually runs (verified)
`sptb_analysis.R` is a **top-to-bottom script with NO `main()` guard** — sourcing it executes
every stage. It reads `data/zeros_1.txt` and writes to `numerical_analysis/sptb_out/` via
**relative paths**, so it MUST be invoked from the repo root:
```
Rscript numerical_analysis/sptb_analysis.R
```
(`Rscript` is at `/usr/local/bin/Rscript`.) Stage → output CSV map:
- §6/§10 variance sweep → `variance_table.csv` (cols: Tmax, Delta, lambda, blocks, sum_eps2,
  sum_d2, F_lambda, T_logT_loglogT, T_log2T). Grid is `cfg$variance_T_grid =
  c(2e2,5e2,1e3,2e3,5e3,1e4)` — **max T = 1e4, there is NO T=5e4 row.**
- §9 bias sweep → `bias_summary.csv`, `bias_blocks.csv` (over `bias_eta_list=c(1e-4,5e-4,1e-3)`
  × `bias_T_list=c(1e4,2e4,5e4)`; note η=1e-5 is NOT in the grid).
- §13 λ-robustness → `lambda_robustness.csv`. §14 multi-zero → `multi_zero_test.csv`.
- §15 `verify_c_der` → `c_der_verification.csv` (poincare_lb=12/Δ², sharp_constant=4π²/Δ²;
  committed `empirical_min ≈ 39.48 ≈ 4π²`, NOT 1/12). §16 → `C0_bootstrap.csv`.
- Plots saved by `cfg$save_plots`; tables by `cfg$save_tables`.
Key config knobs (`cfg$...`): `zeros_path`, `zeros_n_max`, `sigma=0.5`, `alpha=1.0`,
`gamma_max=200`, `zeros_cap=6000`, `variance_T_grid`, `bias_eta_list`, `bias_T_list`,
`precision` (double|mpfr), `use_huber=TRUE`, `seed=7`.

### Build (verified)
From `paper_clean/`: `make all` runs `latexmk -pdf -interaction=nonstopmode -shell-escape main.tex`.
Also `make clean`, `make distclean`. Check `paper_clean/main.log` for undefined refs / errors.

## The four review dimensions (dispatch as subagents IN PARALLEL)
Send all four Agent/Task calls in a SINGLE message. Each returns a structured report
(findings list, each tagged severity + file:line + suggested fix). Paste these prompts:

### Agent 1 — Mathematical rigor (skeptical analytic number theorist)
> Review `paper_clean/parts/part1.tex`..`part4.tex` and `appendices/appA-C.tex` of an RH paper for
> mathematical rigor. Be adversarial. Specifically determine: (a) Is the kernel
> `H_σ(t)=Σ_ρ |ρ|^{-α} e^{(β−σ)t} cos(γt)` ever DERIVED from ζ (explicit formula / Hadamard
> product) or just posited ad hoc? (b) Is the "detection theorem" CIRCULAR — i.e. is the
> off-line term injected by hand so its exponential growth is tautological rather than forced
> by ζ? (c) Convergence: at σ=1/2, α=1 the amplitude sum behaves like Σ log n / n which
> DIVERGES; the derivative H_σ' needs α>2. Check whether the paper acknowledges this.
> (d) Montgomery–Vaughan: the quoted MV large-sieve inequality is for n^{-it} with
> well-separated frequencies; ζ ordinates are NOT well-separated — is MV transplanted
> illegitimately? (e) Locate the cited "Lemma 3, Part 1" — it does not appear to exist.
> (f) The constant `c_der = 1/12`: this is the variance-of-a-linear-function constant; the
> sharp affine derivative-Rayleigh constant is ~4π²/Δ² (≈π²/Δ² scaling), NOT 1/12 — verify and
> flag. (g) Check the Step-5 "assembly" for sign errors / dropped cross-terms. Report each
> issue with file:line, severity (blocker/major/minor), and a concrete fix.

### Agent 2 — Literature & novelty (RH-literature expert; use WebSearch)
> Assess novelty and citations of an RH paper claiming an SPTB "detection theorem". Use web
> search. Check: (a) The core mechanism (off-line zero ⇒ exponential oscillation growth) is
> essentially the classical **Turán–Pintz Ω-theorem** repackaged — confirm and cite. (b)
> Missing antecedents to engage: Nyman–Beurling / Báez-Duarte L² criterion, Li's criterion,
> Lagarias, Speiser, Bombieri–Lagarias, Rudnick–Sarnak. (c) The "Horocycle Conjecture" name
> collides with Sarnak (1981) and Zagier closed-horocycle work — the paper must engage or
> rename. (d) Bibliography has only 6 entries and the body cites almost nothing — list the
> claims that need citations. (e) Balazs–Vörös appears mischaracterized; Berry–Keating is
> cited to the wrong venue (the canonical reference is SIAM Review 41 (1999) 236–266). Report
> with severity and the correct citation where applicable.

### Agent 3 — Numerical reproducibility (run the R code from repo ROOT if R is available)
> Audit reproducibility of an RH paper's numerics. From the repo root run
> `Rscript numerical_analysis/sptb_analysis.R` and compare outputs in
> `numerical_analysis/sptb_out/*.csv` against the numbers stated in `paper_clean/parts/*.tex`.
> Check: (a) The bias test is CIRCULAR — it injects e^{ηt} then regresses log-energy and
> "recovers" slope 2η; that's an algebraic identity, not validation. (b) The variance table in
> Part 3 does NOT match `variance_table.csv` — wrong magnitude (~5 orders off), wrong
> monotonicity, and the paper has a T=5×10⁴ row the grid never produces (max T=1e4). (c)
> `c_der_verification.csv` has empirical_min ≈ 39.48 ≈ 4π², NOT 1/12 — the paper's 1/12 is
> wrong. (d) Methodology mismatch: paper says "Clenshaw–Curtis to 1e-10" but the code uses the
> trapezoid rule (`trapz`); paper says "first 10⁵ zeros" but the code uses only γ≤200
> (`gamma_max=200`). (e) The "0.001%" headline cherry-picks the η=10⁻⁴ row; η=10⁻⁵ is ~1.5%
> off and isn't even in the CSVs. (f) `C0_bootstrap.csv` CI spans ~2 orders of magnitude.
> Report each as a table-vs-CSV diff with the exact CSV value.

### Agent 4 — Exposition & journal-readiness (math journal editor)
> Copy-edit an RH paper for submission. Check `paper_clean/main.tex` and inputs. Flag: (a) appD.tex
> is ORPHANED (not `\input`) so the title's "for L-Functions" claim is unsupported; (b)
> `thm:bias` has NO proof block; (c) broken `\ref`s (`eq:worst`, `fi1`, `fig4`); (d) duplicate
> `\label{sec:geom-preview}`; (e) manual `\tag{8.x}` conflicts with `\numberwithin`; (f) the
> forced `\clearpage` per section produces ~30 blank pages; (g) abstract appears duplicated;
> (h) author email mismatch (gmail vs thealliance.health); (i) any cover letter targets Adv.
> Math. (unrealistic) and repeats 0.001%; (j) α is never defined; (k) the σ-range is stated
> inconsistently ([1/2,1) vs (1/2,1) vs ≥1/2); (l) 6-entry bibliography is uncited and
> `references.bib` is unused. Report each with file:line and the minimal fix.

## Full pipeline (orchestrate this)
1. **Review** — dispatch all 4 agents in parallel (one message, four calls). Collect reports.
2. **Synthesize** — merge into ONE prioritized report, tiered:
   - **Tier 0 (reframe/conceptual)**: derivation gap, circularity of the detection theorem,
     Turán–Pintz framing, non-equivalence honesty.
   - **Tier 1 (integrity)**: variance table vs CSV, the 0.001% claim, c_der=1/12 vs 4π²,
     honesty about the circular bias test.
   - **Tier 2 (math gaps)**: convergence at α=1, MV misapplication, missing "Lemma 3".
   - **Tier 3 (literature + mechanics)**: citations, broken refs, orphaned appD, blank pages,
     email/cover-letter, bib.
   Present this to the user and **PAUSE for confirmation before editing.**
3. **Revise** — apply the agreed mechanical/Tier-1 fixes by editing `paper_clean/parts/*`,
   `paper_clean/appendices/*`, `paper_clean/main.tex`. Keep `paper_clean/` in sync or flag the divergence.
   **NEVER fabricate numbers** — pull every table value from the committed CSVs, not from prose.
4. **Re-run numerics** — from repo root: `Rscript numerical_analysis/sptb_analysis.R`.
   Regenerate CSVs + PNGs. If a needed package is missing, report it; do not fake outputs.
5. **Rebuild PDF** — from `paper_clean/`: `make all`. Read `paper_clean/main.log` for undefined refs /
   errors. (Mirror to `paper_clean/` only if the user wants it.)
6. **Loop** — re-verify that paper tables match the regenerated CSVs and all `\ref`s resolve.
   Repeat 3–5 until clean. Report a diff of what changed each pass.

## Known failure modes (project-specific checklist)
- [ ] Kernel `H_σ` is posited, never derived from ζ (explicit formula / Hadamard).
- [ ] Detection theorem is circular: off-line term injected by hand ⇒ growth is tautological.
- [ ] Convergence ignored: σ=1/2, α=1 ⇒ Σ log n/n diverges; H_σ' needs α>2.
- [ ] Montgomery–Vaughan transplanted from n^{-it} to non-separated ζ ordinates.
- [ ] Cited "Lemma 3, Part 1" does not exist.
- [ ] `c_der = 1/12` is WRONG; CSV shows ≈39.48 ≈ 4π² (sharp constant 4π²/Δ²).
- [ ] Step-5 assembly: check sign / cross-term validity.
- [ ] Core result = Turán–Pintz Ω-theorem repackaged; cite it.
- [ ] Missing: Nyman–Beurling/Báez-Duarte, Li, Lagarias, Speiser, Bombieri–Lagarias, Rudnick–Sarnak.
- [ ] "Horocycle" collides with Sarnak 1981 / Zagier — engage or rename.
- [ ] Bibliography 6 entries, body cites nothing; `references.bib` unused.
- [ ] Balazs–Vörös mischaracterized; Berry–Keating wrong venue (SIAM Review 1999 is canonical).
- [ ] Bias test circular (inject e^{ηt}, recover slope 2η — algebraic identity).
- [ ] Variance table in Part 3 ≠ `variance_table.csv` (~5 orders off, wrong monotonicity,
      phantom T=5×10⁴ row; grid max is T=1e4).
- [ ] Methodology: paper "Clenshaw–Curtis 1e-10" but code uses trapezoid; "first 10⁵ zeros"
      but code uses γ≤200.
- [ ] "0.001%" cherry-picks η=10⁻⁴; η=10⁻⁵ is ~1.5% off and not in the CSVs.
- [ ] `C0_bootstrap.csv` CI spans ~2 orders of magnitude.
- [ ] appD.tex orphaned ⇒ title "for L-Functions" unsupported.
- [ ] `thm:bias` has no proof block.
- [ ] Broken refs: `eq:worst`, `fi1`, `fig4`. Duplicate label `sec:geom-preview`.
- [ ] `\tag{8.x}` conflicts with `\numberwithin`. Forced `\clearpage` ⇒ ~30 blank pages.
- [ ] Abstract duplicated. Author email mismatch (gmail vs thealliance.health).
- [ ] Cover letter targets Adv. Math. (unrealistic) and repeats 0.001%.
- [ ] α never defined. σ-range inconsistent ([1/2,1) vs (1/2,1) vs ≥1/2).
- [ ] `data/readme.md` overstates zeros file as ~50MB/300k (actual ~1.8MB/~100k, γ≤200 used).

## Guardrails
- This is the user's own thesis. Rigorous but constructive; offer fixes, not just verdicts.
- **Integrity first**: every table/number in the paper must be regenerable from the committed
  CSVs/code. Flag mismatches; NEVER invent or "smooth" a number.
- Faithfully distinguish **PROVED vs CONJECTURED vs HEURISTIC**. Preserve the paper's honest
  non-equivalence disclaimer (it claims a detector, not an RH equivalence) — do not overclaim.
- Do NOT assert the PDF builds or the numerics pass unless you actually ran them and saw it
  (check `paper_clean/main.log`; read the R console output / regenerated CSVs).
- Pause for user confirmation after synthesis (step 2) before applying edits.
