# Paper III — information dimension and noncommuting observation limits (iteration 16)

The consolidated PDF and LaTeX addendum are authored outside the repository; this note records the exact theorem and its scope.

Let p₀,…,p_d be pairwise distinct primes, d≥2. Define L_j = log(p_j/p₀) and Φ_d(t)=(exp(-itL_j))_j∈T^d. Let ν_B be the pushforward of uniform measure on [0,B], and H_B(n) its Shannon entropy in the n^d equal angular cells.

**Theorem (written proof, NOT Lean formalized).**

* For every fixed 0<B<∞, H_B(n)=log(n)+O_B(1) as n→∞, hence D₁(ν_B)=1.
* For fixed n≥2, ν_B converges weakly to Haar measure on T^d and H_B(n)→d log(n) as B→∞.
* The ordered limits are 1 and d:
  lim_{B→∞} lim_{n→∞} H_B(n)/log n = 1,
  lim_{n→∞} lim_{B→∞} H_B(n)/log n = d.
* Any fixed finite real β introduces positive coordinate radii exp(-βL_j), a bi-Lipschitz map; the information dimensions remain 1 and d.
* KL divergence to uniform cell distribution is d log(n) − H_B(n).

**Proof mechanism.** Unique prime factorization implies the L_j are rationally independent. At fixed B the smooth trajectory crosses O_B(n) grid cells; each cell receives at most O_B(1/n) mass, so H_B(n)=log n+O_B(1). For long B, every nontrivial integer Fourier character has average bounded by 2/[B|k·L|], giving Haar equidistribution and convergence of masses of fixed cells.

**Limits and nonclaims.** These are classical information-dimension/Kronecker-Weyl techniques applied to prime-phase observations. The continuous parameter t is NOT a sample of discrete zeta-zero ordinates. D₁ does not distinguish on-line from off-line real spectral parts. No compatible arithmetic spectral metric or RH result follows. The entropy-limit proof is NOT Lean checked. A finite-dimensional radial scaling distinction lemma (`nonzero_radial_scaling_preserves_phase_distinction`) passed Lean compilation and standard-axiom audit in [Actions run 44881142589](https://github.com/aakbarie/SPTB-RiemannHypothesis/actions/runs/44881142589). This does not formalize the Shannon entropy limits.

**Relation to Papers I–III.** Paper I: distinction changes with resolution and reference. Paper II: phase coherence. Paper III: finite-horizon observability and compatibility. The dimensional noncommutation isolates the horizon-resolution operation and does not replace the missing arithmetic operator-symmetry law.
