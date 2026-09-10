# Independent numerical implementation review

Reviewed files: `numerics/run_experiments.py`, all numerical CSVs, `numerics/numerical_results.md`, `numerics/environment.json`, and `numerics/requirements.txt`.

Review date: 10 September 2026. This is an AI review with independent numerical spot checks, not interval certification or formal verification.

## Verdict

The numerical implementation and reported conclusions are accepted within their stated ordinary floating-point scope. No substantive defect was found. The computations illustrate finite fields and finite horizons; they do not establish the infinite-field asymptotic hypothesis of the analytic theorem.

## Projection and integration

The block construction uses the exact requested formula `h = kappa / log(T)`, anchors the mesh at zero, and clips the final block at the horizon. For the centered reference coordinate `x` in `[-1,1]`, the affine Gram matrix is diagonal with entries `2` and `2/3`. Therefore the code's coefficients

`c0 = sum(w*f)/2`, `c1 = 3*sum(w*x*f)/2`

are the correct quadrature approximations to the L2 affine projection coefficients. The physical Jacobian is constant on each block and cancels from these projection equations. Multiplication by the half-width in the subsequent norm integrals is correct. Direct integration of the squared residual avoids subtracting nearly equal large global energies.

For the real fields used here, expanding the cosines into exponentials gives coefficients `a/2` at exponents `eta +/- i*gamma`. The pairwise sum of exponents, rather than their conjugate differences, correctly integrates the square of this real field. The `expm1(s*T)/s` formula is valid, with limit `T` at `s=0`. All actual experiments have exactly zero or safely nonzero pairwise exponents for the branch used in the code. The routine is not being certified as a general-purpose near-zero complex-exponent integrator.

## Independent checks performed

- Parsed 48 critical-line field rows, 45 synthetic-mode rows, nine mesh-sensitivity rows, and 1,001 local matrix rows.
- Checked nonnegative residuals bounded by the total computed energy for all field, synthetic, and sensitivity rows.
- Checked energy normalization by the horizon.
- Recomputed the reported quadrature-order comparisons directly from the CSVs. Differences in the last printed digits below arise from equivalent floating-point formulas for relative differences.
- Recomputed the sampled eigenvalue minimum and determinant agreement from the matrix CSV.
- Independently integrated the most oscillatory critical-line case at the shortest horizon, `N=64`, `alpha=1`, `T=8`, using adaptive `scipy.integrate.quad_vec` on every block. This used physical centered moments and the separate formula `I2 - I0^2/l - 12*I1^2/l^3`, rather than the production Gauss-Legendre residual formula. The result was total energy `0.07708429842447534` and residual `0.07368751945583894`. Relative differences from the stored order-128 values were at most `7.84e-14`.
- Independently checked every order-128 synthetic analytic energy against the real single-mode identity obtained from `cos^2(gamma*t) = (1+cos(2*gamma*t))/2`. Maximum relative discrepancy was `2.23e-16`.

| Check | Independently recomputed result |
|---|---:|
| Largest critical-line residual difference, orders 32 and 128 | 0.00110603808045 |
| Largest critical-line residual difference, orders 64 and 128 | approximately 7.06e-14 |
| Largest synthetic residual difference, orders 64 and 128 | approximately 7.76e-14 |
| Largest mesh-sensitivity residual difference, orders 128 and 256 | approximately 4.4e-15 |
| Largest order-128 energy discrepancy from analytic formula | 6.98531579e-14 |
| Smallest sampled matrix eigenvalue | 0.000161683405534 |
| Largest sampled determinant relative discrepancy | 2.45e-15 or less |

The original report correctly identifies the order-32 computation as insufficiently resolved in the most oscillatory case. Its accepted results use higher orders. Agreement between higher quadrature orders is a numerical consistency check, and the report correctly avoids calling it a rigorous error bound.

## Interpretation and labels

The critical-line examples retain a fixed finite zero list across horizons within each run. The code uses real part one half and zero exponential rate in those examples, so bounded finite-field energy density is expected. The report explicitly states that doubling the cutoff does not estimate or bound the entire omitted infinite tail.

The synthetic tests modify only the exponential rate of a single mode while retaining its frequency and amplitude. They are labelled synthetic and are not represented as newly found off-line zeta zeros. Their contrasting damped, neutral, and growing behavior supports only the illustrated finite-mode mechanism. It is not an every-horizon lower bound for an infinite superposition.

The local matrix scan is correctly distinguished from the exact continuum coercivity proof. The sampled eigenvalue minimum is not presented as a certified global minimum.

The mesh-sensitivity comparisons are correctly treated as finite observations rather than a monotonicity theorem for arbitrary partitions. The documented package versions match the environment record.

## Remaining limitations

No numerical repair is required for the present claims. The computations do not use interval arithmetic, independently certify a complete zero count, bound the infinite tail, test arbitrarily large horizons, or numerically establish the complementary-mesh reconstruction inequality for all L2 functions. They should remain labelled finite-field illustrations, exactly as the current report states.

## Final manuscript addition

The added fixed-finite-field Taylor estimate has the correct coefficient `1/320`. With `M = ||H_N''||_infinity`, midpoint linear approximation has pointwise error at most `M |t-m|^2 / 2`. Squaring and integrating on a block of length `l` gives `M^2 l^5 / 320`. Since the optimal affine projection has no larger error, and `sum l^5 <= h^4 sum l = h^4 T`, the claimed total estimate follows. The derivative norm depends on the cutoff, as the manuscript explicitly states; no uniform infinite-field bound follows.

## Final production check

The final PDF was rebuilt after layout review. Table 1 is now kept together and Table 2 remains with its caption. All 14 final pages were rendered and visually inspected, with no remaining clipping or overlap observed. The Taylor coefficient 1/320 was separately checked during final assembly.
