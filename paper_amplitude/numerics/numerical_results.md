# Numerical experiments

These are reproducible finite-field illustrations of the amplitude residual. They do not estimate the infinite-field truncation error, establish an asymptotic bound, locate unknown zeros, or prove RH. The synthetic growing mode is not claimed to be a zeta zero.

## Method

The first 64 positive critical-line zeros are computed using mpmath.zetazero at 40 decimal digits and saved in zeta_zeros.csv, with numerical |zeta(rho)| checks. The experiments use N=32,64, alpha=1,2, sigma=1/2, kappa=1, and T=8,16,32,64. These ordinates are the zeros returned by that routine; this calculation is not an independent certified zero count.

For each T the origin-anchored blocks have width 1/log(T), including a shortened last block. Gauss-Legendre quadrature independently computes the affine L2 projection on each block in the centered basis 1,x. Squared residuals are integrated directly, avoiding subtraction of two nearly equal total energies. Orders 32,64,128 are all saved. Finite-field total energy is also computed by analytically integrating the pairwise exponential expansion. The same fixed finite field is used across T for each N; no moving cutoff is used within a run.

The largest relative difference in residual between orders 32 and 128 is 1.106e-03 for the zeta fields, showing that the coarsest order is insufficient for the most oscillatory case. The largest relative change in residual from order 64 to 128 is 7.062e-14 for the zeta fields and 7.764e-14 for the synthetic modes. The largest relative discrepancy between order-128 energy and the independent analytic energy is 6.985e-14. These are numerical consistency checks, not rigorous quadrature error certificates.

## Critical-line finite fields

| N | alpha | T | energy / T | residual / T | residual / energy |
|---:|---:|---:|---:|---:|---:|
| 32 | 1 | 8 | 0.0087672826 | 0.0083403512 | 0.95130402 |
| 32 | 1 | 16 | 0.0086993294 | 0.0068227134 | 0.78428038 |
| 32 | 1 | 32 | 0.0087014332 | 0.0060918163 | 0.70009344 |
| 32 | 1 | 64 | 0.0086804624 | 0.0053353887 | 0.61464337 |
| 32 | 2 | 8 | 1.826086e-05 | 1.6128643e-05 | 0.88323569 |
| 32 | 2 | 16 | 1.8451889e-05 | 1.05872e-05 | 0.57377325 |
| 32 | 2 | 32 | 1.8522581e-05 | 7.7791483e-06 | 0.41998187 |
| 32 | 2 | 64 | 1.8468565e-05 | 5.5890661e-06 | 0.30262591 |
| 64 | 1 | 8 | 0.0096355373 | 0.0092109399 | 0.95593423 |
| 64 | 1 | 16 | 0.0095764049 | 0.0076796032 | 0.80192967 |
| 64 | 1 | 32 | 0.0095662129 | 0.0069352109 | 0.72496932 |
| 64 | 1 | 64 | 0.0095416852 | 0.0061880721 | 0.64853032 |
| 64 | 2 | 8 | 1.8312357e-05 | 1.6179481e-05 | 0.88352809 |
| 64 | 2 | 16 | 1.8502559e-05 | 1.062781e-05 | 0.57439675 |
| 64 | 2 | 32 | 1.857237e-05 | 7.8188967e-06 | 0.42099618 |
| 64 | 2 | 64 | 1.8518251e-05 | 5.636686e-06 | 0.30438543 |

### Mesh constant sensitivity

At T=64, N=64 and alpha=1, all other parameters fixed, the following values use order 256. Orders 64,128,256 are retained in kappa_sensitivity.csv. Different mesh constants define different observations; comparisons are finite numerical values, not monotonicity claims for non-nested partitions.

| kappa | h | residual / T |
|---:|---:|---:|
| 0.5 | 0.12022459 | 0.0034560577 |
| 1 | 0.24044917 | 0.0061880721 |
| 2 | 0.48089835 | 0.0087798806 |

Maximum relative residual difference between orders 128 and 256 in the mesh sensitivity experiment: 4.347e-15.

All tested modes have real part 1/2, so neutral exponential factors are built into these finite examples. Their bounded energy density is expected; it provides no evidence about possible zeros outside the sampled set. Differences between N=32 and N=64 measure only this finite cutoff change and are not bounds on the omitted infinite tail.

## Synthetic stress test

Use one mode |1/2+i gamma_1|^{-1} exp(eta t) cos(gamma_1 t), with eta=-0.08,0,+0.08. Only the exponential rate is modified. No synthetic beta is asserted to be the real part of an actual zeta zero.

| eta | T | energy / T | residual / T |
|---:|---:|---:|---:|
| -0.08 | 8 | 0.0014097234 | 0.0010903329 |
| -0.08 | 16 | 0.00090088288 | 0.00039934319 |
| -0.08 | 32 | 0.00048527809 | 0.00011563664 |
| -0.08 | 64 | 0.00024409034 | 3.236089e-05 |
| -0.08 | 128 | 0.00012204954 | 9.5184414e-06 |
| +0.00 | 8 | 0.0024990627 | 0.0018852423 |
| +0.00 | 16 | 0.0024990631 | 0.0010967475 |
| +0.00 | 32 | 0.0024990644 | 0.0005945995 |
| +0.00 | 64 | 0.0024990696 | 0.00033438501 |
| +0.00 | 128 | 0.0024990902 | 0.00019551874 |
| +0.08 | 8 | 0.005069146 | 0.0036927753 |
| +0.08 | 16 | 0.011648507 | 0.005000082 |
| +0.08 | 32 | 0.081132628 | 0.019169081 |
| +0.08 | 64 | 6.8228997 | 0.95952248 |
| +0.08 | 128 | 95377.031 | 7073.4397 |

The growing mode eventually dominates the shrinking block width in the sampled range; the experiment illustrates a mechanism, not an every-T lower bound for an infinite superposition.

## Local matrix check

For 1001 equally spaced q in [1/8,7/8], the smallest observed eigenvalue of M(q) is 0.0001616834055. The conservative analytic lower bound (7/64)^4/12 is 1.192589601e-05. The sampled minimum is not a certified continuum minimum; positivity on the continuum follows from the exact determinant proof in the paper.

## Reproduction and limitations

Run `python run_experiments.py` with numpy, scipy, mpmath, and matplotlib installed. Exact package and platform versions are recorded in environment.json. CSVs retain all quadrature orders; the figures show order 128. No randomness is used. The computations are ordinary floating point after high-precision zero generation. No interval arithmetic, full-field tail certification, asymptotic rate fit, or numerical proof is claimed.
