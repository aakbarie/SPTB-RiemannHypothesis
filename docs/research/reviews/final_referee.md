# Final adversarial mathematical referee report

Audit date: 10 September 2026. Object: `complementary_spline_observations_revised.md`, including the endpoint and constant-dependence clarifications incorporated during this audit. This is an independent AI mathematical check, not human peer review or formal verification. I read the argument directly and inspected the uploaded SPTB manuscript, especially Definitions 3.1, Sections 4 and 6, and Conjecture 13.3. I did not rely on other agents' audit conclusions.

## Verdict

**I find no fatal mathematical gap remaining in the revised theorem under its stated assumptions:** a fixed full zero field, alpha at least 1, independent affine projections, the exact origin-anchored mesh h(T) = kappa/log T, and polynomial residual bounds for every sufficiently large real T. The four-way equivalence in Section 1 is supported by the proof.

**The original manuscript's broader canonical-regime statement is not fully proved by this argument.** Its Section 4 permits any block width between kappa_1/log T and kappa_2/log T. The revised theorem proves the Horocycle implication for a fixed exact choice kappa/log T within that range. It also proves the implication if the original premise is explicitly uniform over all such admissible choices, since one can then select that exact choice. It does not establish the implication from a bound along an arbitrary single admissible width function. This is a scope limitation, not a fatal gap in the narrowed theorem.

**No independent RH upper bound is supplied.** At sigma = 1/2 this is a criterion equivalent to RH. The reverse bound assumes confinement and cannot provide an unconditional proof of the hypothesis.

## Findings ranked by severity

### High: original-theorem coverage requires an explicit mesh interpretation

The two observations use T and T + kappa/2. Their complementarity follows from the exact identity involving log(1 + kappa/(2T)); comparability of mesh widths does not imply it. The source manuscript's general canonical inequalities do not themselves guarantee the required second observation.

The draft's piecewise-linear example is valid. Interpolate the values exp(n^2) at the nonnegative integers. Every grid of width 1/m, anchored at zero, refines its affine intervals, including the final shortened block. The residual therefore vanishes, although the function has superpolynomial growth. The sparse exact horizons T_m = exp(kappa m), and the all-real width function 1/ceil(log T/kappa), respectively defeat the proposed general reconstruction implication without the exact-mesh/all-real assumptions. These are counterexamples for general functions, not counterexamples for the zeta-specific conjecture.

Required presentation fix: retain the existing scope warning and do not describe this as a proof of the entire originally stated canonical-regime conjecture without specifying its quantifiers. Likewise, the argument does not recover the manuscript's exponential lower bound for every sufficiently large T.

### High: full-field and derivative distinctions must remain explicit

The comparison P_T f against P_(T+kappa/2) f requires the same f. A moving truncation of the amplitude changes f between observations and invalidates the triangle-inequality argument unless additional error estimates are supplied. Truncating only the derivative contribution is harmless for this converse, provided the functional is well-defined and dominates the full amplitude residual.

The uploaded manuscript's Remark 3.2 explicitly truncates the derivative, whereas Definition 3.1 defines the full harmonic field. Thus the revised amplitude interpretation does match that distinction. The source's original absolute-convergence assertion at alpha = 1 is incorrect, but the revised Section 5 supplies a valid replacement. No derivative-variance theorem from the source is needed for the revised converse or the four-way equivalence for E.

### Low, fixed during audit: endpoint and constants

The initial affine piece in the jump reconstruction must intersect J_T in positive length. A piece that touches only its left endpoint could choose the wrong initial polynomial. The revised wording now specifies positive length. Jump traces are correctly assigned to the finite spline P, requiring no pointwise traces of f. The draft now also states that the ultimate growth constant can depend on the initial compact norm and affine coefficients; residuals cannot control a freely added global affine function. These clarifications resolve the issues.

## Verification of the proof transitions

1. **Local jump detection.** Subtracting the left polynomial leaves a step plus a hinge. Projection using the displayed inverse Gram matrix yields the stated matrix, including the cross term. Its determinant is q^4(1-q)^4/12. For q in [1/8,7/8], its trace is positive and below 1, so its minimum eigenvalue is at least its determinant and therefore at least c_*. Under rescaling, the coefficient vector is (d, h'e), and the integration contributes h'. This gives exactly h'|d|^2 + (h')^3|e|^2. Complex coefficients cause no change to this real symmetric Hermitian quadratic form.

2. **Mesh placement.** For a first knot x_j in [T/2,3T/4], its second coordinate is j + theta_j with theta_j between 1/8 and 3/8 for large T. Thus the relevant second block has an interior break bounded away from its endpoints. Since h' < h, it has no other first-mesh knot. Distinct first knots occupy distinct second blocks, and these blocks stay inside [0,T] once T is large. The final shortened blocks at either horizon are irrelevant to these interior observations.

3. **Jump sum and reconstruction.** On each crossing block Q is an affine competitor for the best affine error of P. The square of the norm of P-Q is at most twice the sum of residual energies, restricted to [0,T] when necessary. Summing the local lower bounds over disjoint blocks proves (12), because h'/h is bounded below. The telescoping formula (13) correctly adds a value jump and a slope jump at each crossed knot. With O(T/h) knots, Cauchy-Schwarz and integration over length O(T) give the displayed T^2/h^2 and T^4/h^4 losses. No cancellation, sign, frequency separation, or regularity of f is assumed here.

4. **Affine coefficient propagation.** The intervals for T_n = (4/3)^n T_0 overlap in an interval of length T_n/12. The approximation errors bound the norm of the difference of consecutive affine fits there. Rescaling to the fixed interval [2/3,3/4] proves both coefficient estimates (17), including the intercept estimate despite the interval being far from zero. Taking r = p+5 absorbs the fourth power of log T eventually and makes both geometric-summation exponents positive. Hence the affine fits, then f, have polynomial L2 mass on each interval. These overlapping intervals cover the tail. Summing through the last interval needed to reach R costs only a constant multiple of R^r. The initial compact interval supplies a finite additive constant.

5. **Full field at alpha = 1.** The zero count in each unit ordinate interval is O(log V), with multiplicity. As |rho| is comparable to |gamma| at large height, the total absolute coefficient in frequency bin k is O(log(2+|k|)/(1+|k|)^alpha). The squares of these bin sums are summable for alpha at least 1. Two integrations by parts against a fixed smooth cutoff are uniform in the bounded real exponents eta. Binning replaces the frequency-difference kernel by a constant multiple of (1+|k-l|)^(-2), including neighboring bins. The convolution estimate then bounds every finite tail by the sum of its squared bin masses, which tends to zero. This establishes the required L2-local limit without a minimum zero gap. The proof does not mistake coefficient square-summability alone for sufficient convergence.

6. **Primitive and continuation.** The sum of a_rho/|z_rho| converges by the same zero count, since its high-height terms have exponent alpha+1 at least 2. This gives absolute uniform convergence of the integrated series on compact time intervals; L2-local convergence identifies it with the primitive of H. For Re w sufficiently large, absolute domination permits termwise integration and yields (26). On every compact subset avoiding 0 and the locally finite zero poles, its tails converge uniformly because |w-z_rho| is comparable to |gamma|. Thus the expression genuinely defines a meromorphic function, not just a formal expansion.

   At a putative right-half-plane pole z_0, only zeros with exactly the same beta and positive ordinate contribute to its residue. Their positive weights add. Conjugate poles have negative ordinate and cannot coincide. The residue is therefore nonzero. Polynomial growth of the primitive makes its actual Laplace transform holomorphic throughout Re w > 0. The right half-plane minus a locally finite discrete set is connected, so the identity theorem transports equality from the initial far-right half-plane to a punctured neighborhood of z_0. The holomorphic transform then contradicts the nonzero residue. No maximal real part or isolated spectral edge is needed; only each individual finite-plane zero is isolated.

7. **Reverse O(T) direction.** Under eta <= 0, translation to each unit time interval multiplies each coefficient by e^(eta n) times a unit phase, so its magnitude does not increase. On the fixed cutoff support, all needed derivatives of the real exponential factors remain uniformly bounded, even on the negative part of that support. The same bin argument gives a uniform L2 bound on every unit interval for finite truncations. Passing to the established local limit preserves that bound. Summation gives O(1+T), and orthogonal projection gives E <= the full L2 norm. This proves the reverse direction without using the manuscript's derivative penalty estimates.

## Final assessment

The revised exact-mesh theorem withstands this adversarial check. The essential outstanding limitation is the mismatch between that proved formulation and a possible arbitrary-width reading of the original canonical regime. A new argument would be required to remove that limitation. Establishing any equivalent polynomial bound for the actual field independently of confinement remains a separate task; nothing in this proof establishes RH unconditionally. No novelty assessment is made.
