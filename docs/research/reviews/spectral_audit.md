# Independent spectral audit and proposed spectral lemmas

Audit target: Sections 5–6 of `output/complementary_spline_observations.md`.
The geometric reconstruction estimate is assumed here; this audit does not certify that separate argument.

## Verdict

The local convergence, primitive, and Laplace-pole argument are valid under the stated fixed-full-field assumptions. I found no fatal gap in Sections 5–6. The claimed alpha >= 1 range is conservative: these arguments work for alpha > 1/2. Retaining alpha >= 1 in the main paper avoids changing the paper's scope unnecessarily.

The proof does not require a positive minimum gap between ordinates, a rightmost zero, a spectral edge attained by a zero, or RH. It does require a locally finite zero multiset and a nonzero aggregate coefficient at each distinct spectral exponent. For the stated positive weights, multiplicities satisfy the latter automatically.

## 1. Standard arithmetic input

Write positive-ordinate nontrivial zeros as rho=beta+i gamma, counted with multiplicity. The unconditional zero count gives

    N(V+1)-N(V) <= C log(V+2).

For a directly verified primary source, Timothy Trudgian, *An improved upper bound for the argument of the Riemann zeta-function on the critical line II*, arXiv:1208.5846v2, Corollary 1, gives the Riemann–von Mangoldt main term with an explicit O(log V) error. Subtracting its estimates at V and V+1 gives the displayed unit-bin bound. Its theorem does not assume RH. Endpoint conventions are harmless: either use a slightly enlarged unit interval or include the endpoint multiplicities, which the same estimate controls.

Source: https://arxiv.org/pdf/1208.5846

## 2. Local L2 construction

Set eta_rho=beta-sigma and a_rho=|rho|^(-alpha). Expand each cosine into frequencies +gamma and -gamma with coefficients a_rho/2. For signed frequency bins [k,k+1), let A_k be the sum of absolute coefficients. Then

    A_k <= C log(2+|k|)/(1+|k|)^alpha.

Consequently sum_k A_k^2 is finite for alpha>1/2. All eta_rho are in one bounded interval.

Fix chi smooth and compactly supported on the real line, equal to one on the compact observation interval. For any two signed modes of frequencies omega and omega', integration by parts twice bounds their localized inner product by

    | integral chi(t)^2 exp((eta+eta')t) exp(i(omega-omega')t) dt |
      <= C_chi (1+|omega-omega'|)^(-2).

The constant is uniform over the zeros because eta and eta' remain bounded. If omega and omega' belong to bins k and l, then

    (1+|omega-omega'|)^(-2) <= 4(1+|k-l|)^(-2).

For |k-l|>=2 this follows from |omega-omega'|>=|k-l|-1; for the remaining cases increase the constant to four. Thus every finite tail satisfies

    ||chi H_tail||_2^2 <= C sum_{k,l} A_k^tail A_l^tail (1+|k-l|)^(-2)
                       <= C ||(1+|k|)^(-2)||_ell1 sum_k (A_k^tail)^2.

The final bound is Young's inequality on ell2 (or Cauchy–Schwarz at each fixed difference). The right side tends to zero when the lower ordinate cutoff goes to infinity, including partially filled endpoint bins. Hence the ordinate truncations converge in L2 on every compact interval. This establishes the fixed field H without a spacing assumption.

## 3. Primitive and meromorphic transform

Let z_rho=eta_rho+i gamma. Since gamma>0, z_rho is nonzero. The unit-bin count implies

    sum_rho a_rho/|z_rho| < infinity

for every alpha>0. Therefore

    K(t)=1/2 sum_rho a_rho [ (exp(z_rho t)-1)/z_rho
                            +(exp(conj(z_rho)t)-1)/conj(z_rho) ]

converges absolutely and uniformly on compact t-intervals. For finite truncations, K_X(t)=integral_0^t H_X(u)du. L2 convergence gives uniform convergence of these primitives on [0,R], since

    sup_{0<=t<=R} |integral_0^t(H_X-H)(u)du| <= sqrt(R)||H_X-H||_L2(0,R).

Thus K is exactly the primitive of the L2 field. No pointwise convergence or termwise differentiation is required.

For Re w>max(1-sigma,0), Tonelli's absolute integrability criterion, using the above summability, permits integration term by term. Elementary integration gives

    L K(w) = (1/(2w)) sum_rho a_rho [1/(w-z_rho)+1/(w-conj(z_rho))].

The series on the right is normally convergent on compact sets avoiding zero and the spectral points: sufficiently high ordinates have |w-z_rho|>=gamma/2 uniformly on such a compact set, and the remaining finite terms are elementary rational functions. It therefore defines a meromorphic function there.

At a spectral point z_0 with Re z_0>0 its residue is

    (1/(2z_0)) sum_{rho:z_rho=z_0} a_rho != 0.

Only finitely many terms share z_0. Other points cannot accumulate locally, and their sum is analytic nearby. The negative-frequency terms have negative imaginary part and cannot contribute to this positive-imaginary pole. Repeated zeros add positive weights. Zeros with the same ordinate but different real parts have different spectral points.

If integral_0^R|H|^2 <= C(1+R)^p, then |K(t)|<=C'(1+t)^((p+1)/2). Hence its actual Laplace transform is holomorphic on Re w>0: on every compact sub-half-plane, the integrand and all its w derivatives have an integrable exponential majorant. The punctured right half-plane is connected because the removed set is locally finite. The identity theorem identifies the actual transform with the meromorphic expression throughout that punctured domain. A holomorphic extension over z_0 contradicts the computed nonzero residue. Therefore beta<=sigma for all zeros.

This proves the spectral implication needed by the draft.

## 4. Useful converse: a linear L2 bound under confinement

Assume eta_rho<=0 for all zeros. Fix chi in C_c^infinity((-1,2)), equal to one on [0,1]. On [n,n+1] put t=n+s. Each mode has coefficient (a_rho/2)exp(eta_rho n)exp(+-i gamma n). Its absolute value is at most a_rho/2 for n>=0. The same cutoff estimate in s has one uniform constant, since eta_rho ranges over a fixed bounded interval and s is restricted to [-1,2]. Therefore

    integral_n^(n+1) |H(t)|^2 dt <= C sum_k A_k^2,

uniformly for all nonnegative integers n. Passing from finite sums to H uses local L2 convergence. Summing the unit intervals gives

    integral_0^R |H(t)|^2 dt <= C(1+R).

Because zero is an allowed affine fit, E_H(T)<=integral_0^T|H|^2. Thus confinement implies E_H(T)=O(1+T) for *any* partition; the hard direction needs the mesh theorem.

Combining this with the assumed complementary-mesh reconstruction yields the clean criterion:

    all beta<=sigma
      <=> integral_0^R|H|^2=O(1+R)
      <=> E_H(T)=O(1+T)
      <=> E_H(T)=O((1+T)^p) for some finite p and all large real T,

where the E statements use the specified canonical partition. The middle-to-left implications use the mesh theorem and the spectral argument above. This is a conditional equivalence criterion, not an unconditional proof of confinement or RH.

## 5. Optional exact exponential rate for the full-field energy

Let B=sup_rho eta_rho. The same translated-window estimate, now using exp(eta_rho n)<=exp(Bn), gives

    integral_n^(n+1)|H|^2 <= C exp(2Bn).

Consequently

    limsup_{R->infinity} [log(1+integral_0^R|H|^2)]/R <= 2 max(B,0).

For B>0, suppose this limsup were smaller than 2B. Choose b>0 with the limsup <2b<2B. Then the integrated energy is O(exp(2bt)), and |K(t)|<=C sqrt(t)exp(bt). Its Laplace transform is holomorphic on Re w>b. By the definition of supremum there exists a zero with eta_rho>b. The same nonzero-pole contradiction applies on this shifted half-plane. Thus

    limsup_{R->infinity} [log(1+integral_0^R|H|^2)]/R = 2 max(B,0).

For B<=0, the previous upper bound and nonnegativity immediately give equality with zero. This is an exact *limsup* theorem for the full-field integrated energy. It is not a lower bound for every large R and is not yet an exact exponential-rate theorem for the spline residual.

## 6. Limits of this audit

- No moving amplitude cutoff is covered.
- No full untruncated derivative at alpha=1 is claimed to be L2.
- Positivity is used only to ensure nonzero aggregate residues at coincident exponents; generic amplitude kernels with opposite coefficients do not inherit this conclusion.
- The affine-mesh reconstruction is an independent proof obligation.
- Novelty has not been established by a literature search.
- The argument remains a mathematical derivation audited by another model in the same system, not external peer review.
