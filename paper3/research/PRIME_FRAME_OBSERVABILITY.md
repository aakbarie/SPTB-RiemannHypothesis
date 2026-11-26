# Paper III: prime-sampled arithmetic observability (iteration 14)

**Status:** Mathematical addendum with written proofs, **not** a proof of RH. The prime-sampling results below are not yet Lean-checked. Six prior generic metric-rigidity lemmas passed Lean build and standard-axiom audit in [Actions run 44878323066](https://github.com/aakbarie/SPTB-RiemannHypothesis/actions/runs/44878323066).

For distinct complex parameters \(s_1,\ldots,s_m\), observe
\[
  ({\cal O}_P c)_p=\sum_{j=1}^m c_j e^{-s_j\log p},\quad p\in P,
  \qquad G_P={\cal O}_P^*{\cal O}_P.
\]
The Gram is automatically positive semidefinite. Strict positivity requires injectivity.

## Three-prime theorem
For \(s_1\ne s_2\), the samples \(P=\{2,3,5\}\) yield \(G_P\succ0\). Otherwise a nonzero kernel vector would give \(e^{-(s_1-s_2)\log(3/2)}=e^{-(s_1-s_2)\log(5/2)}=1\). Taking moduli forces equal real parts; the imaginary part then forces \(\log(5/2)/\log(3/2)\in\mathbb Q\), contradicting unique prime factorization.

## Finite-spectrum extension
For any distinct finite \(s_1,\dots,s_m\), there exists a finite prime subset giving full rank. If a nontrivial exponential polynomial vanished at every log prime, normalize by the smallest real part. Its leading trigonometric polynomial would tend to zero along log primes. The PNT gives gaps \(\log p_{n+1}-\log p_n\to0\), and bounded derivatives extend this convergence to all large real arguments. That contradicts the positive mean-square \(\sum |c_j|^2\) of a nonzero trigonometric polynomial.

## Quantitative obstruction
Set \(L_1=\log(3/2)\), \(\alpha=\log(5/2)/L_1\), \(b_k=2\pi k/L_1\) at convergent denominators of irrational \(\alpha\). For the pair \(s_0,s_0+i b_k\),
\[
 \lambda_{\min}(G_{\{2,3,5\}}) \le
  2\pi^2 5^{-2\Re s_0}\|k\alpha\|^2 < 2\pi^2 5^{-2\Re s_0}/k^2.
\]
Thus exact separation does not yield a uniform frame bound over unbounded frequencies.

## Compatibility obstruction
For a diagonal algebraic Casimir model \(A_{ii}=s_i(1-s_i)=\lambda_i\), the symmetry defect is
\[
 (A^*G_P-G_P A)_{ij}=(\bar\lambda_i-\lambda_j)(G_P)_{ij}.
\]
Even if all \(\lambda_i\) are real, off-diagonal elements typically remain nonzero. **Prime-sampled positivity is not an arithmetic invariant metric**, and no identification with Zagier's representation is claimed.

**Next formal step:** a kernel-checked elementary frame lemma assuming a nonzero observation minor, followed separately by formalization of the irrational logarithm ratio. This must not be encoded as an unproved axiom.


## Iteration 15: prime phase and frequency-window results

A standalone three-page research addendum develops four results:

1. **Real 2-prime theorem:** distinct real exponents have linearly independent evaluation columns at primes 2 and 3. Both its nonzero-minor and positive-energy statements are Lean-checked in Actions run [44879252060](https://github.com/aakbarie/SPTB-RiemannHypothesis/actions/runs/44879252060).
2. **Arithmetic nonaliasing:** for all naturals m and n>0, 5^n 2^m != 3^m 2^n. This number-theoretic result is Lean-checked in [44879158825](https://github.com/aakbarie/SPTB-RiemannHypothesis/actions/runs/44879158825).
3. **Complex 3-prime theorem:** use the irrational gap ratio log(5/2)/log(3/2), so exact phase aliases cannot coincide across both prime ratios. Its formalization is underway; do not equate the Nat-power lemma with a formal proof of the full complex result.
4. **Bounded-frequency quantitative frame:** with a fixed base exponent and separation ε<=|δ|<=B, the minimum eigenvalue of the 3-prime Gram is bounded below by c(ε,B)>0 by continuity and compactness. This does not contradict near-aliasing for |δ|→∞.

### Exact trace-energy diagnostic and its limitation

For a finite diagonal model A=diag(λ_j), λ_j=s_j(1-s_j), G=O_P* O_P, and U(t)=exp(-itA), set Q_T=(2T)^(-1) integral_{-T}^T tr(U(t)* G U(t)) dt. One obtains

  Q_T = sum_j G_jj sinh(2 Im(λ_j) T)/(2 Im(λ_j) T),

with the continuous value 1 when Im(λ_j)=0. Since all G_jj>0, Q_T stays bounded iff each Im(λ_j)=0. This is an **equivalence in the finite model**, not an independent arithmetic conservation law. It cannot prove RH without proving a uniform bound from separate arithmetic premises.

**Open bridge:** identify the zeta-associated operator and a genuine arithmetic conservation/trace identity, then prove a bound independent of the location of its zeros. Mere positivity of prime evaluation Grams cannot accomplish this.
