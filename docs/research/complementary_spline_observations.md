# Complementary spline observations and the full-field Horocycle converse

**Research draft prepared for Akbar Esfahani, 10 September 2026.**

**Status:** A proposed, self-contained proof of a precisely stated converse for the fixed full harmonic field. The argument has not been independently reviewed. Exact local algebra and finite mesh-placement checks accompany the analytic reasoning. This is not a proof of the Riemann Hypothesis and does not establish an unconditional SPTB upper bound.

## 1. The statement being addressed

Fix positive constants \(\kappa\) and \(\alpha\ge1\), and \(\sigma\in[1/2,1)\). Let the positive-ordinate nontrivial zeta zeros, counted with multiplicity, be

\[
\rho=\beta+i\gamma,\qquad \gamma>0,
\qquad \eta_\rho=\beta-\sigma,
\qquad a_\rho=|\rho|^{-\alpha}>0.
\]

Use the fixed full field

\[
H_\sigma(t)=\lim_{X\to\infty}
\sum_{0<\gamma\le X}a_\rho e^{\eta_\rho t}\cos(\gamma t)
\quad\text{in }L^2_{\mathrm{loc}}([0,\infty)). \tag{1}
\]

Existence of this limit is justified in Section 5. Summing over both conjugates, as in the manuscript, multiplies this field by two and does not change any conclusion.

For every sufficiently large **real** \(T\), set

\[
h_T=\frac{\kappa}{\log T}.
\]

Partition \([0,T]\) at the multiples of \(h_T\), with the final interval shortened if necessary. Let \(P_T f\) be the independent \(L^2\)-best affine fit on each block; no continuity between fits is imposed. Define the amplitude residual

\[
E_f(T)=\|f-P_Tf\|_{L^2(0,T)}^2. \tag{2}
\]

**Proposed theorem.** If, for some constants \(C,p\ge0\),

\[
E_{H_\sigma}(T)\le C(1+T)^p
\quad\text{for every sufficiently large real }T, \tag{3}
\]

then every nontrivial zeta zero satisfies \(\beta\le\sigma\).

In particular, the conjectured implication follows for any well-defined SPTB functional whose amplitude term is (2), whose other terms are nonnegative, and for which

\[
F_\lambda(H_\sigma;T,h_T)\ll T\log T\log\log T
\quad\text{uniformly for all sufficiently large real }T. \tag{4}
\]

Indeed, (4) implies (3), for example with \(p=2\). The derivative penalty is not used in this proposed converse.

**Scope matters.** Both observations below must act on the same field. This proof does not automatically apply if the amplitude itself is replaced by \(H_\sigma^{(X(T))}\), if the hypothesis is known only on a sparse sequence of \(T\)'s, or if \(h_T\) is an arbitrary irregular choice merely comparable to \(1/\log T\). It does apply to the exact canonical mesh above. Truncating only a nonnegative derivative penalty does not affect the argument.

## 2. The complementary observation

The FFT butterfly identity motivating the construction is

\[
|a+\omega b|^2+|a-\omega b|^2
=2(|a|^2+|b|^2),\qquad |\omega|=1.
\]

One output can vanish while the other retains the signal. For splines, the analogous issue is that a jump in value or slope can hide at a block boundary. A second partition with that boundary inside one of its blocks detects the jump.

This is an analogy between two concrete linear observations. FFT unitarity does not itself prove the estimates that follow.

### 2.1 Exact local calculation

On \([0,1]\), consider a piecewise-affine function with one break at \(q\in(0,1)\). After subtracting its left affine part it has the form

\[
g(t)=d\,\mathbf1_{t>q}+e(t-q)_+.
\]

Here \(d\) is the value jump and \(e\) the slope jump. Projecting onto \(\operatorname{span}\{1,t\}\) gives

\[
\inf_{\ell\ \mathrm{affine}}\|g-\ell\|_2^2
=\begin{pmatrix}\bar d&\bar e\end{pmatrix}
M(q)\begin{pmatrix}d\\e\end{pmatrix}, \tag{5}
\]

where, writing \(r=1-q\),

\[
M(q)=
\begin{pmatrix}
qr(1-3qr)&\frac12q^2r^2(2q-1)\\
\frac12q^2r^2(2q-1)&\frac13q^3r^3
\end{pmatrix}.
\]

To verify this, the affine Gram matrix and its inverse are

\[
G=\begin{pmatrix}1&1/2\\1/2&1/3\end{pmatrix},
\qquad G^{-1}=\begin{pmatrix}4&-6\\-6&12\end{pmatrix}.
\]

The moments of the step and hinge are respectively

\[
v_d=\left(r,\frac{1-q^2}{2}\right),
\qquad
v_e=\left(\frac{r^2}{2},\frac13-\frac q2+\frac{q^3}{6}\right).
\]

Subtracting \(v_i^*G^{-1}v_j\) from their raw inner products yields (5). Crucially,

\[
\det M(q)=\frac{q^4(1-q)^4}{12}>0. \tag{6}
\]

For \(q\in[1/8,7/8]\), the trace is less than one, so

\[
\lambda_{\min}(M(q))\ge c_*,
\qquad c_*:=\frac{(7/64)^4}{12}>0. \tag{7}
\]

Scaling to a block of length \(h'\) gives the explicit local bound

\[
\inf_{\ell\ \mathrm{affine}}
\int_{J}|g-\ell|^2
\ge c_*\bigl(h'|d|^2+(h')^3|e|^2\bigr), \tag{8}
\]

provided the break lies between one-eighth and seven-eighths of the way across the block. Thus an interior jump cannot be concealed by one affine fit.

### 2.2 The second mesh is already available in the hypothesis

Set

\[
T_+=T+\frac\kappa2,\qquad
h=h_T,\qquad h'=h_{T_+}<h,
\qquad J_T=[T/2,3T/4].
\]

At a first-mesh knot \(x_j=jh\in J_T\), its coordinate in the second mesh is

\[
\frac{x_j}{h'}
=j+\theta_j,
\qquad
\theta_j=\frac{x_j}{\kappa}
\log\left(1+\frac{\kappa}{2T}\right). \tag{9}
\]

For \(T\ge\kappa/2\), use \(u/2\le\log(1+u)\le u\) for \(0\le u\le1\). Since \(x_j/T\in[1/2,3/4]\),

\[
\frac18\le\theta_j\le\frac38. \tag{10}
\]

Thus each such knot lies uniformly inside a second-mesh block. For large \(T\), \(h/2\le h'<h\). Since the second blocks are shorter than the spacing of the first knots, different first knots occupy different second blocks. All those crossing blocks are contained in \([0,T]\) for sufficiently large \(T\).

This uses the exact mesh formula and the availability of the bound at both \(T\) and \(T_+\).

## 3. Recovering the discarded affine pieces

Let \(P=P_Tf\), and let \(Q=P_{T_+}f\) restricted to \([0,T]\). Write

\[
\mathcal E_T=E_f(T)+E_f(T_+).
\]

The triangle inequality gives

\[
\|P-Q\|_{L^2(0,T)}^2\le2\mathcal E_T. \tag{11}
\]

On a second-mesh block crossing a first-mesh knot, \(P\) has exactly two affine pieces and \(Q\) is affine. Apply (8), then sum over those disjoint crossing blocks. If \(d_j,e_j\) are the value and slope jumps of \(P\) at knots strictly inside \(J_T\), then

\[
\sum_j\bigl(h|d_j|^2+h^3|e_j|^2\bigr)
\le C_1\mathcal E_T, \tag{12}
\]

where \(C_1\) is independent of \(T\) and \(f\).

Let \(\ell_T\) be the affine expression of \(P\) on the first first-mesh piece meeting \(J_T\). At almost every point of \(J_T\),

\[
P(t)-\ell_T(t)
=\sum_{x_j<t}\{d_j+e_j(t-x_j)\}. \tag{13}
\]

There are at most \(C T/h\) such knots. Cauchy-Schwarz, \(|t-x_j|\le T\), and (12) therefore give

\[
\begin{aligned}
\|P-\ell_T\|_{L^2(J_T)}^2
&\le C T\frac{T}{h}
\left(\sum_j|d_j|^2+T^2\sum_j|e_j|^2\right)\\
&\le C\left(\frac{T^2}{h^2}+\frac{T^4}{h^4}\right)\mathcal E_T.
\end{aligned}
\]

Adding \(f-P\) proves the complementary-observation estimate

\[
\boxed{
\inf_{\ell\ \mathrm{affine}}
\|f-\ell\|_{L^2(J_T)}^2
\le C_2(1+T/h_T)^4
\{E_f(T)+E_f(T+\kappa/2)\}.
} \tag{14}
\]

The loss is polynomial. No frequency spacing, spectral sign, or zeta property has entered this estimate. The unavoidable common nullspace consists of global affine functions on the interval.

## 4. From local residual bounds to polynomial growth of the full field

Suppose \(f\in L^2_{\mathrm{loc}}([0,\infty))\) satisfies (3). Since
\((1+T/h_T)^4\ll T^4(\log T)^4\), (14) yields, for \(r=p+5>3\),

\[
\inf_{\ell\ \mathrm{affine}}
\|f-\ell\|_{L^2(J_T)}^2\ll T^r. \tag{15}
\]

Fix a sufficiently large \(T_0\), put \(T_n=(4/3)^nT_0\), and choose the best affine fit
\(\ell_n(t)=a_n+b_nt\) on \(J_{T_n}\). Consecutive intervals overlap on

\[
J_{T_n}\cap J_{T_{n+1}}
=[2T_n/3,3T_n/4],
\]

an interval of length \(T_n/12\). On this overlap, (15) implies

\[
\|\ell_{n+1}-\ell_n\|_2\ll T_n^{r/2}. \tag{16}
\]

For any affine function \(a+bt\) on a fixed-proportion interval at distance comparable to \(T_n\) from zero, a rescaling \(t=T_nu\) and the positive-definite \(2\times2\) Gram matrix give

\[
|a|\le C T_n^{-1/2}\|a+bt\|_2,
\qquad
|b|\le C T_n^{-3/2}\|a+bt\|_2. \tag{17}
\]

Apply this to (16):

\[
|a_{n+1}-a_n|\ll T_n^{(r-1)/2},
\qquad
|b_{n+1}-b_n|\ll T_n^{(r-3)/2}.
\]

Both exponents are positive. Summing geometric series, including the fixed initial coefficients, gives

\[
|a_n|\ll1+T_n^{(r-1)/2},
\qquad
|b_n|\ll1+T_n^{(r-3)/2}.
\]

It follows that \(\|\ell_n\|_{L^2(J_{T_n})}^2\ll T_n^r\), and hence
\(\|f\|_{L^2(J_{T_n})}^2\ll T_n^r\). The overlapping intervals cover \([T_0/2,\infty)\); summing their geometric bounds and including the initial compact interval proves

\[
\boxed{\int_0^R|f(t)|^2\,dt\ll1+R^r.} \tag{18}
\]

In particular, its primitive \(K(t)=\int_0^t f(u)\,du\) obeys

\[
|K(t)|\le\sqrt{t}\,\|f\|_{L^2(0,t)}
\ll1+t^{(r+1)/2}. \tag{19}
\]

Thus the uncontrolled affine pieces cannot accumulate into an exponentially growing full field when both complementary observations remain polynomially bounded.

## 5. Defining the full zero field at alpha = 1

At \(\alpha=1\), absolute convergence of (1) is not available. Square-summability of coefficients alone is also insufficient for arbitrary clustered frequencies. The following argument supplies the missing local convergence using the zeta zero count per unit interval.

The standard Riemann-von Mangoldt count with \(O(\log V)\) remainder implies

\[
N(V+1)-N(V)=O(\log(V+2)). \tag{20}
\]

Include both signs of each ordinate when expanding the cosines. Let \(A_k\) be the sum of the absolute coefficients in the frequency bin \([k,k+1)\). Then, with finitely many low bins absorbed into the constant,

\[
A_k\ll \frac{\log(2+|k|)}{(1+|k|)^\alpha},
\qquad \sum_k A_k^2<\infty\quad(\alpha\ge1). \tag{21}
\]

Fix a compact time interval and a smooth compactly supported function \(\chi\) equal to one on it. The parameters \(\eta_\rho\) lie in a bounded interval. Integration by parts twice therefore gives the uniform bound

\[
\left|\int \chi(t)^2 e^{(\eta_\rho+\eta_{\rho'})t}
e^{i(\gamma-\gamma')t}\,dt\right|
\le C_\chi(1+|\gamma-\gamma'|)^{-2}. \tag{22}
\]

For any finite tail of the series, grouping the squared norm by unit bins gives

\[
\|\chi H_{\mathrm{tail}}\|_2^2
\le C_\chi\sum_{k,l}A_k^{\mathrm{tail}}A_l^{\mathrm{tail}}
(1+|k-l|)^{-2}
\le C_\chi'\sum_k(A_k^{\mathrm{tail}})^2. \tag{23}
\]

The last inequality follows by Cauchy-Schwarz for each difference \(k-l\), then summing the summable kernel. By (21), the right side tends to zero as the lower tail cutoff tends to infinity. Thus the symmetric zero truncations are Cauchy in \(L^2\) on every compact interval, proving (1).

Multiplicity is included in (20). No positive minimum gap between ordinates is assumed.

## 6. Pole preservation through a primitive

Set \(z_\rho=\eta_\rho+i\gamma\). Since \(\gamma>0\), \(z_\rho\ne0\). The primitive of each finite truncation is

\[
K_X(t)=\frac12\sum_{0<\gamma\le X}a_\rho
\left\{\frac{e^{z_\rho t}-1}{z_\rho}
+\frac{e^{\bar z_\rho t}-1}{\bar z_\rho}\right\}.
\]

Here

\[
\sum_{\gamma>0}\frac{a_\rho}{|z_\rho|}<\infty \tag{24}
\]

for \(\alpha\ge1\), by the zero-counting estimate. The primitive series consequently converges absolutely and uniformly on compact time intervals. Local \(L^2\) convergence of \(H_X\) identifies its limit with

\[
K(t)=\int_0^t H_\sigma(u)\,du. \tag{25}
\]

For \(\Re w\) sufficiently large, absolute integration of the primitive series yields

\[
\mathcal LK(w)
=\frac1{2w}\sum_{\gamma>0}a_\rho
\left\{\frac1{w-z_\rho}+\frac1{w-\bar z_\rho}\right\}. \tag{26}
\]

For example, termwise integration is justified by bounding each primitive term by a constant times
\((a_\rho/|z_\rho|)(e^{(1-\sigma)t}+1)\) and using (24).

Away from \(w=0\) and the displayed poles, the series in (26) converges locally uniformly: for large \(\gamma\), denominators are comparable to \(\gamma\), uniformly on each compact set, and (24) applies. The poles form a locally finite set because zeta has only finitely many zeros in a bounded region.

Fix an off-line-right zero \(\rho_0\), so \(\Re z_{\rho_0}>0\). At \(w=z_{\rho_0}\), the residue in (26) is

\[
\operatorname{Res}_{w=z_{\rho_0}}\mathcal LK(w)
=\frac1{2z_{\rho_0}}
\sum_{\rho:z_\rho=z_{\rho_0}}a_\rho\ne0. \tag{27}
\]

Other distinct zeros contribute an analytic function in a small neighborhood of this point. Repeated zeros add positive coefficients, so their common residue does not cancel. Notice that no rightmost zero or attained spectral edge was needed.

But (3), Sections 3-4, and (25) imply polynomial growth of \(K\). Its Laplace transform is therefore holomorphic on the whole half-plane \(\Re w>0\). Identity with (26) in a far-right half-plane extends throughout the connected right half-plane with its locally finite poles removed. Holomorphy of the actual transform would force the pole in (27) to be removable, contradicting its nonzero residue.

Therefore there is no zero with \(\beta>\sigma\). This completes the proposed argument for the theorem in Section 1.

## 7. What the argument would establish, and what it does not

For the specified full field and exact canonical mesh, the proposed implication is

\[
\text{polynomial SPTB bound for every large real }T
\Longrightarrow
\text{polynomial full-field }L^2\text{ growth}
\Longrightarrow
\text{no positive-real-part Laplace poles}
\Longrightarrow
\beta\le\sigma.
\]

The new candidate mechanism in this draft is the estimate (14), followed by the overlap reconstruction. It replaces a uniform lower bound for the Gram matrix of all zero modes with a comparison of two physical observation meshes. The primitive argument avoids assuming absolute convergence of the original alpha = 1 field.

At \(\sigma=1/2\), zero symmetry would turn the conclusion into RH **if the polynomial bound were independently established for the actual field**. This draft does not establish that premise. The variance-regime assumption \(\beta\le\sigma\) cannot be used to establish it in an RH proof.

The draft also does not establish:

- a sharp exponential lower bound at every sufficiently large \(T\);
- recovery of individual zero coefficients with a uniform condition number;
- equivalence between the surrogate and an unsmoothed prime-counting remainder;
- the result for a moving spectral cutoff in the amplitude term;
- any general claim that flat phase fields, affordability, or complex-plane geometry alone forbid cancellation;
- novelty relative to the approximation-theory or analytic-number-theory literature.

The existing detection estimates and the attained-edge sketch in the manuscript are not used in this proof.

## 8. Checks and review targets

The local projection formula and determinant were independently recomputed from the moment integrals using exact rational polynomial arithmetic. The reproducible checks are in [check_complementary_mesh.py](check_complementary_mesh.py). The following checks passed:

1. All three entries of \(M(q)\) agree identically as polynomials with the direct Gram-projection calculation.
2. The determinant identity (6) agrees identically as a rational polynomial.
3. A unit step at the center of a block has exact affine residual \(1/16\), while it is exactly representable on a mesh split at that point.
4. Nine finite configurations, using three values of \(\kappa\) and three values of \(T\), satisfy the mesh-placement inequalities and the coordinate identity (9).

These finite checks do not validate the complete theorem. The proof above supplies the general arguments. Before treating this draft as a settled result, independent review should particularly examine:

- the quantifier over all large real \(T\), the fixed full field, and the manuscript's actual mesh convention;
- the use of disjoint crossing blocks in (12) and the propagation of affine coefficients in Section 4;
- the local convergence construction in Section 5 and the analytic-continuation argument in Section 6;
- any intended replacement of the full amplitude by a \(T\)-dependent truncated amplitude.

Adding this research note does not revise the manuscript or compiled PDFs. No R simulation was rerun. No computation of zeta zeros was used as proof evidence.

## 9. Source alignment

The target is the Horocycle Conjecture and harmonic-field definition in Akbar Esfahani's SPTB repository, inspected at commit `f8dbb8802cd6dc605e3f06d2286ba5c9a4ae42d7`:

- [Harmonic field, canonical regime, and truncation remark](https://github.com/aakbarie/SPTB-RiemannHypothesis/blob/f8dbb8802cd6dc605e3f06d2286ba5c9a4ae42d7/paper_clean/parts/part1.tex).
- [Horocycle conjecture and rigidity discussion](https://github.com/aakbarie/SPTB-RiemannHypothesis/blob/f8dbb8802cd6dc605e3f06d2286ba5c9a4ae42d7/paper_clean/parts/part2.tex).

The complementary-output motivation comes from page 6, Sections 4.2-4.3, of the uploaded *A Geometric Semantics of the Fast Fourier Transform: Schedule Invariance and Structure*, dated 12 January 2026. Only its displayed butterfly matrix is needed here.

Standard analytic inputs are the Riemann-von Mangoldt zero count with its classical logarithmic remainder, elementary finite-dimensional projection algebra, and Laplace-transform analyticity. [DLMF Section 25.10](https://dlmf.nist.gov/25.10) records the critical strip, zero symmetry, and zero-counting context. [DLMF Section 1.14(iii)](https://dlmf.nist.gov/1.14#iii) records the Laplace-transform framework. These references support the standard background, not the proposed complementary-mesh theorem.
