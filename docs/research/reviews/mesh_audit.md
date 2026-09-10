# Independent audit: complementary spline meshes

Scope: Sections 2–4 of `output/complementary_spline_observations.md`. This is a proof audit of the general reconstruction argument, independent of any assertion about zeta zeros. No existing draft or repository file was modified.

## Verdict

The reconstruction argument is valid. I found no fatal mathematical gap. Its assumptions must specify the mesh exactly, the common fixed function, and the observations available at complementary horizons. The full continuum of horizons is sufficient but not necessary: a geometrically spaced sequence together with its fixed-offset partners suffices. Arbitrary sparse observations do not suffice.

The proof can be made independent of the computer-checked projection formula: finite-dimensional injectivity and compactness prove the local coercivity directly. The displayed matrix and its determinant are also consistent with the moment calculation.

## Precise theorem

Fix κ>0, p≥0, and f in L²_loc([0,∞)). At each horizon t>1, partition [0,t] at multiples of h_t=κ/log t, with the last block shortened. Let P_t be the blockwise L² orthogonal projection onto affine functions and E_f(t)=||f-P_tf||²_{L²(0,t)}.

Choose a sufficiently large T₀ depending only on κ, put λ=4/3 and T_n=λⁿT₀. Suppose there is A<∞ such that

    E_f(T_n)+E_f(T_n+κ/2) ≤ A T_n^p  (n≥0).

Then there is B<∞, depending on f on an initial compact interval, A,p,κ,T₀, such that, for all R≥2,

    ∫₀ᴿ |f(t)|² dt ≤ B [1+R^(p+4)(log R)^4].

In particular, polynomial residual bounds for every sufficiently large real horizon imply polynomial full-field growth. The draft's weaker exponent p+5 is valid.

## Proof

### 1. Local coercivity

On [0,1], define s_q(x)=1_{x>q} and r_q(x)=(x-q)_+. Let Π project onto affine functions. For q in [1/8,7/8], the map

    (d,e) ↦ (I-Π)(d s_q+e r_q)

is injective. Indeed, if d s_q+e r_q equals an affine function almost everywhere, that affine function vanishes on the left interval of positive length and is therefore identically zero. On the right interval, d+e(x-q)=0 almost everywhere, forcing d=e=0. Its squared norm is a continuous positive-definite quadratic form in (d,e), continuous in q. Compactness in q and the unit sphere gives an absolute c>0 such that

    dist²(d s_q+e r_q, Aff) ≥ c(|d|²+|e|²).

This argument works for real or complex coefficients. Scaling an interval of length b transforms a physical slope jump e into coefficient be, and yields

    dist²(g, Aff; interval) ≥ c[b|d|²+b³|e|²].

For comparison with the draft's explicit computation, setting r=1-q gives

    M11=qr(1-3qr), M12=q²r²(2q-1)/2, M22=q³r³/3,
    det M=q⁴r⁴/12.

At q=1/2 these are 1/16, 0, 1/192, respectively. The draft's explicit lower constant follows from det/trace and trace<1.

### 2. Complementary geometry

Fix large T. Write h=κ/log T, h'=κ/log(T+κ/2), and I=[T/2,3T/4]. For a first-mesh knot x_j=jh strictly inside I,

    x_j/h' = j + θ_j,
    θ_j=(x_j/κ) log(1+κ/(2T)).

If T≥κ/2, then 1/8≤θ_j≤3/8, using u/2≤log(1+u)≤u for u=κ/(2T)≤1. Thus x_j is strictly inside the second block [jh',(j+1)h'] with the required relative position. Since h'<h, that block contains no other first knot. Blocks corresponding to distinct first knots are distinct and have disjoint interiors. For large T, h'≥h/2 and all these blocks are contained in [0,T]. None is the truncated terminal block. This proves the needed geometry for every sufficiently large T, without numerical mesh checks.

### 3. Sum jump bounds

Let P=P_Tf and Q=P_{T+κ/2}f restricted to [0,T]. Then

    ||P-Q||²_{L²(0,T)} ≤ 2[E_f(T)+E_f(T+κ/2)] = 2D_T.

On a crossing second block, P has one break, Q is affine, and the preceding coercivity applies. If d_j=P(x_j+)-P(x_j-) and e_j=P'(x_j+)-P'(x_j-), summing over disjoint crossing blocks gives

    Σ_j [h|d_j|²+h³|e_j|²] ≤ C D_T.

The constant is absolute once the large-T mesh comparability holds. No regularity of f beyond local L² is used: P is a finite piecewise-affine function, and all traces here belong to P, not f.

### 4. Recover one affine approximation on I

Take ℓ to be P's affine expression on the first first-mesh block having positive-length intersection with I. At almost every t in I,

    P(t)-ℓ(t)=Σ_{x_j∈int(I), x_j<t}[d_j+e_j(t-x_j)].

The positive-length convention resolves the harmless endpoint ambiguity when T/2 itself is a mesh knot. There are N≤C(1+T/h) terms. Cauchy–Schwarz, |t-x_j|≤T, and the jump bound give

    ||P-ℓ||²_{L²(I)} ≤ C T N [Σ|d_j|²+T²Σ|e_j|²]
                        ≤ C[(T/h)²+(T/h)⁴]D_T

for sufficiently large T. Adding f-P yields

    dist²(f,Aff; I) ≤ C(1+T/h)⁴D_T.

This is equation (14) of the draft. The common nullspace of the two observations on this region is the affine space, but the argument provides the necessary quantitative estimate rather than relying only on nullspace intersection.

### 5. Control affine drift

Write q=p+4≥4. At the paired horizons in the theorem, the last estimate implies

    dist²(f,Aff; I_n) ≤ C T_n^q(log T_n)^4,
    I_n=[T_n/2,3T_n/4].

Let ℓ_n(t)=a_n+b_nt be the best affine fit on I_n. The overlap

    O_n=I_n∩I_{n+1}=[2T_n/3,3T_n/4]

has length T_n/12. The triangle inequality on O_n gives

    ||ℓ_{n+1}-ℓ_n||_{L²(O_n)} ≤ C T_n^(q/2)(log T_n)².

For g(t)=a+bt, rescale t=T_nu on O_n. The Gram matrix of 1,u on [2/3,3/4] is fixed and positive definite, so

    |a|≤C T_n^(-1/2)||g||₂,
    |b|≤C T_n^(-3/2)||g||₂.

Therefore

    |a_{n+1}-a_n|≤C T_n^((q-1)/2)(log T_n)²,
    |b_{n+1}-b_n|≤C T_n^((q-3)/2)(log T_n)².

Both powers of T_n are strictly positive. Summing along the geometric sequence yields the same bounds for a_n and b_n, up to fixed initial coefficients. Consequently

    ||ℓ_n||²_{L²(I_n)}≤C T_n^q(log T_n)^4,
    ||f||²_{L²(I_n)}≤C T_n^q(log T_n)^4.

The intervals I_n overlap and cover [T₀/2,∞). For a given R, sum only through the first interval whose right endpoint reaches R; its T_n is bounded by a constant multiple of R. The resulting geometric sum is at most C R^q(log R)^4. Add the finite integral on [0,T₀/2] to finish the proof.

## Counterexamples delimiting the theorem

Let f be the continuous linear interpolant of f(m)=exp(m²) at nonnegative integers. This f is locally square integrable and has superpolynomial L² growth. For example, on [m,m+1/2] it is at least exp(m²), so its mass up to m+1 exceeds exp(2m²)/2.

1. On the exact canonical mesh only at horizons T_m=exp(κm), h_{T_m}=1/m. Every integer breakpoint is a mesh boundary and f is affine on every mesh block, hence E_f(T_m)=0. Thus an arbitrary sparse sequence is insufficient.
2. With the altered mesh h_T=1/ceil(log T/κ), which satisfies h_T~κ/log T, every integer breakpoint remains a mesh boundary at every horizon. Again E_f(T)=0 for every large T. Thus asymptotic mesh comparability alone is insufficient.

The first counterexample does not contradict the paired-horizon theorem. Adding T_m+κ/2 destroys the common alignment and exposes its large slope jumps.

## Editorial recommendations

- Replace “every large real T is essential” by “the stated all-real-T condition supplies the complementary observations; arbitrary sparse observations are insufficient.” The paired sequence gives a strictly weaker valid hypothesis.
- State explicitly that first-mesh pieces meeting I are selected by positive-length intersection, avoiding singleton endpoint ambiguity.
- State constants may depend on κ, the onset threshold, and the initial affine coefficients of f. No bound on those coefficients can follow from residual observations alone, since every global affine f has identically zero residual.
- Separate the general reconstruction theorem from any zeta application; it does not require spectral facts, positive coefficients, or a derivative penalty.
- Keep the fixed-field condition. This proof compares P_T f and P_{T+κ/2} f and does not justify substituting two different T-dependent truncated fields.
