#!/usr/bin/env python3
"""Finite-field numerical illustrations; never an RH test or full-field tail bound.
Run: python run_experiments.py   (numpy, scipy, matplotlib, mpmath required).
All output paths are relative to this script, and random numbers are not used.
"""
from pathlib import Path
import csv, json, platform, sys
import numpy as np
import scipy
from scipy.special import roots_legendre
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT = Path(__file__).resolve().parent
KAPPA = 1.0
ORDERS = (32,64,128)
HORIZONS = (8.,16.,32.,64.)

def write_csv(name, rows):
    with (OUT/name).open('w', newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def field(t, gammas, amplitudes, etas):
    return np.sum(amplitudes*np.exp(t[...,None]*etas)*np.cos(t[...,None]*gammas),axis=-1)

def observe(T, gammas, amplitudes, etas, order, kappa=KAPPA):
    h=kappa/np.log(T)
    left=np.arange(int(np.ceil(T/h)),dtype=float)*h
    right=np.minimum(left+h,T)
    keep=right>left; left=left[keep]; right=right[keep]
    half=(right-left)/2;mid=(right+left)/2
    x,w=roots_legendre(order)
    t=mid[:,None]+half[:,None]*x
    f=field(t,gammas,amplitudes,etas)
    # In centered coordinates, 1 and x have Gram diagonal (2,2/3).
    c0=np.sum(f*w,axis=1)/2
    c1=3*np.sum(f*x*w,axis=1)/2
    residual=f-c0[:,None]-c1[:,None]*x
    energy=np.sum(half*np.sum(w*f*f,axis=1))
    error=np.sum(half*np.sum(w*residual*residual,axis=1))
    return float(energy),float(error),len(left),h

def exact_energy(T,gammas,amplitudes,etas):
    # Finite exponential expansion integrated analytically; independent of quadrature.
    z=np.concatenate((etas+1j*gammas,etas-1j*gammas))
    c=np.concatenate((amplitudes/2,amplitudes/2))
    s=z[:,None]+z[None,:]
    integral=np.full(s.shape,complex(T))
    nonzero=np.abs(s)>1e-14
    integral[nonzero]=np.expm1(s[nonzero]*T)/s[nonzero]
    value=np.sum(c[:,None]*c[None,:]*integral)
    if abs(value.imag)>1e-10*max(1.,abs(value.real)):
        raise ArithmeticError('Unexpected imaginary energy')
    return float(value.real)

def main():
    mp.mp.dps=40
    zeros=[mp.zetazero(k) for k in range(1,65)]
    write_csv('zeta_zeros.csv',[dict(index=k+1,beta=mp.nstr(z.real,40),gamma=mp.nstr(z.imag,40),abs_zeta=mp.nstr(abs(mp.zeta(z)),6)) for k,z in enumerate(zeros)])
    gamma=np.array([float(z.imag) for z in zeros])
    rows=[]
    for n in (32,64):
        for alpha in (1.,2.):
            g=gamma[:n];a=np.hypot(.5,g)**(-alpha);eta=np.zeros(n)
            for T in HORIZONS:
                exact=exact_energy(T,g,a,eta)
                for q in ORDERS:
                    energy,error,blocks,h=observe(T,g,a,eta,q)
                    rows.append(dict(N=n,alpha=alpha,sigma=.5,kappa=KAPPA,T=T,quadrature_order=q,blocks=blocks,h=h,energy=energy,residual=error,energy_per_T=energy/T,residual_per_T=error/T,residual_fraction=error/energy,analytic_energy=exact,energy_relative_error=abs(energy-exact)/exact))
    write_csv('zeta_field_results.csv',rows)
    synth=[]
    # One damped/neutral/growing mode, with exactly the same frequency and amplitude.
    g=gamma[:1];a=np.hypot(.5,g)**(-1.)
    for eta_value in (-.08,0.,.08):
        for T in (*HORIZONS,128.):
            exact=exact_energy(T,g,a,np.array([eta_value]))
            for q in ORDERS:
                energy,error,blocks,h=observe(T,g,a,np.array([eta_value]),q)
                synth.append(dict(eta=eta_value,gamma=float(g[0]),amplitude=float(a[0]),T=T,quadrature_order=q,h=h,energy=energy,residual=error,energy_per_T=energy/T,residual_per_T=error/T,analytic_energy=exact,energy_relative_error=abs(energy-exact)/exact))
    write_csv('synthetic_mode_results.csv',synth)
    sensitivity=[]
    g=gamma;a=np.hypot(.5,g)**(-1.);eta=np.zeros(64)
    for kappa in (.5,1.,2.):
        for order in (64,128,256):
            energy,error,blocks,h=observe(64.,g,a,eta,order,kappa)
            sensitivity.append(dict(N=64,alpha=1.,T=64.,kappa=kappa,quadrature_order=order,h=h,blocks=blocks,energy=energy,residual=error,residual_per_T=error/64.))
    write_csv('kappa_sensitivity.csv',sensitivity)
    # Direct numerical eigenvalues supplement the exact symbolic matrix identity.
    matrix=[]
    for q in np.linspace(1/8,7/8,1001):
        r=1-q
        M=np.array([[q*r*(1-3*q*r),q*q*r*r*(2*q-1)/2],[q*q*r*r*(2*q-1)/2,q**3*r**3/3]])
        matrix.append(dict(q=q,min_eigenvalue=np.linalg.eigvalsh(M)[0],determinant=np.linalg.det(M),exact_determinant=q**4*r**4/12))
    write_csv('local_matrix_scan.csv',matrix)
    env=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__,matplotlib=matplotlib.__version__,mpmath_dps=mp.mp.dps)
    (OUT/'environment.json').write_text(json.dumps(env,indent=2)+'\n')
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.25})
    fig,axes=plt.subplots(1,2,figsize=(9,3.5))
    for alpha,ax in zip((1.,2.),axes):
        for n in (32,64):
            rr=[r for r in rows if r['N']==n and r['alpha']==alpha and r['quadrature_order']==128]
            ax.plot([r['T'] for r in rr],[r['residual_per_T'] for r in rr],'o-',label=f'N={n}')
        ax.set(xlabel='T',ylabel='Affine residual E_N(T) / T',title=f'Critical-line truncations, alpha={alpha:g}')
        ax.legend()
    fig.tight_layout()
    for ext in ('pdf','png'):fig.savefig(OUT/f'finite_zeta_residuals.{ext}',dpi=180)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(6,3.7))
    for eta in (-.08,0.,.08):
        rr=[r for r in synth if r['eta']==eta and r['quadrature_order']==128]
        ax.semilogy([r['T'] for r in rr],[r['residual_per_T'] for r in rr],'o-',label=f'eta={eta:+.2f}')
    ax.set(xlabel='T',ylabel='Affine residual / T',title='Synthetic mode: damping, neutrality, growth')
    ax.legend();fig.tight_layout()
    for ext in ('pdf','png'):fig.savefig(OUT/f'synthetic_mode_residuals.{ext}',dpi=180)
    plt.close(fig)
    def convergence(data,key,lower=64):
        groups={}
        for r in data:groups.setdefault(tuple(r[k] for k in key),{})[r['quadrature_order']]=r
        return max(abs(v[lower]['residual']-v[128]['residual'])/max(v[128]['residual'],1e-300) for v in groups.values())
    conv=convergence(rows,('N','alpha','T'));sconv=convergence(synth,('eta','T'))
    coarse=convergence(rows,('N','alpha','T'),32)
    if conv>1e-9 or sconv>1e-9:raise ArithmeticError('Residual quadrature did not converge')
    max_energy=max(r['energy_relative_error'] for r in rows+synth if r['quadrature_order']==128)
    lines=['# Numerical experiments','',
    'These are reproducible finite-field illustrations of the amplitude residual. They do not estimate the infinite-field truncation error, establish an asymptotic bound, locate unknown zeros, or prove RH. The synthetic growing mode is not claimed to be a zeta zero.','',
    '## Method','',
    'The first 64 positive critical-line zeros are computed using mpmath.zetazero at 40 decimal digits and saved in zeta_zeros.csv, with numerical |zeta(rho)| checks. The experiments use N=32,64, alpha=1,2, sigma=1/2, kappa=1, and T=8,16,32,64. These ordinates are the zeros returned by that routine; this calculation is not an independent certified zero count.','',
    'For each T the origin-anchored blocks have width 1/log(T), including a shortened last block. Gauss-Legendre quadrature independently computes the affine L2 projection on each block in the centered basis 1,x. Squared residuals are integrated directly, avoiding subtraction of two nearly equal total energies. Orders 32,64,128 are all saved. Finite-field total energy is also computed by analytically integrating the pairwise exponential expansion. The same fixed finite field is used across T for each N; no moving cutoff is used within a run.','',
    f'The largest relative difference in residual between orders 32 and 128 is {coarse:.3e} for the zeta fields, showing that the coarsest order is insufficient for the most oscillatory case. The largest relative change in residual from order 64 to 128 is {conv:.3e} for the zeta fields and {sconv:.3e} for the synthetic modes. The largest relative discrepancy between order-128 energy and the independent analytic energy is {max_energy:.3e}. These are numerical consistency checks, not rigorous quadrature error certificates.','',
    '## Critical-line finite fields','',
    '| N | alpha | T | energy / T | residual / T | residual / energy |','|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        if r['quadrature_order']==128:lines.append(f"| {r['N']} | {r['alpha']:g} | {r['T']:g} | {r['energy_per_T']:.8g} | {r['residual_per_T']:.8g} | {r['residual_fraction']:.8g} |")
    lines+=['','### Mesh constant sensitivity','',
    'At T=64, N=64 and alpha=1, all other parameters fixed, the following values use order 256. Orders 64,128,256 are retained in kappa_sensitivity.csv. Different mesh constants define different observations; comparisons are finite numerical values, not monotonicity claims for non-nested partitions.','',
    '| kappa | h | residual / T |','|---:|---:|---:|']
    for r in sensitivity:
        if r['quadrature_order']==256:lines.append(f"| {r['kappa']:g} | {r['h']:.8g} | {r['residual_per_T']:.8g} |")
    sg={}
    for r in sensitivity:sg.setdefault(r['kappa'],{})[r['quadrature_order']]=r['residual']
    sensconv=max(abs(v[128]-v[256])/v[256] for v in sg.values())
    if sensconv>1e-9:raise ArithmeticError('Sensitivity quadrature did not converge')
    lines += ['',f'Maximum relative residual difference between orders 128 and 256 in the mesh sensitivity experiment: {sensconv:.3e}.']
    lines+=['','All tested modes have real part 1/2, so neutral exponential factors are built into these finite examples. Their bounded energy density is expected; it provides no evidence about possible zeros outside the sampled set. Differences between N=32 and N=64 measure only this finite cutoff change and are not bounds on the omitted infinite tail.','',
    '## Synthetic stress test','',
    'Use one mode |1/2+i gamma_1|^{-1} exp(eta t) cos(gamma_1 t), with eta=-0.08,0,+0.08. Only the exponential rate is modified. No synthetic beta is asserted to be the real part of an actual zeta zero.','',
    '| eta | T | energy / T | residual / T |','|---:|---:|---:|---:|']
    for r in synth:
        if r['quadrature_order']==128:lines.append(f"| {r['eta']:+.2f} | {r['T']:g} | {r['energy_per_T']:.8g} | {r['residual_per_T']:.8g} |")
    lines+=['','The growing mode eventually dominates the shrinking block width in the sampled range; the experiment illustrates a mechanism, not an every-T lower bound for an infinite superposition.','',
    '## Local matrix check','',f"For 1001 equally spaced q in [1/8,7/8], the smallest observed eigenvalue of M(q) is {min(r['min_eigenvalue'] for r in matrix):.10g}. The conservative analytic lower bound (7/64)^4/12 is {(7/64)**4/12:.10g}. The sampled minimum is not a certified continuum minimum; positivity on the continuum follows from the exact determinant proof in the paper.",'',
    '## Reproduction and limitations','',
    'Run `python run_experiments.py` with numpy, scipy, mpmath, and matplotlib installed. Exact package and platform versions are recorded in environment.json. CSVs retain all quadrature orders; the figures show order 128. No randomness is used. The computations are ordinary floating point after high-precision zero generation. No interval arithmetic, full-field tail certification, asymptotic rate fit, or numerical proof is claimed.','']
    (OUT/'numerical_results.md').write_text('\n'.join(lines))
    print(json.dumps(dict(zeta_residual_relative_change_64_128=conv,synthetic_residual_relative_change_64_128=sconv,max_energy_relative_error=max_energy,files=len(list(OUT.iterdir()))),indent=2))
if __name__=='__main__':main()
