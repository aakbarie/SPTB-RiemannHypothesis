"""Exact local algebra and finite geometry checks for the proposed mesh proof.
These checks supplement, and do not replace, the written analytic argument.
"""
from fractions import Fraction as F
from math import log, floor, ceil

def trim(p):
    while len(p)>1 and p[-1]==0: p.pop()
    return p
def add(a,b):
    c=[F(0)]*max(len(a),len(b))
    for i,x in enumerate(a):c[i]+=x
    for i,x in enumerate(b):c[i]+=x
    return trim(c)
def scale(a,k):return trim([x*k for x in a])
def sub(a,b):return add(a,scale(b,-1))
def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim(c)
def power(a,n):
    r=[F(1)]
    for _ in range(n):r=mul(r,a)
    return r
def proj(a,b):
    return add(sub(scale(mul(a[0],b[0]),4),scale(add(mul(a[0],b[1]),mul(a[1],b[0])),6)),scale(mul(a[1],b[1]),12))

one=[F(1)]; q=[F(0),F(1)]; ell=sub(one,q)
vd=[ell,scale(sub(one,power(q,2)),F(1,2))]
ve=[scale(power(ell,2),F(1,2)),add(sub([F(1,3)],scale(q,F(1,2))),scale(power(q,3),F(1,6)))]
m11=sub(ell,proj(vd,vd))
m12=sub(scale(power(ell,2),F(1,2)),proj(vd,ve))
m22=sub(scale(power(ell,3),F(1,3)),proj(ve,ve))
assert m11==mul(mul(q,ell),sub(one,scale(mul(q,ell),3)))
assert m12==scale(mul(mul(power(q,2),power(ell,2)),sub(scale(q,2),one)),F(1,2))
assert m22==scale(mul(power(q,3),power(ell,3)),F(1,3))
assert sub(mul(m11,m22),power(m12,2))==scale(mul(power(q,4),power(ell,4)),F(1,12))
print('PASS: local projection matrix and determinant, exact rational polynomial identities')

def ev(p,x):return sum(c*x**i for i,c in enumerate(p))
for qv in [F(1,8),F(1,4),F(1,2),F(3,4),F(7,8)]:
    assert ev(m11,qv)>0 and ev(m11,qv)*ev(m22,qv)-ev(m12,qv)**2>0
# Centered step: invisible on a partition with a knot at 1/2;
# observed on the entire [0,1] block, its optimal affine residual is 1/16.
assert ev(m11,F(1,2))==F(1,16)
print('PASS: an aligned step has zero first-mesh error but exact second-mesh error 1/16')

for kappa in [1.0,6.283185307179586,20.0]:
    for T in [100.0,1000.0,10000.0]:
        h=kappa/log(T); hp=kappa/log(T+kappa/2)
        phases=[]
        for j in range(ceil(T/(2*h)),floor(3*T/(4*h))+1):
            x=j*h
            phase=x/kappa*log(1+kappa/(2*T))
            assert 1/8-1e-12<=phase<=3/8+1e-12
            assert abs((x/hp-j)-phase)<1e-9
            phases.append(phase)
        assert h/2<=hp<h
print('PASS: complementary mesh placement in nine finite configurations')
print('All checks passed; no zeta-zero simulation or RH verification was performed.')
