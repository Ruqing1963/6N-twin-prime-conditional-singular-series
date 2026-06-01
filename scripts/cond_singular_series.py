#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
First-order conditional singular series  S_{2,omega}(d)   (Part III)
================================================================================
Derivation (first order, from the congruence-lockdown rule of Part II).
For a gap d between two twin centres N and N+d, condition per prime q>3 on
whether q | N:
  - q | N : left centre auto-safe; right centre N+d safe iff d != +-6^{-1} (mod q),
            equivalently q does NOT divide 6d-1 or 6d+1. So q contributes factor
            0 if q | (6d +-1) (the gap is locked out), else removes the usual
            suppression.
  - q ! N : usual unconditional factor (1 - nu_q(d)/q)/(1 - 2/q)^2,
            nu_q(d) = |{-1,1,6d-1,6d+1} mod q|.
Collecting:
      S_{2,omega}(d) ~ C2(d) * L_omega(d)
  C2(d)      = product over q>3 of (1 - nu_q(d)/q)/(1 - 2/q)^2  (unconditional)
  L_omega(d) = Pr over stratum-omega centres N that q ! N for every q | (6d +-1)
             = the "lockdown-allowance" of gap d in the stratum.
The admissibility of d is governed by the factorisation of 6d +-1:
  6*35 -+1 = 209,211 with 209 = 11*19  -> 210 locked by the common primes 11,19
  6*7  -+1 = 41,43 (both prime)        -> 42 locked only by the rare primes 41,43

This script evaluates L_omega(d) directly from the real twin centres of a shell
S_K (computing each centre's small-prime factor set), forms the normalised
prediction ~ C2(d)*L_omega(d), and compares it to the observed relative
preference r(d|omega). It prints obs vs pred vs residual and writes css_plot_S{K}.csv.

USAGE
  MAXK=9  python cond_singular_series.py     # S9  (~1-2 min)
  MAXK=10 python cond_singular_series.py     # S10 (~15 min)
Requires: numpy.

NOTE: the first-order series reproduces 30/42/60 and the direction of 210; a
systematic high-omega residual remains on 210 (observed falls below the
first-order prediction), attributed to the second-order correlation between the
factorisations of N and N+d. That residual is the Part III open problem.
================================================================================
"""
import numpy as np, math, time, os, csv
from collections import defaultdict

def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(math.isqrt(n))+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0].astype(np.int64)

MAXK=int(os.environ.get("MAXK",9))
LO=10**(MAXK-1)//6+1; HI=10**MAXK//6; SEG=4_000_000
PB=int(math.isqrt(10**MAXK))+1; BP=primes_upto(PB)
GMAX=60
SMALLQ=[5,7,11,13,17,19,23,29,31,37,41,43,47]
# for each d: which small primes lock it (divide 6d-1 or 6d+1)
LOCK={d:[q for q in SMALLQ if ((6*d-1)%q==0 or (6*d+1)%q==0)] for d in range(1,GMAX+1)}

twN=[]; twOm=[]; twMask=[]
n=LO; t0=time.time()
print(f"[setup] S{MAXK}: N in [{LO:,}, {HI:,}]")
while n<=HI:
    nh=min(n+SEG,HI+1); sz=nh-n
    rem=np.arange(n,nh,dtype=np.int64); ob=np.zeros(sz,np.int16); mask=np.zeros(sz,np.int32)
    for p in BP:
        if p*p>nh-1: break
        f=((n+p-1)//p)*p
        if f>=nh: continue
        idx=np.arange(f-n,sz,p)
        if idx.size==0: continue
        sub=rem[idx]; m=(sub%p)==0
        while m.any(): sub[m]//=p; m=(sub%p)==0
        rem[idx]=sub
        if p>3:
            ob[idx]+=1
            if p in SMALLQ: mask[idx]|=(1<<SMALLQ.index(int(p)))
    ob[rem>1]+=1
    vlo=6*n-1; vhi=6*(nh-1)+1; span=vhi-vlo+1
    comp=np.zeros(span,bool); sq=int(math.isqrt(vhi))+1
    for p in BP:
        if p>sq: break
        st=max(p*p,((vlo+p-1)//p)*p)
        if st>vhi: continue
        comp[st-vlo:span:p]=True
    Narr=np.arange(n,nh,dtype=np.int64)
    tw=(~comp[(6*Narr-1)-vlo])&(~comp[(6*Narr+1)-vlo])
    pos=np.nonzero(tw)[0]
    twN.append(Narr[pos]); twOm.append(ob[pos]); twMask.append(mask[pos])
    n=nh
twN=np.concatenate(twN); twOm=np.concatenate(twOm); twMask=np.concatenate(twMask)
print(f"[done] S{MAXK} twins {len(twN):,}; scan {time.time()-t0:.0f}s")

# observed gap distribution by omega (left-centre attribution)
g=np.diff(twN); omL=twOm[:-1]
keep=(g>=1)&(g<=GMAX); g=g[keep]; omL=omL[keep]
obs=defaultdict(lambda: np.zeros(GMAX+1))
for i in range(len(g)): obs[int(omL[i])][int(g[i])]+=1
overall=np.zeros(GMAX+1)
for om in obs: overall+=obs[om]
base_obs=overall/overall[1:GMAX+1].sum()

# C2(d)
QP=[q for q in primes_upto(200000) if q>3]
def C2(d):
    p=1.0
    for q in QP:
        nu=len({(-1)%q,1%q,(6*d-1)%q,(6*d+1)%q})
        if nu==q: return 0.0
        p*=(1.0-nu/q)/(1.0-2.0/q)**2
    return p
Cg=np.array([0.0]+[C2(d) for d in range(1,GMAX+1)])
lockbits={d:sum(1<<SMALLQ.index(q) for q in LOCK[d]) for d in range(1,GMAX+1)}

omegas=[om for om in sorted(obs) if obs[om][1:GMAX+1].sum()>=10000]
pred=defaultdict(lambda: np.zeros(GMAX+1))
for om in omegas:
    sel=twMask[twOm==om]
    for d in range(1,GMAX+1):
        lb=lockbits[d]; allow=1.0 if lb==0 else ((sel&lb)==0).mean()
        pred[om][d]=Cg[d]*allow
    pred[om][1:GMAX+1]/=pred[om][1:GMAX+1].sum()
predbase=np.zeros(GMAX+1)
for d in range(1,GMAX+1):
    lb=lockbits[d]; allow=1.0 if lb==0 else ((twMask&lb)==0).mean()
    predbase[d]=Cg[d]*allow
predbase[1:GMAX+1]/=predbase[1:GMAX+1].sum()

print(f"\nFirst-order S_2,omega(d): predicted vs observed r(d|omega), S{MAXK}")
for d in [5,7,10,35]:   # 30,42,60,210
    print(f"\n--- 6dN={6*d} (d={d}); lockdown primes of d: {LOCK[d] if LOCK[d] else 'none'} ---")
    print(f"  {'omega':>5}{'obs_r':>9}{'pred_r':>9}{'resid':>9}")
    for om in omegas:
        tot=obs[om][1:GMAX+1].sum()
        obs_r=(obs[om][d]/tot)/base_obs[d] if base_obs[d]>0 else 0
        pred_r=pred[om][d]/predbase[d] if predbase[d]>0 else 0
        print(f"  {om:>5}{obs_r:>9.3f}{pred_r:>9.3f}{obs_r-pred_r:>9.3f}")

with open(f'css_plot_S{MAXK}.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['6dN','d','omega','obs_r','pred_r','resid'])
    for d in [5,7,10,35]:
        for om in omegas:
            tot=obs[om][1:GMAX+1].sum()
            obs_r=(obs[om][d]/tot)/base_obs[d] if base_obs[d]>0 else 0
            pred_r=pred[om][d]/predbase[d] if predbase[d]>0 else 0
            w.writerow([6*d,d,om,f'{obs_r:.4f}',f'{pred_r:.4f}',f'{obs_r-pred_r:.4f}'])
print(f"\n[ok] wrote css_plot_S{MAXK}.csv")
