"""Stage CM 2: how much of the (1-e)W term the local route needs.  In the local form the term is (y/2) E_rho[(s^2/tau) R'^2]-weighted instead of
(y/2) E_rho[R'^2]; the factor kappa = E[(s^2/tau) R'^2]/E[R'^2] can be << 1 (a bright ripple at the focus centre).  Here the x-free local function is
minimised with the W-term G sigma^2/(16 b) multiplied by 1, 1/2 and 0, with the proved local brightness bound Psi(mu, b).  Python 3.11.
-> results/cm_wterm.json
"""
import json, math
import numpy as np, cm_local_bound as L
QB=0.25; LD=np.longdouble
def env_local(ell, wfac, n=40, rounds=8):
    sig = np.exp(LD(-ell)); sell = np.sqrt(LD(ell))
    lo = np.array([math.log(1e-3), -2*ell + math.log(1e4), -3*ell]); hi = np.array([math.log(4*ell+40), 0.0, math.log(0.25)])
    best=None; cen=None; span=hi-lo
    for k in range(rounds):
        if cen is None: a,c = lo,hi
        else: a=np.maximum(lo,cen-span/2); c=np.minimum(hi,cen+span/2)
        g1=[np.linspace(a[i],c[i],n) for i in range(3)]
        LM,LB,LR=np.meshgrid(*g1,indexing='ij'); mu=np.exp(LM.astype(LD)); B=np.exp(LB.astype(LD)); rr=np.exp(LR.astype(LD))
        lP=np.array([[min(math.log(2*m*QB*QB), L.Psi(m,bb)) for bb in np.exp(g1[1])] for m in np.exp(g1[0])])
        P=np.exp(np.asarray(lP,dtype=LD))[:,:,None]; G=P*np.sqrt(rr/(2+rr))
        kap=2*QB*np.sqrt(rr/(1+rr)); cb=1.184
        alpha=rr*2*(1-kap)+G/(mu*(2+rr)); beta=B/(8*rr)+G/(16*(2+rr))
        num=2*np.sqrt(alpha*beta)+3*sig/8+wfac*G*sig/(16*B)-0.023*sig
        val=num/(cb*B+G/2)*sell
        j=np.unravel_index(int(np.argmin(val)),val.shape); v=float(val[j]); cen=np.array([LM[j],LB[j],LR[j]],dtype=float)
        best=v if best is None else min(best,v); span=(c-a)/4
    return best
res={}
for ell in (55.26204223, 110.52408446, 221.04816893, 960.8):
    T10=1+12.1/ell+2.5e-4
    row={'T10':T10,'x1':env_local(ell,1.0),'x0.5':env_local(ell,0.5),'x0':env_local(ell,0.0)}
    res[str(round(ell,2))]=row
    print(round(ell,2),'T10 %.4f'%T10, 'W-term x1: %.4f  x0.5: %.4f  x0: %.4f'%(row['x1'],row['x0.5'],row['x0']), flush=True)
json.dump(res, open('results/cm_wterm.json','w'), indent=1)
