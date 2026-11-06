"""Stage CO: a second, independent minimiser of the local envelope with the proved brightness bound Psi_2 (the grid of co_envpsi2 found a lower
value on a COARSER grid at ell = 55.26, so the refined grid is not trusted on its own).  Multi-start Nelder-Mead in (log mu, log b, log r, log kappa),
starts = the grid argmins of co_envpsi2/cn_envloc plus random points, long double evaluation of the same function.  Python 3.11.
-> results/co_envmin.json
"""
import json, math, sys, time
import numpy as np
from scipy.optimize import minimize
import cn_lemmas as CL
QB = 0.25; LD = np.longdouble


def Psi2_opt(mu, b, kap):
    """Lemma CN-1's corollary with (K, L) optimised continuously (any real K > 1, L > 0 with L K kappa <= 1 is admissible; cn_lemmas.Psi2 used the
    grids K in {2..64}, L in {1..16}): coarse grid, then three local refinements in (log K, log L).  Returns P (not capped)."""
    import cn_search as S
    best = (S.psi_raw(mu, b), None)
    Ks = np.exp(np.linspace(math.log(1.2), math.log(1e6), 22)); Ls = np.exp(np.linspace(math.log(0.05), math.log(1e4), 22))
    for it in range(3):
        for K in Ks:
            for L in Ls:
                if L * K * kap > 1: continue
                lo, hi = -700.0, 3.0
                if CL.rhs2(math.exp(hi), mu, b, kap, K, L) >= mu: continue
                for _ in range(60):
                    m = 0.5 * (lo + hi)
                    if CL.rhs2(math.exp(m), mu, b, kap, K, L) >= mu: lo = m
                    else: hi = m
                if math.exp(lo) < best[0]: best = (math.exp(lo), (K, L))
        if best[1] is None: return best[0]
        K0, L0 = best[1]; w = 1.2 / (it + 1)
        Ks = np.exp(np.linspace(math.log(K0) - w, math.log(K0) + w, 7)); Ls = np.exp(np.linspace(math.log(L0) - w, math.log(L0) + w, 7))
        Ks = Ks[Ks > 1.0]; Ls = Ls[Ls > 0]
    return best[0]


def f(z, ell, kmin=1e-30, opt=False):
    lm, lb, lr, lk = z
    if not (math.log(1e-3) <= lm <= math.log(4 * ell + 40) and -2 * ell + math.log(1e4) <= lb <= 0 and -3 * ell <= lr <= math.log(0.25) and math.log(kmin) <= lk <= math.log(3)): return 10.0
    mu, bb, kk = math.exp(lm), math.exp(lb) if lb > -700 else 0.0, math.exp(lk)
    P = LD(min(2 * mu * QB * QB, Psi2_opt(mu, bb, kk) if opt else CL.Psi2(mu, bb, kk)))
    sig = np.exp(LD(-ell)); B = np.exp(LD(lb)); rr = np.exp(LD(lr)); kap = LD(kk); mu = LD(mu)
    G = P * np.sqrt(rr / (2 + rr)); kL = 2 * QB * np.sqrt(rr / (1 + rr)); cb = LD(1.184)
    alpha = rr * 2 * (1 - kL) + G / (mu * (2 + rr)); beta = B / (8 * rr) + G / (16 * (2 + rr))
    num = 2 * np.sqrt(alpha * beta) + 3 * sig / 8 + kap * G * sig / (16 * B) - LD(0.023) * sig
    return float(num / (cb * B + G / 2) * np.sqrt(LD(ell)))


def run(ell, key, nstart=40, seed=1, opt=False):
    rng = np.random.default_rng(seed)
    starts = []
    for fn, k in (('results/co_envpsi2.json', 'argmin_Psi2'), ('results/cn_envloc.json', 'argmin_C=50'), ('results/cn_envloc.json', 'argmin_C=inf'), ('results/cn_envloc.json', 'argmin_C=1.41421')):
        try:
            a = json.load(open(fn))[key][k]; starts.append([math.log(a['mu']), a['log_b_over_sigma2'] - 2 * ell, a['log_r'], math.log(max(a["kappa"], 1e-30))])
        except Exception: pass
    lo = np.array([math.log(1e-3), -2 * ell + math.log(1e4), -3 * ell, math.log(1e-30)]); hi = np.array([math.log(4 * ell + 40), 0.0, math.log(0.25), math.log(3.0)])
    for _ in range(nstart): starts.append(list(lo + rng.uniform(0, 1, 4) * (hi - lo)))
    best = (math.inf, None); hist = []
    for s in starts:
        r = minimize(f, np.array(s), args=(ell, 1e-30, opt), method='Nelder-Mead', options={'xatol': 1e-4, 'fatol': 1e-7, 'maxfev': 2500 if not opt else 600, 'initial_simplex': None})
        hist.append(float(r.fun))
        if r.fun < best[0]: best = (float(r.fun), [float(t) for t in r.x])
    z = best[1]
    return {'min': best[0], 'argmin': {'mu': math.exp(z[0]), 'log_b_over_sigma2': z[1] + 2 * ell, 'log_r': z[2], 'kappa': math.exp(z[3])}, 'starts': len(starts), 'local_minima_sorted': sorted(hist)[:8]}


if __name__ == '__main__':
    ells = {'55.26': 55.26204223, '110.52': 110.52408446, '221.05': 221.04816893, '455.26': 455.26, '690.78': 690.77552790, '960.8': 960.8}
    opt = 'opt' in sys.argv; fn = 'results/co_envmin.json' if not opt else 'results/co_envmin_opt.json'
    try: res = json.load(open(fn))
    except FileNotFoundError: res = {}
    for key, ell in ells.items():
        if key in res and 'quick' not in sys.argv: continue
        t0 = time.time(); r = run(ell, key, nstart=int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 40, opt=opt); r['T0'] = 1 + 2.1 / ell + 2.5e-4; r['seconds'] = round(time.time() - t0)
        res[key] = r; print(key, 'min %.5f (T0 %.4f) argmin %s  local minima %s  %ds' % (r['min'], r['T0'], r['argmin'], ['%.4f' % v for v in r['local_minima_sorted']], r['seconds']), flush=True)
        json.dump(res, open(fn, 'w'), indent=1)
