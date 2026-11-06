"""Stage CO 2.2: the local x-free envelope (cn_envloc.env_kappa) with the PROVED joint bound Psi_2 of Lemma CN-1's corollary in place of the
conjectured min(1, C sqrt kappa), and the "no joint lemma" envelope with the kappa grid extended downwards.  Python 3.11.
  G_max = min(2 mu Qbar^2, Psi_2(mu, b, kappa)) g_s   (Psi_2 <= Psi; Psi_2 = Psi when no (K, L) is admissible),   W-term x kappa.
Targets: T_D = 1 + (2.1 + D)/ell + 2.5e-4 for D = 0, 2, 5 and T_10.
Finding (CO section 1): the C = inf envelope of CN depends on the lower edge kappa >= 1e-4 of the grid; with kappa -> 0 the W-term vanishes and the
envelope falls to the W-free value (0.81 at ell = 55.26, cm_wterm.json).  So some kappa-dependent brightness bound is needed at every D;
the proved one (Psi_2) is enough iff the column 'Psi2' exceeds T_D.
-> results/co_envpsi2.json
"""
import json, math, sys, time
import numpy as np
import cn_lemmas as CL
import cn_envloc as EK
QB = 0.25; LD = np.longdouble


def env_psi2(ell, n=16, rounds=6, kmin=1e-12, cache=None):
    cache = {} if cache is None else cache
    sig = np.exp(LD(-ell)); sell = np.sqrt(LD(ell))
    lo = np.array([math.log(1e-3), -2 * ell + math.log(1e4), -3 * ell, math.log(kmin)]); hi = np.array([math.log(4 * ell + 40), 0.0, math.log(0.25), math.log(3.0)])
    best = None; cen = None; span = hi - lo
    for k in range(rounds):
        if cen is None: a, c = lo, hi
        else: a = np.maximum(lo, cen - span / 2); c = np.minimum(hi, cen + span / 2)
        g1 = [np.linspace(a[i], c[i], n) for i in range(4)]
        LM, LB_, LR, LK = np.meshgrid(*g1, indexing='ij')
        mu = np.exp(LM.astype(LD)); B = np.exp(LB_.astype(LD)); rr = np.exp(LR.astype(LD)); kap = np.exp(LK.astype(LD))
        lP = np.empty((n, n, n))
        for i, lm in enumerate(g1[0]):
            for j, lb in enumerate(g1[1]):
                for l, lk in enumerate(g1[3]):
                    key = (round(lm, 9), round(lb, 9), round(lk, 9))
                    if key not in cache:
                        m, bb, kk = math.exp(lm), math.exp(lb), math.exp(lk)        # bb may underflow to 0.0 (b = 0 in Psi_2: harmless)
                        cache[key] = math.log(min(2 * m * QB * QB, CL.Psi2(m, bb, kk)))
                    lP[i, j, l] = cache[key]
        P = np.exp(np.asarray(lP, dtype=LD))[:, :, None, :]
        G = P * np.sqrt(rr / (2 + rr))
        kL = 2 * QB * np.sqrt(rr / (1 + rr)); cb = 1.184
        alpha = rr * 2 * (1 - kL) + G / (mu * (2 + rr)); beta = B / (8 * rr) + G / (16 * (2 + rr))
        num = 2 * np.sqrt(alpha * beta) + 3 * sig / 8 + kap * G * sig / (16 * B) - 0.023 * sig
        val = num / (cb * B + G / 2) * sell
        j = np.unravel_index(int(np.argmin(val)), val.shape)
        v = float(val[j]); cen = np.array([LM[j], LB_[j], LR[j], LK[j]], dtype=float)
        best = v if best is None else min(best, v)
        span = (c - a) / 4
    return best, {'mu': math.exp(cen[0]), 'log_b_over_sigma2': float(cen[1] + 2 * ell), 'log_r': float(cen[2]), 'kappa': math.exp(cen[3])}


if __name__ == '__main__':
    quick = 'quick' in sys.argv[1:]
    ells = (55.26204223, 110.52408446, 221.04816893, 455.26, 690.77552790, 960.8) if not quick else (55.26204223,)
    fn = 'results/co_envpsi2.json' if not quick else 'results/co_envpsi2_quick.json'
    try: res = json.load(open(fn))
    except FileNotFoundError: res = {}
    for ell in ells:
        if str(round(ell, 2)) in res: continue
        t0 = time.time(); row = {'T10': 1 + 12.1 / ell + 2.5e-4}
        for Dd in (0.0, 2.0, 5.0): row['T_D=%g' % Dd] = 1 + (2.1 + Dd) / ell + 2.5e-4
        for kmin in (1e-4, 1e-8, 1e-12):
            v, at = EK.env_kappa(ell, float('inf'), n=16, rounds=6, kmin=kmin); row['Cinf_kmin=%g' % kmin] = v; row['argmin_Cinf_kmin=%g' % kmin] = at
        v, at = env_psi2(ell, n=12 if quick else 16, rounds=5 if quick else 6); row['Psi2'] = v; row['argmin_Psi2'] = at
        res[str(round(ell, 2))] = row
        print(round(ell, 2), ' '.join('%s:%.4f' % (k, v) for k, v in row.items() if not k.startswith('argmin')), 'argmin_Psi2', at, '%.0fs' % (time.time() - t0), flush=True)
        json.dump(res, open(fn, 'w'), indent=1)
