"""Stage CN 2.4: the local x-free envelope with kappa as a fourth variable.  Python 3.11.
Main function of Lemma R_theta (cc_core, (L1)-(L3) at their worst as in ci_envloc), with
  G <= G_max(mu, b, kappa) = min(2 mu Qbar^2, Psi(mu, b) * g(kappa)) * g_s,   g(kappa) = min(1, C sqrt(kappa))   (the joint lemma, constant C;
  the weak bound 2 mu Qbar^2 is a bound on E0 alone and carries no kappa factor),
  and the W-term G sigma^2/(16 b) multiplied by kappa (CM section 3.2: the (1-e)W term is kappa-weighted in the local form).
The minimum over (log mu, log b, log r, log kappa) is found on a grid refined 7 times (long double).  C = inf means "no joint lemma".
-> results/cn_envloc.json
"""
import json, math, sys
import numpy as np
import cm_local_bound as LB
QB = 0.25; LD = np.longdouble
from cn_search import psi_raw


def env_kappa(ell, Cj, n=20, rounds=7, kmin=1e-4, wfac=None):
    sig = np.exp(LD(-ell)); sell = np.sqrt(LD(ell))
    lo = np.array([math.log(1e-3), -2 * ell + math.log(1e4), -3 * ell, math.log(kmin)]); hi = np.array([math.log(4 * ell + 40), 0.0, math.log(0.25), math.log(3.0)])
    best = None; cen = None; span = hi - lo
    for k in range(rounds):
        if cen is None: a, c = lo, hi
        else: a = np.maximum(lo, cen - span / 2); c = np.minimum(hi, cen + span / 2)
        g1 = [np.linspace(a[i], c[i], n) for i in range(4)]
        LM, LB_, LR, LK = np.meshgrid(*g1, indexing='ij')
        mu = np.exp(LM.astype(LD)); B = np.exp(LB_.astype(LD)); rr = np.exp(LR.astype(LD)); kap = np.exp(LK.astype(LD))
        lPsi = np.array([[math.log(psi_raw(m, bb)) for bb in np.exp(g1[1])] for m in np.exp(g1[0])])
        Praw = np.exp(np.asarray(lPsi, dtype=LD))[:, :, None, None]; Wk = (2 * mu * QB * QB)
        gk = np.minimum(1.0, Cj * np.sqrt(kap)) if np.isfinite(Cj) else 1.0 + 0 * kap
        P = np.minimum(Wk, Praw * gk)                      # G_max = min(2 mu Qbar^2, Psi(mu,b) g(kappa)): the kappa factor only on the CM-1 part
        G = P * np.sqrt(rr / (2 + rr))
        kL = 2 * QB * np.sqrt(rr / (1 + rr)); cb = 1.184
        alpha = rr * 2 * (1 - kL) + G / (mu * (2 + rr))
        beta = B / (8 * rr) + G / (16 * (2 + rr))
        wk = kap if wfac is None else wfac
        num = 2 * np.sqrt(alpha * beta) + 3 * sig / 8 + wk * G * sig / (16 * B) - 0.023 * sig
        val = num / (cb * B + G / 2) * sell
        j = np.unravel_index(int(np.argmin(val)), val.shape)
        v = float(val[j]); cen = np.array([LM[j], LB_[j], LR[j], LK[j]], dtype=float)
        best = v if best is None else min(best, v)
        span = (c - a) / 4
    return best, {'mu': math.exp(cen[0]), 'log_b_over_sigma2': cen[1] + 2 * ell, 'log_r': cen[2], 'kappa': math.exp(cen[3])}


if __name__ == '__main__':
    res = {}
    ells = (55.26204223, 110.52408446, 221.04816893, 455.26, 690.77552790, 960.8)
    Cs = (1.0, 1.41421356, 2.0, 4.0, 10.0, 20.0, 50.0, float('inf'))
    for ell in ells:
        T10 = 1 + 12.1 / ell + 2.5e-4; row = {'T10': T10}
        for Cj in Cs:
            v, at = env_kappa(ell, Cj); row['C=%s' % ('inf' if not np.isfinite(Cj) else '%g' % Cj)] = v; row['argmin_C=%s' % ('inf' if not np.isfinite(Cj) else '%g' % Cj)] = at
        lo, hi = 0.3, 30.0                      # the largest C of g = min(1, C sqrt kappa) with envelope >= T10
        for _ in range(10):
            m = math.sqrt(lo * hi)
            if env_kappa(ell, m, n=16, rounds=6)[0] >= T10: lo = m
            else: hi = m
        row['C_max'] = lo
        for Dd in (0.0, 2.0, 5.0): row['T_D=%g' % Dd] = 1 + (2.1 + Dd) / ell + 2.5e-4
        res[str(round(ell, 2))] = row
        print(round(ell, 2), 'T10 %.4f' % T10, ' '.join('C=%s:%.4f' % (k[2:], v) for k, v in row.items() if k.startswith('C=')), 'C_max %.2f' % lo, flush=True)
    json.dump(res, open('results/cn_envloc.json', 'w'), indent=1)
