"""Stage CI 2.1(b): the LOCAL envelope -- the x-free main function of Lemma R_theta with its worst local corrections.  Python 3.11.
Lemma R_theta (report CI 2.1) has the main terms of Lemma R with local (theta_tau-weighted) quantities and three local corrections:
  (L1) w-information of the background: Z J_ww = 1/L^2 + (focus term) + <theta h_w^2/h>, |focus term| <= 2 Qbar sqrt(r/(1+r))/L^2 =: kappa/L^2
       (int Gamma E = sqrt(s^2/(s^2+L^2)) = sqrt(r/(1+r)), |2/L^2 - w^2/L^4| <= 2/L^2 on the E-tilted Gaussian whose variance is < L^2);
  (L2) curvature term of W: |x int Gamma E sum n^2 q_n e_n cos phi_n| <= (2/e) Qbar b = 0.184 b  (x n^2 e_n = b (n^2 tau) e^{-n^2 tau/2} <= 2b/e);
  (L3) local brightness: P_theta <= C Pi*(mu_theta) (local Z3' up to a factor C; measured C <= 2.13 in ci_local.py).
Worst case in the main function (cc_core): N: 2Ar -> 2Ar(1 - kappa), G sigma^2/(16b) -> G sigma^2/(16b) - 0.023 sigma^2 ((1-e)/2 * 0.184 b with 1 - e ~ sigma^2/(4b));
D: Ar -> Ar(1 + kappa), b -> 1.184 b.  Minimised over A exactly, over mu, b (>= 1e4 sigma^2; Lemma BE below), r on grids (long double); largest gamma with phi >= 0 by bisection.
[derivation + numbers; grid minimum, not a certificate]  -> results/ci_envloc.json"""
import json, math
import numpy as np
import cc_core as C
QB = 0.25

def env(ell, Cf, local=True, n=48, rounds=9):
    """min over (log mu, log b, log r) of the main-function lower bound of g sqrt(ell): grid of n^3 points, then zoom around the best point.
    The gamma sigma terms inside alpha are dropped (relative size <= sigma <= e^-55); the A-minimisation is exact (2 sqrt(alpha beta))."""
    LD = np.longdouble; sig = np.exp(LD(-ell)); sell = np.sqrt(LD(ell))
    lo = np.array([math.log(1e-3), -2 * ell + math.log(1e4), -3 * ell]); hi = np.array([math.log(4 * ell + 40), 0.0, math.log(0.25)])
    best = None; cen = None; span = hi - lo
    for k in range(rounds):
        if cen is None: a, c = lo, hi
        else: a = np.maximum(lo, cen - span / 2); c = np.minimum(hi, cen + span / 2)
        g1 = [np.linspace(a[i], c[i], n) for i in range(3)]
        LM, LB, LR = np.meshgrid(*g1, indexing='ij')
        mu = np.exp(LM.astype(LD)); B = np.exp(LB.astype(LD)); rr = np.exp(LR.astype(LD))
        lPi = np.vectorize(C.log_Pistar)(np.exp(g1[0]))
        P = Cf * np.exp(lPi.astype(LD))[:, None, None]
        G = P * np.sqrt(rr / (2 + rr))
        kap = 2 * QB * np.sqrt(rr / (1 + rr)) if local else 0 * rr
        cb = 1.184 if local else 1.0
        alpha = rr * 2 * (1 - kap) + G / (mu * (2 + rr))
        beta = B / (8 * rr) + G / (16 * (2 + rr))
        num = 2 * np.sqrt(alpha * beta) + 3 * sig / 8 + G * sig / (16 * B) - (0.023 * sig if local else 0)
        den = sig * (cb * B + G / 2) + sig * 0                     # Phi/sigma = num - g (cb b + G/2) - g sigma(...)A-terms (dropped, <= sigma)
        val = num / (cb * B + G / 2) * sell                         # g sqrt(ell) at which Phi = 0 (to relative accuracy sigma)
        j = np.unravel_index(int(np.argmin(val)), val.shape)
        v = float(val[j]); cen = np.array([LM[j], LB[j], LR[j]], dtype=float)
        best = v if best is None else min(best, v)
        span = (c - a) / 4
    return best, {'log_mu': cen[0], 'log_b_over_sigma2': cen[1] + 2 * ell, 'log_r': cen[2]}


if __name__ == '__main__':
    res = {}
    for ell in (55.26204223, 110.52408446, 221.04816893, 455.26, 690.77552790, 870.8, 2302.585093, 3453.877639):
        T10 = 1 + 12.1 / ell + 2.5e-4
        v, at = env(ell, 1.0, local=False); row = {'T10': T10, 'period_C1': v, 'argmin_period': at}
        for c in (1.0, 2.13, 3.0, 10.0):
            row['local_C%g' % c] = env(ell, c, local=True)[0]
        res[str(round(ell, 2))] = row
        print(round(ell, 2), {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    json.dump(res, open('results/ci_envloc.json', 'w'), indent=1)
