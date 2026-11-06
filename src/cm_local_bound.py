"""Stage CM 1: the local Z3' in its PROVED form (positivity -> analytic continuation -> log-Bernstein), its numerical check, and the local
envelope with it.  Python 3.11.

Lemma (local brightness; proof in the report).  For every window of the class (f >= 0 of mean 1, K >= x, c_b <= 1, (P3)) and every focus
(t0, tau = x/b, b <= 1), with E0 = <theta_tau R^2> = S0_theta/2, E1 = <theta_tau R'^2> = S2_theta/2, E0' = <theta_{2 tau} R^2>:
    x E1 <= (1 + b/2) E0 ln(2 B0/E0) + (sqrt2/e) b E0',     B0 = (1 + sqrt(E0) + (E0/2)^{1/(2+b)})^2,   E0' <= Qbar^2 = 1/16.
In (mu_theta, P_theta) = (x E1/E0, 2 x E1):
    mu <= (1 + b/2) ln(4 B0 mu/P) + 2 (sqrt2/e) b mu E0'/P.                                                 (*)
Psi(mu, b) := the largest P with (*) (E0' = 1/16); the brightness bound used in the cells is min(2 mu Qbar^2, Psi(mu, b)).
As b -> 0 and for small E0, Psi ~ 4 B0 mu e^{-mu} with B0 -> 1: the form the coherent-packet counterexample needs, with constant ~4-5.
-> results/cm_local_bound.json
"""
import json, math
import numpy as np
import cc_core as C

QB = 0.25
K2 = math.sqrt(2) / math.e          # 0.5203


def rhs(P, mu, b, E0p=1.0 / 16):
    E0 = P / (2 * mu)
    B0 = (1 + math.sqrt(E0) + (E0 / 2) ** (1.0 / (2 + b))) ** 2
    return (1 + b / 2) * math.log(4 * B0 * mu / P) + 2 * K2 * b * mu * E0p / P


def Psi(mu, b, E0p=1.0 / 16):
    """largest P with mu <= rhs(P); rhs is decreasing in P.  Returns log P."""
    lo, hi = -700.0, math.log(2 * mu * QB * QB)
    if rhs(math.exp(hi), mu, b, E0p) >= mu: return hi
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if rhs(math.exp(m), mu, b, E0p) >= mu: lo = m
        else: hi = m
    return lo


def local_quantities(ns, q, al, x, tau):
    """exact S0, S2 at focus variance tau for R = sum q_n cos(n t - phi_n) seen at t0 (al_n = n t0 - phi_n); returns E0, E1."""
    D = ns[:, None] - ns[None, :]; S = ns[:, None] + ns[None, :]
    cd = np.cos(al[:, None] - al[None, :]); cs = np.cos(al[:, None] + al[None, :])
    qq = q[:, None] * q[None, :]; nm = ns[:, None] * ns[None, :]
    ed = np.exp(-D ** 2 * tau / 2); es = np.exp(-S ** 2 * tau / 2)
    S0 = float(np.sum(qq * (ed * cd + es * cs))); S2 = float(np.sum(nm * qq * (ed * cd - es * cs)))
    return S0 / 2, S2 / 2


def check_window(ns, araw, phi, x, b, t0, K=None):
    """araw: raw |a_n| (f = 1 + sum araw cos(n t - phi)); returns the slack of the lemma: RHS - x E1 (>= 0 if the lemma holds)."""
    K = x if K is None else K
    q = araw * np.exp(-ns ** 2 * K / 2); al = ns * t0 - phi; tau = x / b
    E0, E1 = local_quantities(ns, q, al, x, tau); E0p, _ = local_quantities(ns, q, al, x, 2 * tau)
    if E0 <= 0: return None
    B0 = (1 + math.sqrt(E0) + (E0 / 2) ** (1.0 / (2 + b))) ** 2
    R = (1 + b / 2) * E0 * math.log(2 * B0 / E0) + K2 * b * E0p
    return {'xE1': x * E1, 'rhs': R, 'ratio': x * E1 / R, 'mu': x * E1 / E0, 'P': 2 * x * E1, 'E0': E0, 'E0p': E0p}


def random_checks(rng, n_win=400):
    """windows with f >= 0 (checked on a grid), random coefficients and phases, (P2)+(P3) enforced by scaling; several foci and b."""
    worst = 0.0; rows = []
    T = np.linspace(0, 2 * math.pi, 1 << 14, endpoint=False)
    for _ in range(n_win):
        kind = rng.integers(4)
        if kind == 0:     # packet at a node
            N = int(rng.integers(20, 400)); W = int(rng.integers(1, 12)); ns = np.concatenate(([1], np.arange(N, N + W)))
            t0 = math.pi / 2; phi = np.concatenate(([0.0], ns[1:] * t0)); araw = np.concatenate(([rng.uniform(0.05, 0.3)], np.full(W, 1.0)))
        elif kind == 1:   # random band
            N = int(rng.integers(5, 300)); W = int(rng.integers(2, 30)); ns = np.concatenate(([1], np.arange(N, N + W)))
            phi = rng.uniform(0, 2 * math.pi, len(ns)); araw = np.concatenate(([rng.uniform(0.05, 0.3)], rng.uniform(0, 1, W))); t0 = rng.uniform(0, 2 * math.pi)
        elif kind == 2:   # comb-like (many harmonics, equal amplitudes)
            M = int(rng.integers(5, 200)); ns = np.arange(1, M + 1); phi = ns * 0.0 + rng.uniform(0, 2 * math.pi) * 0; araw = np.full(M, 1.0); t0 = rng.uniform(0, 2 * math.pi)
        else:             # few harmonics, random
            M = int(rng.integers(2, 8)); ns = np.arange(1, M + 1); phi = rng.uniform(0, 2 * math.pi, M); araw = rng.uniform(0, 1, M); t0 = rng.uniform(0, 2 * math.pi)
        araw = np.asarray(araw, float); ns = np.asarray(ns)
        # (P2): fundamental largest
        araw[0] = max(araw[0], araw.max())
        # f >= 0 by scaling the whole ripple (f -> 1 + eps (f - 1))
        f = 1 + sum(a * np.cos(k * T - p) for k, p, a in zip(ns, phi, araw))
        m = f.min()
        if m < 0: araw = araw / (1 - m) * 0.999
        x = float(np.exp(rng.uniform(math.log(1e-6), math.log(0.003))))
        # (P3) by scaling again if needed (keeps f >= 0)
        s = float(np.sum(araw * np.exp(-ns ** 2 * x / 2)))
        if s > QB: araw = araw * QB / s * 0.999
        for b in (1.0, 0.3, 0.03, 1e-3):
            if x / b > 25: continue
            r = check_window(ns, araw, phi, x, b, t0)
            if r is None: continue
            worst = max(worst, r['ratio']); rows.append({'kind': int(kind), 'x': x, 'b': b, **r})
    return worst, rows


def env_local(ell, n=48, rounds=9):
    """min of g sqrt(ell) from the local main function with G_max = min(2 mu Qbar^2, Psi(mu, b)) g_s (depends on b).  long double."""
    LD = np.longdouble; sig = np.exp(LD(-ell)); sell = np.sqrt(LD(ell))
    lo = np.array([math.log(1e-3), -2 * ell + math.log(1e4), -3 * ell]); hi = np.array([math.log(4 * ell + 40), 0.0, math.log(0.25)])
    best = None; cen = None; span = hi - lo
    for k in range(rounds):
        if cen is None: a, c = lo, hi
        else: a = np.maximum(lo, cen - span / 2); c = np.minimum(hi, cen + span / 2)
        g1 = [np.linspace(a[i], c[i], n) for i in range(3)]
        LM, LB, LR = np.meshgrid(*g1, indexing='ij')
        mu = np.exp(LM.astype(LD)); B = np.exp(LB.astype(LD)); rr = np.exp(LR.astype(LD))
        lP = np.array([[min(math.log(2 * m * QB * QB), Psi(m, bb)) for bb in np.exp(g1[1])] for m in np.exp(g1[0])])
        P = np.exp(np.asarray(lP, dtype=LD))[:, :, None]
        G = P * np.sqrt(rr / (2 + rr))
        kap = 2 * QB * np.sqrt(rr / (1 + rr)); cb = 1.184
        alpha = rr * 2 * (1 - kap) + G / (mu * (2 + rr))
        beta = B / (8 * rr) + G / (16 * (2 + rr))
        num = 2 * np.sqrt(alpha * beta) + 3 * sig / 8 + G * sig / (16 * B) - 0.023 * sig
        val = num / (cb * B + G / 2) * sell
        j = np.unravel_index(int(np.argmin(val)), val.shape)
        v = float(val[j]); cen = np.array([LM[j], LB[j], LR[j]], dtype=float)
        best = v if best is None else min(best, v)
        span = (c - a) / 4
    return best, {'mu': math.exp(cen[0]), 'log_b_over_sigma2': cen[1] + 2 * ell, 'log_r': cen[2]}


if __name__ == '__main__':
    res = {}
    # 1. the proved bound as a function of mu (b -> 0 and b = 1), against mu e^-mu and Pi*
    tab = []
    for mu in (2, 3, 4, 5, 5.8, 6, 8, 10, 15, 20, 30):
        row = {'mu': mu}
        for b in (1e-12, 1e-3, 0.1, 1.0):
            lP = min(math.log(2 * mu * QB * QB), Psi(mu, b)); row['c_b%g' % b] = math.exp(lP) / (mu * math.exp(-mu))
        row['Pistar_over_mu_e^-mu'] = math.exp(C.log_Pistar(mu)) / (mu * math.exp(-mu))
        tab.append(row); print({k: round(v, 4) for k, v in row.items()})
    res['table'] = tab
    # 2. numerical check of the lemma
    rng = np.random.default_rng(20261104)
    worst, rows = random_checks(rng)
    res['check'] = {'n': len(rows), 'worst_ratio_xE1_over_rhs': worst, 'max_mu': max(r['mu'] for r in rows)}
    print('lemma check: n', len(rows), 'worst x E1 / RHS =', worst, '(must be <= 1)')
    # 3. the envelope with the proved bound
    env = {}
    for ell in (55.26204223, 110.52408446, 221.04816893, 455.26, 690.77552790, 870.8, 960.8):
        T10 = 1 + 12.1 / ell + 2.5e-4
        v, at = env_local(ell)
        env[str(round(ell, 2))] = {'T10': T10, 'env_proved': v, 'margin': v - T10, 'argmin': at}
        print(round(ell, 2), 'T10 %.4f env %.4f margin %.4f' % (T10, v, v - T10), at, flush=True)
    res['envelope'] = env
    json.dump(res, open('results/cm_local_bound.json', 'w'), indent=1)
