"""Stage CC section 1.1: COMPUTER-ASSISTED PROOF that the main reduced lower bound is at least gamma' = (1 + delta)/sqrt(ell) for every state,
for ell in [ell_lo, ell_hi] (cells in ell), by adaptive cells in (mu, r, b).

The inequality (cc_core.py): g >= gamma'  <=  phi(mu, r, b, G) >= 0 for all mu > 0, r in [sigma^3, 1/4], b in (0, 1], G in [0, Pi*(mu) g_s(r)], with
    phi = 2 sqrt((alpha0 + alpha1 G)(c0 + c1 G)) + 3 sigma/8 + G sigma/(16 b) - gamma b - gamma G/2,
    alpha0 = r (2 - gamma sigma), alpha1 = (1 - gamma sigma/2)/(mu (2 + r)), c0 = b/(8 r), c1 = 1/(16 (2 + r)).
Cell lower bound (every term monotone in every variable; positive terms at their cell minimum, negative terms at their cell maximum, gamma at the
cell's smallest ell, sigma at the cell's largest ell for positive terms and smallest ell inside alpha0, alpha1):
    phi_cell(G) = 2 sqrt((a0 + a1 G)(k0 + k1 G)) + 3 sigma_lo/8 + G sigma_lo/(16 b_hi) - gamma b_hi - gamma G/2,
concave in G, so it suffices to check G = 0 and G = G_cell = Pi*_max(mu cell) sqrt(r_hi/(2 + r_hi)).  Pi*_max over a mu-cell: Pi* is evaluated at
both ends and on 64 interior points; the cell is accepted only if the bound uses the largest of these times (1 + 1e-6) and Pi* is checked to be
monotone (decreasing) on the cell (else the cell is split).  All in logarithms (G ~ e^{-mu} underflows).
Tails: b in (0, b_min] is one cell (G sigma/(16 b) >= G sigma/(16 b_min)); mu in [mu_max, inf): alpha1 >= 0 dropped, G <= G(mu_max); r >= sigma^3.
In-cell audit: random interior points of every accepted cell are evaluated with the exact phi (all variables, G at its exact maximum and at 0).
-> results/cc_cells.jsonl (one line per ell cell, resumable).  Python 3.11."""
import json, math, os, sys
import numpy as np
import cc_core as C

FN = 'results/cc_cells.jsonl'
LN2 = math.log(2.0)


def lae2(a, b):
    m = np.maximum(a, b); return m + np.log(np.exp(a - m) + np.exp(b - m))


def lPistar_vec(lmu):
    """vectorised cc_core.log_Pistar for mu >= 2 (the only use)."""
    mu = np.exp(np.atleast_1d(lmu)); lt0 = math.log(math.e / 2) + np.log(mu - 0.003) + (mu - 0.003) / 2
    la = math.log(C.CS) - mu / 2
    val = 2 * math.log(C.QB) + np.log(mu) + lae2(np.zeros_like(mu), la + lt0) - 2 * lae2(np.zeros_like(mu), lt0)
    return np.minimum(val, np.log(mu) + 2 * math.log(C.QB))


def cell_ok(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, delta):
    """vectorised: arrays of cell bounds (logs).  Returns boolean array (True = phi_cell >= 0 at both G end points)."""
    gam = (1 + delta) / math.sqrt(ell_a); lgam = math.log(gam)
    lsig_lo = -ell_b; sig_hi = math.exp(-ell_a)
    r_hi = np.exp(lr_b)
    la0 = lr_a + math.log(2 - gam * sig_hi)
    la1 = math.log(1 - gam * sig_hi / 2) - lmu_b - np.log(2 + r_hi)
    lk0 = lb_a - LN2 * 3 - lr_b
    lk1 = -math.log(16) - np.log(2 + r_hi)
    # G_cell: max of Pi* over the mu cell (Pi* decreasing for mu >= 2: use the left end; for mu < 2 use mu_b Qbar^2)
    lPi = np.where(lmu_a >= math.log(2.0), lPistar_vec(np.maximum(lmu_a, math.log(2.0))), lmu_b + 2 * math.log(C.QB)) + 1e-6
    lG = lPi + 0.5 * (lr_b - np.log(2 + r_hi))
    out = np.ones(len(lmu_a), dtype=bool)
    for Gz in (False, True):
        if Gz:
            lT1 = LN2 + 0.5 * (lae2(la0, la1 + lG) + lae2(lk0, lk1 + lG))
            lT3 = lG + lsig_lo - lb_b - math.log(16)
            lpos = lae2(lae2(lT1, np.full_like(lT1, lsig_lo + math.log(3 / 8))), lT3)
            lneg = lae2(lgam + lb_b, lgam + lG - LN2)
        else:
            lT1 = LN2 + 0.5 * (la0 + lk0)
            lpos = lae2(lT1, np.full_like(lT1, lsig_lo + math.log(3 / 8)))
            lneg = lgam + lb_b
        out &= lpos >= lneg + 1e-12
    return out


def certify(ell_a, ell_b, delta, max_cells=60_000_000):
    """adaptive bisection in (log mu, log r, log b) for the ell cell [ell_a, ell_b]."""
    lsig3 = -3 * ell_b
    mu_lo, mu_hi = 1e-4, max(4.0 * ell_b, 50.0)
    # initial grid
    lm = np.linspace(math.log(mu_lo), math.log(mu_hi), 17); lr = np.linspace(lsig3, math.log(0.25), 25); lb = np.linspace(-3 * ell_b, 0.0, 25)
    cells = [(a, b, c, d, e, f) for a, b in zip(lm[:-1], lm[1:]) for c, d in zip(lr[:-1], lr[1:]) for e, f in zip(lb[:-1], lb[1:])]
    A = np.array(cells)
    # tails: mu < mu_lo (alpha1 large, G <= mu Qbar^2 small), mu > mu_hi, b < b_min, as extra cells
    tails = [(-30.0, math.log(mu_lo), lsig3, math.log(0.25), -3 * ell_b, 0.0), (math.log(mu_hi), math.log(mu_hi) + 60, lsig3, math.log(0.25), -3 * ell_b, 0.0)]
    tails += [(math.log(mu_lo), math.log(mu_hi), lsig3, math.log(0.25), -700.0 * 3 - 3 * ell_b, -3 * ell_b)]
    A = np.vstack([A, np.array(tails)])
    accepted = 0; total = 0; minside = 1e9; audit = []
    while len(A):
        total += len(A)
        if total > max_cells: return {'ell_a': ell_a, 'ell_b': ell_b, 'delta': delta, 'closed': False, 'cells': total, 'reason': 'cell budget'}
        ok = cell_ok(A[:, 0], A[:, 1], A[:, 2], A[:, 3], A[:, 4], A[:, 5], ell_a, ell_b, delta)
        good = A[ok]; accepted += len(good)
        if len(good):
            idx = np.random.default_rng(len(audit)).integers(0, len(good), min(20, len(good)))
            audit.extend(good[idx].tolist())
        bad = A[~ok]
        if not len(bad): break
        w = np.stack([bad[:, 1] - bad[:, 0], (bad[:, 3] - bad[:, 2]) * 0.7, (bad[:, 5] - bad[:, 4]) * 0.7], axis=1)
        k = np.argmax(w, axis=1); minside = min(minside, float(np.min(np.max(w, axis=1))))
        if minside < 1e-7: return {'ell_a': ell_a, 'ell_b': ell_b, 'delta': delta, 'closed': False, 'cells': total, 'reason': 'cells too small'}
        new = []
        for j in range(3):
            sel = bad[k == j]
            if not len(sel): continue
            mid = (sel[:, 2 * j] + sel[:, 2 * j + 1]) / 2
            L = sel.copy(); R = sel.copy(); L[:, 2 * j + 1] = mid; R[:, 2 * j] = mid
            new += [L, R]
        A = np.vstack(new)
    # in-cell audit with the exact phi (cc_core.phi_over_G at G = G_max) at random interior points
    rng = np.random.default_rng(7); viol = 0; npts = 0
    for c in audit[:3000]:
        for _ in range(3):
            m0, m1 = max(c[0], -30), min(c[1], 30); b0, b1 = max(c[4], -3 * ell_b - 50), c[5]
            if m1 <= m0 or b1 <= b0: continue
            lmu = rng.uniform(m0, m1); lr_ = rng.uniform(c[2], c[3]); lb_ = rng.uniform(b0, b1)
            ell = rng.uniform(ell_a, ell_b); gam = (1 + delta) / math.sqrt(ell)
            try:
                v = C.phi_over_G(lmu, lr_, lb_, -ell, gam)
            except OverflowError:
                v = 1.0
            npts += 1; viol += v < -1e-9
    return {'ell_a': ell_a, 'ell_b': ell_b, 'delta': delta, 'closed': True, 'cells_evaluated': total, 'cells_accepted': accepted, 'audit_points': npts, 'audit_violations': viol}


if __name__ == '__main__':
    ell_from = float(sys.argv[1]); ell_to = float(sys.argv[2]); step = float(sys.argv[3])
    dfix = float(sys.argv[4]) if len(sys.argv) > 4 else 2.5e-4               # delta = 2.1/ell + dfix (eta <= 1/ell, tau >= 20, L^2 <= 1e-4, Lemma R remainders)
    have = set()
    if os.path.exists(FN):
        for l in open(FN):
            if l.strip(): d = json.loads(l); have.add((round(d['ell_a'], 6), d.get('dfix')))
    ea = ell_from
    while ea < ell_to - 1e-9:
        eb = min(ell_to, ea + step)
        if (round(ea, 6), dfix) not in have:
            d = certify(ea, eb, 2.1 / ea + dfix)
            d['dfix'] = dfix
            open(FN, 'a').write(json.dumps(d) + '\n')
            print('ell [%.3f, %.3f] delta %.2e: %s  cells %s  audit %s/%s' % (ea, eb, d['delta'], 'CLOSED' if d['closed'] else 'OPEN (' + d.get('reason', '') + ')',
                  d.get('cells_evaluated', d.get('cells')), d.get('audit_violations'), d.get('audit_points')), flush=True)
        ea = eb
