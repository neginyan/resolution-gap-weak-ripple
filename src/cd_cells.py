"""Stage CD 1.3 (and the split class of 1.1): COMPUTER-ASSISTED PROOF of the main reduced inequality at the RAISED target
    g sqrt(ell) >= T_D(ell) = 1 + 2.1/ell + 2.5e-4 + D/ell      (D = 10),
for ell in [55.26, 6907.8] (sigma = 1e-24 .. 1e-3000), by the cells of cc_cells.py (same cell test cc_cells.cell_ok, same tails, same in-cell
audit).  Only the bookkeeping differs: the bad cells are processed in chunks (a stack) so that the memory stays bounded at large ell, and the
cell budget is larger.  The room D/ell is what Lemma S's corollary needs (cd_split.py): g sqrt(ell) >= Lambda_S T_D kappa_tau ..., and the
CC counterexample windows have 1 - Lambda_S <= 0.008 < D/ell.  T_D >= T_0, so these cells also extend the CC certificate (D = 0).
-> results/cd_cells.jsonl (one line per ell cell, resumable by (ell_a, D)).  CLI: python cd_cells.py ELL_FROM ELL_TO STEP [D].  Python 3.11."""
import json, math, os, sys, time
import numpy as np
import cc_cells as CC, cc_core as C

FN = 'results/cd_cells.jsonl'
CHUNK = 2_000_000


def certify(ell_a, ell_b, delta, max_cells=1_500_000_000):
    lsig3 = -3 * ell_b
    mu_lo, mu_hi = 1e-4, max(4.0 * ell_b, 50.0)
    lm = np.linspace(math.log(mu_lo), math.log(mu_hi), 17); lr = np.linspace(lsig3, math.log(0.25), 25); lb = np.linspace(-3 * ell_b, 0.0, 25)
    A = np.array([(a, b, c, d, e, f) for a, b in zip(lm[:-1], lm[1:]) for c, d in zip(lr[:-1], lr[1:]) for e, f in zip(lb[:-1], lb[1:])])
    tails = [(-30.0, math.log(mu_lo), lsig3, math.log(0.25), -3 * ell_b, 0.0), (math.log(mu_hi), math.log(mu_hi) + 60, lsig3, math.log(0.25), -3 * ell_b, 0.0),
             (math.log(mu_lo), math.log(mu_hi), lsig3, math.log(0.25), -700.0 * 3 - 3 * ell_b, -3 * ell_b)]
    stack = [np.vstack([A, np.array(tails)])]
    total = 0; accepted = 0; minside = 1e9; audit = []; rng_a = np.random.default_rng(11)
    while stack:
        X = stack.pop()
        if len(X) > CHUNK:
            stack.extend(np.array_split(X, int(math.ceil(len(X) / CHUNK)))); continue
        total += len(X)
        if total > max_cells: return {'ell_a': ell_a, 'ell_b': ell_b, 'delta': delta, 'closed': False, 'cells': total, 'reason': 'cell budget'}
        ok = CC.cell_ok(X[:, 0], X[:, 1], X[:, 2], X[:, 3], X[:, 4], X[:, 5], ell_a, ell_b, delta)
        good = X[ok]; accepted += len(good)
        if len(good) and len(audit) < 6000:
            audit.extend(good[rng_a.integers(0, len(good), min(10, len(good)))].tolist())
        bad = X[~ok]
        if not len(bad): continue
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
        stack.append(np.vstack(new))
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
    ell_from = float(sys.argv[1]); ell_to = float(sys.argv[2]); step = float(sys.argv[3]); Dd = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0
    have = set()
    if os.path.exists(FN):
        for l in open(FN):
            if l.strip(): d = json.loads(l); have.add((round(d['ell_a'], 6), d['D']))
    ea = ell_from
    while ea < ell_to - 1e-9:
        eb = min(ell_to, ea + step)
        if (round(ea, 6), Dd) not in have:
            t0 = time.time(); d = certify(ea, eb, 2.1 / ea + 2.5e-4 + Dd / ea); d['D'] = Dd; d['dfix'] = 2.5e-4; d['seconds'] = round(time.time() - t0, 1)
            open(FN, 'a').write(json.dumps(d) + '\n')
            print('ell [%.3f, %.3f] delta %.3e: %s cells %s audit %s/%s  %.0fs' % (ea, eb, d['delta'], 'CLOSED' if d['closed'] else 'OPEN (' + d.get('reason', '') + ')',
                  d.get('cells_evaluated', d.get('cells')), d.get('audit_violations'), d.get('audit_points'), d['seconds']), flush=True)
        ea = eb
