"""Stage CF 2.3: COMPUTER-ASSISTED PROOF for the 17 'thin block' windows that Corollary S did not cover.  Python 3.11.
For each window (exact sums: S0, S2, eta_lo, eta_hi, sup R'_hi, gap or not) the Lemma S / S' factors low_S0, low_S2, up_S2 are put INTO the cells
(ce_link.FAC: low_S0 on the S0 main term of the numerator, low_S2 on the F-main and W terms, up_S2 on the S2 term of den; every remainder keeps its
etabar = 1/2 inflation and the S0 part of den its factor 2), and the claim
    Num_R* - gamma sigma den_R* >= 0  with  gamma = 1/sqrt(ell)        (i.e. g sqrt(ell) >= 1 - epsbar)
is certified for every state (r, b, A, G) of that window: mu and x are the window's, G <= P g_s(r) exactly (P = x S2), b in [1e4 sigma^2, x/20]
(tau >= 20), 2-D cells in (log r, log b) at the window's sigma (zero-width ell cell).  -> results/cf_windows.jsonl (resumable).
"""
import json, math, os, sys, time
import numpy as np
import ce_link as K, ce_split_map as C

FN = 'results/cf_windows.jsonl'


def certify_window(ell, lmu, lP, lx, factors, lbmin_rel=math.log(1e4), max_cells=200_000_000, chunk=500_000):
    K.FAC.update(factors); K.set_tau_min(20.0)
    Dd = -2.1 - 2.5e-4 * ell                                        # T_D = 1 exactly
    lr = np.linspace(-3 * ell, math.log(0.25), 41); lb = np.linspace(-2 * ell + lbmin_rel, lx - math.log(20), 41)
    A = np.array([(c, d, e, f) for c, d in zip(lr[:-1], lr[1:]) for e, f in zip(lb[:-1], lb[1:])])
    stack = [A]; total = 0
    while stack:
        X = stack.pop()
        if len(X) > chunk: stack.extend(np.array_split(X, int(math.ceil(len(X) / chunk)))); continue
        total += len(X)
        if total > max_cells: return {'closed': False, 'cells': total, 'reason': 'budget'}
        n = len(X); mu = np.full(n, lmu)
        with np.errstate(all='ignore'):
            lG = lP + 0.5 * (X[:, 1] - np.log(2 + np.exp(X[:, 1])))
            v0 = K.phi_cells(mu, mu, X[:, 0], X[:, 1], X[:, 2], X[:, 3], ell, ell, Dd, -np.inf)
            v1 = K.phi_cells(mu, mu, X[:, 0], X[:, 1], X[:, 2], X[:, 3], ell, ell, Dd, lG)
        ok = (v0 > 1e-10) & (v1 > 1e-10)
        bad = X[~ok]
        if not len(bad): continue
        w = np.stack([(bad[:, 1] - bad[:, 0]), (bad[:, 3] - bad[:, 2])], axis=1); k = np.argmax(w, axis=1)
        if np.min(np.max(w, axis=1)) < 1e-6:
            j = int(np.argmin(np.max(w, axis=1)))
            return {'closed': False, 'cells': total, 'reason': 'cells too small', 'stuck': bad[j].tolist()}
        new = []
        for jj in range(2):
            sel = bad[k == jj]
            if not len(sel): continue
            mid = (sel[:, 2 * jj] + sel[:, 2 * jj + 1]) / 2
            L = sel.copy(); R = sel.copy(); L[:, 2 * jj + 1] = mid; R[:, 2 * jj] = mid; new += [L, R]
        stack.append(np.vstack(new))
    return {'closed': True, 'cells_evaluated': total}


if __name__ == '__main__':
    have = set()
    if os.path.exists(FN):
        for l in open(FN):
            if l.strip(): d = json.loads(l); have.add((d['sig_exp'], d['eta_target'], d['u_block_over_u_band'], d['block_width_decay_lengths']))
    ref = {48: (-17.49, 3.92), 96: (-40.88, 6.39), 300: (-142.12, 12.46)}
    for r in json.load(open('results/ce_split_uncovered.json')):
        key = (r['sig_exp'], r['eta_target'], r['u_block_over_u_band'], r['block_width_decay_lengths'])
        if key in have: continue
        e = r['sig_exp']; lx = ref[e][0] * math.log(10); x = 10 ** ref[e][0]; u = ref[e][1]; ell = e * math.log(10)
        ws = C.window_general(x, u, r['eta_target'], r['u_block_over_u_band'] * u, r['block_width_decay_lengths'])
        lam, (l0, u0, l2, u2, f0, f2) = C.lam_S(ws, 10.0 ** (-e))
        S0 = ws['S0_lo'] + ws['S0_hi']; S2 = ws['S2_lo'] + ws['S2_hi']
        fac = {'low_S0': l0, 'low_S2': l2, 'up_S2': u2}
        t0 = time.time()
        d = certify_window(ell, math.log(x * S2 / S0), math.log(x * S2), lx, fac)
        out = {'sig_exp': e, 'eta_target': r['eta_target'], 'u_block_over_u_band': r['u_block_over_u_band'], 'block_width_decay_lengths': r['block_width_decay_lengths'],
               'gap': ws['gap'], 'eta': ws['eta_lo'] + ws['eta_hi'], 'mu': x * S2 / S0, 'factors': fac, 'target': 'g sqrt(ell) >= 1 (times 1 - epsbar)',
               'gsl_split_numerical': r['gsl_split'], 'gsl_hi': r['gsl_hi'], **d, 'seconds': round(time.time() - t0, 1)}
        open(FN, 'a').write(json.dumps(out) + '\n')
        print('1e-%d eta %.3f u_c %.2f width %g mu %.1f factors %.4f %.4f %.4f: %s cells %s' % (e, out['eta'], key[2], key[3], out['mu'], l0, l2, u2,
              'CLOSED' if d['closed'] else 'OPEN ' + d.get('reason', ''), d.get('cells_evaluated', d.get('cells'))), flush=True)
