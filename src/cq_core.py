"""Stage CQ: the covering holes of paper H v1.0 closed.  Python 3.11.
  (1) mu tails of certificates (i) and (ii): boxes mu in (0, e^-30] and [mu_hi, inf)  (K1, K2), and for (i) the strip b in (0, e^{-2100-3 ell_b}] (K3);
  (2) the wedge of certificate (ii): boxes in (mu, r, beta, ell), beta = b/sigma^2 in [1e4, 1e4 e^{2(ell_b - ell_a)}] (K4, K5).
The term tables are those of v1.0 (cl_link.phi_cells, cl_cells.cell_ok) and stage CP (cp_link.phi_cells_beta), used through cq_phi.py (generated
by cq_make.py with asserted textual replacements).  The rounding rules are those of CL (cl_link for (ii), cl_cells for (i)).

G TERMS (g_terms).  Every term containing G uses either G alone or G/mu.  With G evaluated at G_cell (concavity in G, as in v1.0):
  body  (finite mu box): G = Pi*(mu_a) g_s(r_b) (mu_a >= 2) or mu_b Qb^2 g_s(r_b) (mu < 2), times e^{1e-6}; G/mu with 1/mu_b (positive terms) and 1/mu_a
        (negative terms) -- exactly v1.0;
  upper ([mu_a, inf), mu_a >= ln 2... in fact mu_a = mu_hi >= 50): G = Pi*(mu_a) g_s(r_b) e^{1e-6} (Pi* decreasing for mu >= 2), the positive G/mu
        term is bounded below by 0 (dropped), the negative G/mu term by Pi*(mu_a) g_s/mu_a  (K1);
  lower ((0, e^-30]): G = mu rho with rho in [0, rho_cell], rho_cell = Qb^2 g_s(r_b) e^{1e-6} (P <= mu Qb^2); for fixed A the bound is affine in
        (G, G/mu), hence affine in mu at fixed rho and in rho at fixed mu: its minimum over the box is at a corner (mu, rho) in {0, e^-30} x {0, rho_cell}.
        Three evaluations: G = 0 (rho = 0), (G = 0, G/mu = rho_cell) (mu -> 0) and (G = e^-30 rho_cell, G/mu = rho_cell)  (K2).
ROUNDING.  The rules of CL apply unchanged: an infinite mu or b end never enters a finite operand (it only produces exact -inf terms, which the
log-sum-exp and the netting drop exactly), so the bound B uses the finite ends only; in the ell-variable boxes every quantity that depends on ell
(gamma, U1, U2, gamma_P, c_chi, sigma ends) is computed by the same formulas from the box's own (exact float) ell ends, with at most 6 further
operations, whose error (< 20u relative) is inside the 1000u allowance of every term.
-> results/cq_tails_ii.jsonl, results/cq_tails_i.jsonl, results/cq_wedge.jsonl (one line per column, resumable).
CLI: python cq_core.py tails_ii | tails_i [WORKER NWORKERS] | wedge [ELL_FROM ELL_TO]"""
import json, math, os, sys, time
import numpy as np
import cl_link as K, cl_cells as LC, cc_cells as CC, cc_core as C
import cq_phi as Q

QB = K.QB; LN2 = math.log(2.0); LQB2 = 2 * math.log(QB)
WRONG = {'W1': False, 'W2': False, 'W5': False}     # deliberately wrong versions (search plan sec. 3); all False in the certificate
BODY, LOWER, UPPER = 0, 1, 2
LBMAX = math.log(0.003 / 20)


def lgs(lr_b):
    return 0.5 * (lr_b - np.log(2 + np.exp(lr_b)))


def g_terms(mode, lmu_a, lmu_b, lr_a, lr_b, corner=1):
    """(lGpos, lGneg, lg0_lo, lg0_hi) for G = G_cell, per row; mode array (BODY / LOWER / UPPER).
    LOWER rows are evaluated at the two corners mu = mu_a (corner 0; mu_a = 0 gives G = 0) and mu = mu_b (corner 1) with G = mu rho_cell and
    G/mu = rho_cell (K2: for fixed A the bound is affine in G and in G/mu, so with G = mu rho it is affine in mu for fixed rho and affine in rho for
    fixed mu; the minimum over A is concave along both axes and its minimum over the rectangle (mu, rho) is at a corner; rho = 0 is the G = 0 case)."""
    r_lo, r_hi = np.exp(lr_a), np.exp(lr_b); n = len(lmu_a)
    lGpos = np.full(n, -np.inf); lGneg = np.full(n, -np.inf); lg0_lo = np.full(n, -np.inf); lg0_hi = np.full(n, -np.inf)
    b = mode == BODY
    if b.any():
        la, lb = lmu_a[b], lmu_b[b]
        lPi = np.where(la >= math.log(2.0), CC.lPistar_vec(np.maximum(la, math.log(2.0))), lb + LQB2) + 1e-6
        lG = lPi + lgs(lr_b[b])
        lGpos[b] = lG; lGneg[b] = lG; lg0_lo[b] = lG - LN2 - lb - np.log(2 + r_hi[b]); lg0_hi[b] = lG - LN2 - la - np.log(2 + r_lo[b])
    lo = mode == LOWER
    if lo.any():
        assert np.all(lmu_b[lo] <= math.log(2.0))
        lrho = LQB2 + 1e-6 + lgs(lr_b[lo])
        lmc = lmu_b[lo] if (corner == 1 or WRONG['W1']) else lmu_a[lo]
        lGpos[lo] = lmc + lrho; lGneg[lo] = lmc + lrho
        lg0_lo[lo] = lrho - LN2 - np.log(2 + r_hi[lo]); lg0_hi[lo] = lrho - LN2 - np.log(2 + r_lo[lo])
    up = mode == UPPER
    if up.any():
        la = lmu_a[up]; assert np.all(la >= math.log(2.0))
        lPi = CC.lPistar_vec(la + (60.0 if WRONG['W2'] else 0.0)) + 1e-6
        lG = lPi + lgs(lr_b[up])
        lGpos[up] = lG; lGneg[up] = lG; lg0_lo[up] = -np.inf; lg0_hi[up] = lG - LN2 - la - np.log(2 + r_lo[up])
    return lGpos, lGneg, lg0_lo, lg0_hi


# ------------------------------------------------------------------ certificate (ii), b variable (v1.0 table): mu tails
def ok_ii(X, ell_a, ell_b, Dd=10.0):
    """X rows (lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, mode)."""
    m = X[:, 6].astype(int); args = [X[:, j] for j in range(6)]
    ok = np.ones(len(X), dtype=bool)
    with np.errstate(all='ignore'):
        ok &= Q.phi_b(*args, ell_a, ell_b, Dd, None) > 1e-10
        for cn in ((0, 1) if (m == LOWER).any() else (1,)):
            ok &= Q.phi_b(*args, ell_a, ell_b, Dd, g_terms(m, X[:, 0], X[:, 1], X[:, 2], X[:, 3], cn)) > 1e-10
    return ok


# ------------------------------------------------------------------ certificate (i) (cl_cells table): mu tails and b -> 0
def lae2s(a, b):
    """log(e^a + e^b) that is exact (-inf) when both are -inf."""
    m = np.maximum(a, b); mm = np.where(np.isfinite(m), m, 0.0)
    return np.where(np.isfinite(m), mm + np.log(np.exp(a - mm) + np.exp(b - mm)), m)


def ok_i(X, ell_a, ell_b, delta):
    """cl_cells.cell_ok with the G terms of g_terms; X rows (lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, mode); lb_a may be -inf (b -> 0 strip, K3)."""
    lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b = (X[:, j] for j in range(6)); m = X[:, 6].astype(int)
    gam = (1 + delta) / math.sqrt(ell_a); lgam = math.log(gam)
    lsig_lo = -ell_b; sig_hi = math.exp(-ell_a)
    r_hi = np.exp(lr_b)
    la0 = lr_a + math.log(2 - gam * sig_hi)
    lk0 = lb_a - LN2 * 3 - lr_b
    lk1 = -math.log(16) - np.log(2 + r_hi)
    fin = lambda v: np.where(np.isfinite(v), np.abs(v), 0.0)
    B = 2 * (fin(lr_a) + fin(lr_b) + fin(lb_a) + fin(lb_b) + fin(lmu_a) + fin(lmu_b)) + ell_b + 100
    out = np.ones(len(lmu_a), dtype=bool)
    def accept(lpos, lneg):
        dl = 60 * LC.U * (fin(lpos) + fin(lneg) + B)
        return lpos >= lneg + np.maximum(1e-12, LC.SAFE * dl)
    lT1 = LN2 + 0.5 * (la0 + lk0)                                                       # G = 0
    out &= accept(lae2s(lT1, np.full_like(lT1, lsig_lo + math.log(3 / 8))), lgam + lb_b)
    for cn in ((0, 1) if (m == LOWER).any() else (1,)):
        lGpos, lGneg, lg0_lo, lg0_hi = g_terms(m, lmu_a, lmu_b, lr_a, lr_b, cn)
        # a1 G = (1 - gamma sigma/2) G/(mu (2 + r)); in the body mode exactly v1.0's la1 + lG, in the tails 2 (1 - gamma sigma_hi/2) [G/(2 mu (2+r))]
        la1G = np.where(m == BODY, math.log(1 - gam * sig_hi / 2) - lmu_b - np.log(2 + r_hi) + lGpos, math.log(1 - gam * sig_hi / 2) + LN2 + lg0_lo)
        lT1 = LN2 + 0.5 * (lae2s(la0, la1G) + lae2s(lk0, lk1 + lGpos))
        lT3 = lGpos + lsig_lo - lb_b - math.log(16)
        lpos = lae2s(lae2s(lT1, np.full_like(lT1, lsig_lo + math.log(3 / 8))), lT3)
        lneg = lae2s(lgam + lb_b, lgam + lGneg - LN2)
        out &= accept(lpos, lneg)
    return out


# ------------------------------------------------------------------ the wedge: (mu, r, beta, ell)
def skip_w(X, ell_col_a):
    """(K5): a box whose states all have b >= 1e4 sigma_top^2 (sigma_top = e^{-ell_col_a}) belongs to the v1.0 b boxes (+ the mu tails of (ii))."""
    le = X[:, 6] if WRONG['W5'] else X[:, 7]
    return X[:, 4] >= math.log(1e4) + 2 * (le - ell_col_a) + 1e-9


def ok_w(X, Dd=10.0):
    """X rows (lmu_a, lmu_b, lr_a, lr_b, lB_a, lB_b, ell_a, ell_b, mode)."""
    m = X[:, 8].astype(int); args = [X[:, j] for j in range(6)]
    ok = np.ones(len(X), dtype=bool)
    with np.errstate(all='ignore'):
        ok &= Q.phi_bl(*args, X[:, 6], X[:, 7], Dd, None) > 1e-10
        for cn in ((0, 1) if (m == LOWER).any() else (1,)):
            ok &= Q.phi_bl(*args, X[:, 6], X[:, 7], Dd, g_terms(m, X[:, 0], X[:, 1], X[:, 2], X[:, 3], cn)) > 1e-10
    return ok


# ------------------------------------------------------------------ adaptive bisection (shared)
def split(bad, cols_w, pseudo=0.0):
    """bisect each bad box along its widest weighted side; a side with an infinite end (mu tail, b -> 0 strip) has width 0: it is never split (the box is refined in the other variables).
    cols_w: list of (index of lower end, weight).  Returns new boxes, smallest widest-side width."""
    W = []
    for j, wt in cols_w:
        a, b = bad[:, j], bad[:, j + 1]
        wdt = np.where(np.isfinite(a) & np.isfinite(b), b - a, pseudo)
        W.append(wdt * wt)
    W = np.stack(W, axis=1); k = np.argmax(W, axis=1); ms = float(np.min(np.max(W, axis=1)))
    new = []
    for ii, (j, wt) in enumerate(cols_w):
        sel = bad[k == ii]
        if not len(sel): continue
        a, b = sel[:, j], sel[:, j + 1]
        mid = np.where(np.isfinite(a) & np.isfinite(b), 0.5 * (a + b), np.where(np.isfinite(a), a + pseudo, b - pseudo))
        L = sel.copy(); R = sel.copy(); L[:, j + 1] = mid; R[:, j] = mid; new += [L, R]
    return np.vstack(new), ms


def run_boxes(X0, okf, cols_w, max_cells=300_000_000, chunk=200_000, minw=1e-6, skipf=None):
    stack = [X0]; total = 0; skipped = 0; t0 = time.time()
    while stack:
        X = stack.pop()
        if len(X) > chunk: stack.extend(np.array_split(X, int(math.ceil(len(X) / chunk)))); continue
        if skipf is not None:
            s = skipf(X); skipped += int(s.sum()); X = X[~s]
            if not len(X): continue
        total += len(X)
        if total > max_cells: return {'closed': False, 'cells': total, 'reason': 'cell budget', 'seconds': round(time.time() - t0, 1)}
        ok = okf(X); bad = X[~ok]
        if not len(bad): continue
        new, ms = split(bad, cols_w)
        if ms < minw:
            W = [np.where(np.isfinite(bad[:, j]) & np.isfinite(bad[:, j + 1]), bad[:, j + 1] - bad[:, j], 0.0) * wt for j, wt in cols_w]
            jj = int(np.argmin(np.max(np.stack(W, axis=1), axis=1)))
            return {'closed': False, 'cells': total, 'reason': 'cells too small', 'stuck_cell': [float(v) for v in bad[jj]], 'seconds': round(time.time() - t0, 1)}
        stack.append(new)
    return {'closed': True, 'cells_evaluated': total, 'skipped_boxes': skipped, 'seconds': round(time.time() - t0, 1)}


def mu_hi_of(ell_b): return max(4.0 * ell_b, 50.0)


def tails_ii_boxes(ell_a, ell_b):
    lsig3 = -3 * ell_b; lr = np.linspace(lsig3, math.log(0.25), 25); lb = np.linspace(math.log(1e4) - 2 * ell_a, LBMAX, 25)
    out = []
    for (m0, m1, md) in ((-np.inf, -30.0, LOWER), (math.log(mu_hi_of(ell_b)), np.inf, UPPER)):
        for c, d in zip(lr[:-1], lr[1:]):
            for e, f in zip(lb[:-1], lb[1:]): out.append((m0, m1, c, d, e, f, md))
    return np.array(out)


def tails_i_boxes(ell_a, ell_b):
    lsig3 = -3 * ell_b; lr = np.linspace(lsig3, math.log(0.25), 25); lb = list(np.linspace(-3 * ell_b, 0.0, 25))
    bstrips = [(e, f) for e, f in zip(lb[:-1], lb[1:])] + [(-700.0 * 3 - 3 * ell_b, -3 * ell_b), (-np.inf, -700.0 * 3 - 3 * ell_b)]
    out = []
    for (m0, m1, md) in ((-np.inf, -30.0, LOWER), (math.log(mu_hi_of(ell_b)), np.inf, UPPER)):
        for c, d in zip(lr[:-1], lr[1:]):
            for e, f in bstrips: out.append((m0, m1, c, d, e, f, md))
    # (K3) the b -> 0 strip for the v1.0 mu boxes: below e^{-2100-3 ell_b} for the 16 body mu boxes (whose v1.0 b tail reaches e^{-2100-3 ell_b}),
    # and below sigma_bot^3 for the two finite v1.0 mu tails [e^-30, 1e-4] and [mu_hi, mu_hi e^60] (v1.0 has no b tail for them: a further gap)
    lm = np.linspace(math.log(1e-4), math.log(mu_hi_of(ell_b)), 17)
    for m0, m1 in zip(lm[:-1], lm[1:]):
        for c, d in zip(lr[:-1], lr[1:]): out.append((m0, m1, c, d, -np.inf, -700.0 * 3 - 3 * ell_b, BODY))
    for m0, m1 in ((-30.0, math.log(1e-4)), (math.log(mu_hi_of(ell_b)), math.log(mu_hi_of(ell_b)) + 60)):
        for c, d in zip(lr[:-1], lr[1:]):
            for e, f in ((-np.inf, -700.0 * 3 - 3 * ell_b), (-700.0 * 3 - 3 * ell_b, -3 * ell_b)): out.append((m0, m1, c, d, e, f, BODY))
    return np.array(out)


def wedge_boxes(ell_a, ell_b, nl=10, nB=10):
    lsig3 = -3 * ell_b; lr = np.linspace(lsig3, math.log(0.25), 25)
    lBtop = math.log(1e4) + 2 * (ell_b - ell_a); lB = np.linspace(math.log(1e4), lBtop, nB + 1); le = np.linspace(ell_a, ell_b, nl + 1)
    lm = np.linspace(math.log(1e-4), math.log(mu_hi_of(ell_b)), 17)
    mus = [(-np.inf, -30.0, LOWER), (-30.0, math.log(1e-4), BODY)] + [(a, b, BODY) for a, b in zip(lm[:-1], lm[1:])] + [(math.log(mu_hi_of(ell_b)), np.inf, UPPER)]
    out = []
    for m0, m1, md in mus:
        for c, d in zip(lr[:-1], lr[1:]):
            for e, f in zip(lB[:-1], lB[1:]):
                for g, h in zip(le[:-1], le[1:]): out.append((m0, m1, c, d, e, f, g, h, md))
    return np.array(out)


def certify_tails_ii(ell_a, ell_b, Dd=10.0):
    r = run_boxes(tails_ii_boxes(ell_a, ell_b), lambda X: ok_ii(X, ell_a, ell_b, Dd), [(0, 1.0), (2, 0.7), (4, 0.7)])
    return {'cert': 'ii', 'ell_a': ell_a, 'ell_b': ell_b, 'D': Dd, **r}


def certify_tails_i(ell_a, ell_b, delta, kind):
    r = run_boxes(tails_i_boxes(ell_a, ell_b), lambda X: ok_i(X, ell_a, ell_b, delta), [(0, 1.0), (2, 0.7), (4, 0.7)], minw=1e-7)
    return {'cert': kind, 'ell_a': ell_a, 'ell_b': ell_b, 'delta': delta, **r}


def certify_wedge(ell_a, ell_b, Dd=10.0, wl=2.0, max_cells=300_000_000):
    X0 = wedge_boxes(ell_a, ell_b)
    r = run_boxes(X0, lambda X: ok_w(X, Dd), [(0, 1.0), (2, 0.7), (4, 0.7), (6, wl)], skipf=lambda X: skip_w(X, ell_a), max_cells=max_cells)
    return {'cert': 'ii-wedge', 'ell_a': ell_a, 'ell_b': ell_b, 'D': Dd, 'initial_boxes': len(X0), **r}


def done(fn, key):
    have = set()
    if os.path.exists(fn):
        for l in open(fn):
            if l.strip(): d = json.loads(l); have.add(key(d))
    return have


if __name__ == '__main__':
    K.set_tau_min(20.0); what = sys.argv[1]
    cols = sorted((d['ell_a'], d['ell_b']) for d in (json.loads(l) for l in open('results/cl_link_cells.jsonl') if l.strip()))
    if what == 'tails_ii':
        fn = 'results/cq_tails_ii.jsonl'; have = done(fn, lambda d: round(d['ell_a'], 6))
        for a, b in cols:
            if round(a, 6) in have: continue
            d = certify_tails_ii(a, b); open(fn, 'a').write(json.dumps(d) + '\n'); print(json.dumps(d), flush=True)
            if not d['closed']: print('STOP', flush=True); break
    elif what == 'tails_i':
        w, nw = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (0, 1)
        rows = sorted([json.loads(l) for l in open('results/cl_main_cells.jsonl') if l.strip()], key=lambda r: (r['cert'], r['ell_a']))
        fn = 'results/cq_tails_i_w%d.jsonl' % w; have = done(fn, lambda d: (d['cert'], round(d['ell_a'], 6)))
        for i, r in enumerate(rows):
            if i % nw != w or (r['cert'], round(r['ell_a'], 6)) in have: continue
            d = certify_tails_i(r['ell_a'], r['ell_b'], r['delta'], r['cert']); open(fn, 'a').write(json.dumps(d) + '\n'); print(json.dumps(d), flush=True)
            if not d['closed']: print('STOP', flush=True); break
    elif what == 'wedge':
        lo, hi = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else (0, 1e9)
        fn = os.environ.get('CQ_WEDGE_FN', 'results/cq_wedge.jsonl'); have = done(fn, lambda d: round(d['ell_a'], 6))
        for a, b in cols:
            if not (lo - 1e-9 <= a < hi) or round(a, 6) in have: continue
            d = certify_wedge(a, b); open(fn, 'a').write(json.dumps(d) + '\n'); print(json.dumps({k: v for k, v in d.items() if k != 'stuck_cell'} | ({'stuck': d['stuck_cell']} if 'stuck_cell' in d else {})), flush=True)
            if not d['closed']: print('STOP', flush=True); break
