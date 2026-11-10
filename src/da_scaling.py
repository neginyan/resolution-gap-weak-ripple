"""Stage DA 2.1: THE SCALING STEP OF PAPER H v1.1 CLOSED BY A CERTIFICATE.  Python 3.11 grammar.
The published code is imported UNCHANGED (from DA_H11REPO/src if set, else the bundled copy v11src_ct/src); only the acceptance function of each
certificate is wrapped (as in stage CT), and the published switch FAC (the sandwich factors per channel, used as it is by certificate (iii)) is set.

WHAT IS CERTIFIED.  Certificate (ii) of paper H (the Lemma R* boxes of cl_link.certify, 92 columns; the wedge of cq_core.certify_wedge, 92 columns;
the mu tails of cq_core.certify_tails_ii, 92 columns) is run again with
    FAC = {low_S0: LAM, low_S2: LAM, up_S2: 1},   LAM = 2/3,   target  gamma' = T_{D'}(ell)/sqrt(ell)  with  D' = (LAM T10(ell_a) - 1 - 2.5e-4) ell_a - 2.1,
ell_a the start of the column.  T_{D'}(ell) >= LAM T10(ell) for every ell >= ell_a (report DA, sec. 1, step 4), so a closed column proves
    T(LAM, 1) := P_bg + LAM P_s - N - LAM gamma sigma (d_bg + d_S2) >= 0      on every state of the column,
in the notation of the report (P_bg the gains without a sandwich factor, P_s those with one, N the remainders, d_bg + d_S2 the den bound).
Together with the published certificate (LAM = 1, which is T(1, 1) >= 0) and the affinity of T in LAM, this gives T(lambda, Upsilon) >= 0 for every
admissible pair (report DA, sec. 1).
Each accepted cell must pass at G = 0, at G = G_cell and at the limit G -> 0+ (the G > 0-only term switched on; stage CT).  On a random 10 % of the
accepted cells a middle G is also evaluated (recorded only): G = G_cell e^-t, t ~ U(0.1, 700) and U(0.1, 40); in the lower mu tail of cq_core also
G further by e^-s, s ~ U(0, 700), at fixed G/mu (as CT).  The smallest normalised value over the accepted cells is recorded (the margin).
Modes (one line per column, resumable, worker files in DA_OUT or results/):
  cl W NW      Lemma R* boxes (cl_link.certify)            -> da_cl_w{W}.jsonl
  wedge W NW   the wedge (cq_core.certify_wedge)           -> da_wedge_w{W}.jsonl
  tails W NW   the mu tails (cq_core.certify_tails_ii)     -> da_tails_w{W}.jsonl
  equiv        the wrapped acceptance with LAM = 1 and without the G -> 0+ part is the published acceptance (bit for bit, on sampled cells)
  power        a target slightly above what the factors allow (LAM' = LAM, target LAM_hi T10 with LAM_hi = 1, i.e. no scaling) must FAIL on some
               column (the run is not vacuous), and the wrapper must be able to reject cells
Environment: DA_H11REPO (published repository folder), DA_OUT (output folder), DA_LAM (default 2/3), DA_MAXCOL (at most so many columns per call),
DA_MAXSEC (no new column is started after so many seconds; a column already started runs to its end).
"""
import json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ['DA_H11REPO'] if os.environ.get('DA_H11REPO') else os.path.dirname(HERE)   # v1.2: the repository itself (src/ is the v1.1 code, unchanged)
sys.path.insert(0, os.path.join(ROOT, 'src'))
RES = os.path.join(ROOT, 'results')
OUTD = os.environ.get('DA_OUT', os.path.join(ROOT, 'results'))
import cl_link as KL, cq_core as CQ, cq_phi as Q

LAM = float(os.environ.get('DA_LAM', repr(2.0 / 3.0)))
ACC = 1e-10
MID_RATE = 0.1


def T10(ell): return 1 + 12.1 / ell + 2.5e-4


def Dprime(ell_a, lam=None):
    lam = LAM if lam is None else lam
    return (lam * T10(ell_a) - 1 - 2.5e-4) * ell_a - 2.1


def set_fac(lam):
    KL.FAC.update(low_S0=lam, low_S2=lam, up_S2=1.0)


class Stat:
    def __init__(self, seed):
        self.rng = np.random.default_rng(seed); self.min_margin = math.inf; self.min_cell = None
        self.g0p_rows = 0; self.g0p_differs = 0; self.g0p_rejected = 0
        self.mid_tested = 0; self.mid_negative = 0; self.mid_min = math.inf; self.mid_worst = None

    def acc(self, vals, ok, X):
        if ok.any():
            m = np.min(np.stack(vals), axis=0); m = np.where(ok, m, np.inf); j = int(np.argmin(m))
            if m[j] < self.min_margin: self.min_margin = float(m[j]); self.min_cell = [float(v) for v in X[j]]

    def mid(self, v, rows):
        if not len(v): return
        self.mid_tested += int(len(v)); self.mid_negative += int(np.sum(~(v > ACC)))
        j = int(np.argmin(v))
        if v[j] < self.mid_min: self.mid_min = float(v[j]); self.mid_worst = [float(x) for x in rows[j]]

    def rec(self):
        return {'min_margin_on_accepted': self.min_margin, 'min_margin_cell': self.min_cell, 'g0p_rows': self.g0p_rows,
                'g0p_rows_differing_from_G0': self.g0p_differs, 'cells_rejected_only_by_G0plus': self.g0p_rejected,
                'mid_tested': self.mid_tested, 'mid_negative': self.mid_negative, 'mid_min': self.mid_min, 'mid_worst_row': self.mid_worst}


# ------------------------------------------------------------------ the Lemma R* boxes (cl_link.certify)
def make_cell_ok(stat, g0plus=True, mid=True):
    """cl_link.cell_ok written out (the same lines) with the G -> 0+ end and the records."""
    def ok_f(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, Dd):
        lPi = np.where(lmu_a >= math.log(2.0), KL.CC.lPistar_vec(np.maximum(lmu_a, math.log(2.0))), lmu_b + 2 * math.log(KL.QB)) + 1e-6
        lGc = lPi + 0.5 * (lr_b - np.log(2 + np.exp(lr_b)))
        v0 = KL.phi_cells(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, Dd, -np.inf)
        v1 = KL.phi_cells(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, Dd, lGc)
        ok = (v0 > ACC) & (v1 > ACC); vals = [v0, v1]
        if g0plus:
            v0p = KL.phi_cells(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, Dd, lGc - math.inf)
            stat.g0p_rows += len(v0p); stat.g0p_differs += int(np.sum(v0p != v0)); stat.g0p_rejected += int(np.sum(ok & ~(v0p > ACC)))
            ok &= v0p > ACC; vals.append(v0p)
        X = np.stack([lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b], axis=1); stat.acc(vals, ok, X)
        if mid:
            rows = np.flatnonzero(ok & (stat.rng.random(len(ok)) < MID_RATE))
            if len(rows):
                A = [x[rows] for x in (lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b)]
                for lo, hi in ((0.1, 700.0), (0.1, 40.0)):
                    t = stat.rng.uniform(lo, hi, len(rows))
                    stat.mid(KL.phi_cells(*A, ell_a, ell_b, Dd, lGc[rows] - t), X[rows])
        return ok
    return ok_f


# ------------------------------------------------------------------ the wedge and the mu tails (cq_core)
def shift_gt(gt, t, s=None, lower=None):
    a, b, c, d = (np.asarray(x, dtype=float) - t for x in gt)
    if s is not None:
        a = np.where(lower, a - s, a); b = np.where(lower, b - s, b)
    return a, b, c, d


def make_cq_ok(kind, stat, g0plus=True, mid=True):
    """cq_core.ok_w / ok_ii written out (the same lines) with the G -> 0+ end (per corner) and the records."""
    def ok_f(X, *a, **k):
        if kind == 'wedge':
            Dd = a[0] if a else k.get('Dd', 10.0); m = X[:, 8].astype(int); ext = (X[:, 6], X[:, 7], Dd); phi = Q.phi_bl
        else:
            ell_a, ell_b = a[0], a[1]; Dd = a[2] if len(a) > 2 else k.get('Dd', 10.0); m = X[:, 6].astype(int); ext = (ell_a, ell_b, Dd); phi = Q.phi_b
        args = [X[:, j] for j in range(6)]
        corners = (0, 1) if (m == CQ.LOWER).any() else (1,)
        with np.errstate(all='ignore'):
            v0 = phi(*args, *ext, None); ok = v0 > ACC; vals = [v0]; okp = np.ones(len(X), dtype=bool)
            for cn in corners:
                gt = CQ.g_terms(m, X[:, 0], X[:, 1], X[:, 2], X[:, 3], cn)
                v = phi(*args, *ext, gt); ok &= v > ACC; vals.append(v)
                if g0plus:
                    vp = phi(*args, *ext, shift_gt(gt, math.inf)); okp &= vp > ACC; vals.append(vp)
                    stat.g0p_rows += len(vp); stat.g0p_differs += int(np.sum(vp != v0))
            if g0plus: stat.g0p_rejected += int(np.sum(ok & ~okp)); ok &= okp
            stat.acc(vals, ok, X)
            if mid:
                rows = np.flatnonzero(ok & (stat.rng.random(len(ok)) < MID_RATE))
                if len(rows):
                    Xs = X[rows]; ms = m[rows]; argss = [Xs[:, j] for j in range(6)]
                    exts = (Xs[:, 6], Xs[:, 7], Dd) if kind == 'wedge' else ext
                    for lo, hi in ((0.1, 700.0), (0.1, 40.0)):
                        t = stat.rng.uniform(lo, hi, len(rows)); s = stat.rng.uniform(0.0, 700.0, len(rows))
                        for cn in corners:
                            gt = CQ.g_terms(ms, Xs[:, 0], Xs[:, 1], Xs[:, 2], Xs[:, 3], cn)
                            stat.mid(phi(*argss, *exts, shift_gt(gt, t, s, ms == CQ.LOWER)), Xs)
        return ok
    return ok_f


# ------------------------------------------------------------------ drivers
PUB = {'cl': 'cl_link_cells.jsonl', 'wedge': 'cq_wedge.jsonl', 'tails': 'cq_tails_ii.jsonl'}


def run_column(kind, p, lam, g0plus=True, mid=True, seed=0):
    KL.set_tau_min(20.0); set_fac(lam); Dd = Dprime(p['ell_a'], lam); st = Stat(seed); t0 = time.time()
    try:
        if kind == 'cl':
            orig = KL.cell_ok; KL.cell_ok = make_cell_ok(st, g0plus, mid)
            try: d = KL.certify(p['ell_a'], p['ell_b'], Dd, p['lbmin_rel'])
            finally: KL.cell_ok = orig
        elif kind == 'wedge':
            orig = CQ.ok_w; CQ.ok_w = make_cq_ok('wedge', st, g0plus, mid)
            try: d = CQ.certify_wedge(p['ell_a'], p['ell_b'], Dd)
            finally: CQ.ok_w = orig
        else:
            orig = CQ.ok_ii; CQ.ok_ii = make_cq_ok('ii', st, g0plus, mid)
            try: d = CQ.certify_tails_ii(p['ell_a'], p['ell_b'], Dd)
            finally: CQ.ok_ii = orig
    finally:
        set_fac(1.0)
    cells = d.get('cells_evaluated', d.get('cells'))
    out = {'kind': kind, 'ell_a': p['ell_a'], 'ell_b': p['ell_b'], 'lam': lam, 'up_S2': 1.0, 'D_prime': Dd, 'target_T_at_ell_a': lam * T10(p['ell_a']),
           'closed': bool(d['closed']), 'cells_evaluated': cells, 'published_cells': p['cells_evaluated'], 'cells_minus_published': cells - p['cells_evaluated'],
           'skipped_boxes': d.get('skipped_boxes'), 'published_skipped': p.get('skipped_boxes'), 'reason': d.get('reason'), 'stuck_cell': d.get('stuck_cell'),
           'seconds': round(time.time() - t0, 1), **st.rec()}
    return out


def run(kind, w, nw):
    os.makedirs(OUTD, exist_ok=True)
    pub = [json.loads(l) for l in open(os.path.join(RES, PUB[kind])) if l.strip()]
    fn = os.path.join(OUTD, 'da_%s_w%d.jsonl' % (kind, w))
    have = {round(json.loads(l)['ell_a'], 6) for l in open(fn) if l.strip()} if os.path.exists(fn) else set()
    maxc = int(os.environ.get('DA_MAXCOL', '100000')); nd = 0
    maxs = float(os.environ.get('DA_MAXSEC', 'inf')); t_start = time.time()      # no new column is started after DA_MAXSEC seconds (a run stays within a time slot)
    for i, p in enumerate(pub):
        if i % nw != w or round(p['ell_a'], 6) in have: continue
        if nd >= maxc or time.time() - t_start > maxs: break
        nd += 1
        d = run_column(kind, p, LAM, seed=int(round(p['ell_a'] * 1000)) + {'cl': 1, 'wedge': 7, 'tails': 3}[kind])
        d['index'] = i
        open(fn, 'a').write(json.dumps(d) + '\n'); print(json.dumps({k: d[k] for k in ('kind', 'index', 'ell_a', 'closed', 'cells_evaluated', 'published_cells', 'min_margin_on_accepted', 'mid_negative', 'seconds')}), flush=True)
        if not d['closed']:
            print('STOP: a column did not close (lam = %g)' % LAM, flush=True); break


def equiv(n=200000, seed=1, save=True):
    """the written-out acceptance (without G -> 0+, LAM = 1, published D) equals the published one on random cells of three columns of each kind."""
    rng = np.random.default_rng(seed); out = {}
    KL.set_tau_min(20.0); set_fac(1.0)
    pub = {k: [json.loads(l) for l in open(os.path.join(RES, f)) if l.strip()] for k, f in PUB.items()}
    for i in (0, 45, 91):
        p = pub['cl'][i]; ea, eb = p['ell_a'], p['ell_b']
        lsig3 = -3 * eb; lbmin = -2 * ea + p['lbmin_rel']; lbmax = math.log(0.003 / 20)
        a = rng.uniform(math.log(1e-4), math.log(max(4 * eb, 50)), n); wm = rng.exponential(0.3, n)
        c = rng.uniform(lsig3, math.log(0.25), n); wr = rng.exponential(3, n); e = rng.uniform(lbmin, lbmax, n); wb = rng.exponential(3, n)
        X = [a, np.minimum(a + wm, math.log(max(4 * eb, 50))), c, np.minimum(c + wr, math.log(0.25)), e, np.minimum(e + wb, lbmax)]
        with np.errstate(all='ignore'):
            r1 = KL.cell_ok(*X, ea, eb, p['D']); r2 = make_cell_ok(Stat(0), g0plus=False, mid=False)(*X, ea, eb, p['D'])
        out['cl column %d' % i] = {'cells': n, 'accepted': int(r1.sum()), 'identical': bool(np.array_equal(r1, r2))}
    for kind in ('wedge', 'tails'):
        for i in (0, 91):
            p = pub[kind][i]
            X0 = CQ.wedge_boxes(p['ell_a'], p['ell_b']) if kind == 'wedge' else CQ.tails_ii_boxes(p['ell_a'], p['ell_b'])
            X = X0[rng.choice(len(X0), min(len(X0), 20000), replace=False)]
            if kind == 'wedge':
                X = X[~CQ.skip_w(X, p['ell_a'])]; r1 = CQ.ok_w(X, p['D']); r2 = make_cq_ok('wedge', Stat(0), False, False)(X, p['D'])
            else:
                r1 = CQ.ok_ii(X, p['ell_a'], p['ell_b'], p['D']); r2 = make_cq_ok('ii', Stat(0), False, False)(X, p['ell_a'], p['ell_b'], p['D'])
            out['%s column %d' % (kind, i)] = {'cells': int(len(X)), 'accepted': int(r1.sum()), 'identical': bool(np.array_equal(r1, r2))}
    out['all_identical'] = all(v['identical'] for v in out.values() if isinstance(v, dict))
    print(json.dumps(out))
    if save: json.dump(out, open(os.path.join(OUTD, 'da_equiv.json'), 'w'), indent=1)
    return out


def power():
    """THE RUN IS NOT VACUOUS: on the first column, with LAM = 2/3 in the factors but the UNSCALED target T10 (what one would need if the factors were
    not restored), the column must not close within a modest budget; and with the factors 2/3 and the target 2/3 T10 it closes."""
    pub = [json.loads(l) for l in open(os.path.join(RES, PUB['cl'])) if l.strip()]; p = pub[0]
    KL.set_tau_min(20.0); set_fac(LAM); st = Stat(0); orig = KL.cell_ok; KL.cell_ok = make_cell_ok(st, True, False)
    try:
        d = KL.certify(p['ell_a'], p['ell_b'], p['D'], p['lbmin_rel'], max_cells=3_000_000)
    finally:
        KL.cell_ok = orig; set_fac(1.0)
    out = {'column': 0, 'factors': LAM, 'target': 'T10 (unscaled)', 'closed': bool(d['closed']), 'reason': d.get('reason'), 'cells': d.get('cells', d.get('cells_evaluated'))}
    out['ok'] = not out['closed']
    print(json.dumps(out)); json.dump(out, open(os.path.join(OUTD, 'da_power.json'), 'w'), indent=1)
    return out


def merge():
    out = {'lam': LAM}
    for kind in PUB:
        R = []
        for f in sorted(os.listdir(OUTD)):
            if f.startswith('da_%s_w' % kind) and f.endswith('.jsonl'): R += [json.loads(l) for l in open(os.path.join(OUTD, f)) if l.strip()]
        R.sort(key=lambda d: d['ell_a'])
        out[kind] = {'columns': len(R), 'closed': sum(d['closed'] for d in R), 'cells': sum(d['cells_evaluated'] for d in R), 'published_cells': sum(d['published_cells'] for d in R),
                     'min_margin': min((d['min_margin_on_accepted'] for d in R), default=None), 'mid_tested': sum(d['mid_tested'] for d in R), 'mid_negative': sum(d['mid_negative'] for d in R),
                     'g0p_rows': sum(d['g0p_rows'] for d in R), 'g0p_rows_differing_from_G0': sum(d['g0p_rows_differing_from_G0'] for d in R),
                     'cells_rejected_only_by_G0plus': sum(d['cells_rejected_only_by_G0plus'] for d in R), 'seconds': sum(d['seconds'] for d in R),
                     'not_closed': [d['ell_a'] for d in R if not d['closed']]}
    out['all_closed'] = all(out[k]['columns'] == 92 and out[k]['closed'] == 92 for k in PUB)
    json.dump(out, open(os.path.join(OUTD, 'da_scaling.json'), 'w'), indent=1); print(json.dumps(out, indent=1))
    return out


if __name__ == '__main__':
    m = sys.argv[1]
    if m in PUB: run(m, int(sys.argv[2]), int(sys.argv[3]))
    elif m == 'equiv': equiv()
    elif m == 'power': power()
    elif m == 'merge': merge()
    else: raise SystemExit(__doc__)
