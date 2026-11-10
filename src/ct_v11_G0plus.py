"""Stage CT 2.1: the G -> 0+ end checked in EVERY published certificate of paper H v1.1 that uses the G-endpoint rule.  Python 3.11.
The published code is imported UNCHANGED (from CT_H11REPO/src if set, else the bundled copy v11src_ct/src); only the acceptance function of each
certificate is wrapped.  Each accepted cell must ALSO pass at G -> 0+ (the limit: every G-term exactly 0 on the G > 0 path; this switches on the one G-free term that
cl_link / cq_phi add only when G > 0, Cn: tnu2 Qb/(16(1+Z)), power 1).  The cells are the published ones, so the cell counts must be the published ones.
On a random fraction of the accepted cells a MIDDLE G is also evaluated (recorded only, never used for acceptance):
  t1 ~ U(0.1, 700) (log-uniform G in [G_cell e^-700, G_cell e^-0.1], as asked) and t2 ~ U(0.1, 40) (denser where the CS band was), G = G_cell e^-t;
  in the lower mu tail (G = mu rho) also mu inside the tail: G shifted by a further s ~ U(0, 700) with G/mu fixed.
Liveness: every record counts the rows on which the G -> 0+ evaluation ran, and how many of them gave a value different from G = 0
(the extra term is really switched on).
Modes (one line per published column/window, resumable, worker files):
  wedge W NW     cq_core.certify_wedge  (certificate (ii), wedge)      -> results/ct_v11_wedge_G0plus_w{W}.jsonl
  tails_ii       cq_core.certify_tails_ii (certificate (ii), mu tails) -> results/ct_v11_tails_ii_G0plus.jsonl
  windows        cf_windows.certify_window with cl_link (certificate (iii), 17 windows) -> results/ct_v11_windows_G0plus.jsonl
  cp_control     cp_link.certify, the two closed control columns (width 1) of results/cp_link_control_w1.jsonl -> results/ct_v11_cp_control_G0plus.jsonl
  cl W NW        cl_link.certify on all 92 published CL columns at the LIMIT G -> 0+ (stage CS checked G_cell e^-700 only), worker W of NW
                 -> results/ct_v11_cl_G0plus_w{W}.jsonl
  cl_mid I,J,..  the same on chosen columns -> results/ct_v11_cl_mid.jsonl
"""
import json, math, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ['CT_H11REPO'] if os.environ.get('CT_H11REPO') else os.path.dirname(HERE)   # v1.2: the repository itself (src/ is the v1.1 code, unchanged)
sys.path.insert(0, os.path.join(ROOT, 'src'))
RES = os.path.join(ROOT, 'results')
OUTD = os.environ.get('CT_OUT', os.path.join(ROOT, 'results'))
import cl_link as KL, cq_core as CQ, cq_phi as Q, cp_link as CP, cf_windows as CF, ce_split_map as SM

SHIFT0 = math.inf        # G -> 0+ is evaluated as the LIMIT: every G-term exactly -inf (exactly 0) on the G > 0 code path (not G0), so the G-free
                         # G > 0-only term is switched on.  A point G_cell e^-700 is not enough: log G_cell can be ~ -10^4 and e^-700 G terms need not be negligible.
MID_RATE = float(os.environ.get('CT_MID_RATE', '0.1'))
ACC = 1e-10


class Stat:
    def __init__(self, seed):
        self.rng = np.random.default_rng(seed); self.g0p_rows = 0; self.g0p_differs = 0; self.g0p_min = math.inf
        self.mid_tested = 0; self.mid_negative = 0; self.mid_min = math.inf; self.mid_worst = None

    def g0p(self, v0, v0p, okall):
        self.g0p_rows += int(len(v0p)); self.g0p_differs += int(np.sum(v0p != v0))
        if okall.any(): self.g0p_min = min(self.g0p_min, float(np.min(v0p[okall])))

    def pick(self, ok):
        return np.flatnonzero(ok & (self.rng.random(len(ok)) < MID_RATE))

    def mid(self, vals, rows, ts):
        if not len(rows): return
        self.mid_tested += int(len(rows)); bad = vals <= ACC; self.mid_negative += int(bad.sum())
        j = int(np.argmin(vals));
        if vals[j] < self.mid_min: self.mid_min = float(vals[j]); self.mid_worst = [float(x) for x in rows[j]] + [float(ts[j])]

    def rec(self):
        return {'g0p_rows': self.g0p_rows, 'g0p_rows_differing_from_G0': self.g0p_differs, 'g0p_min_on_accepted': self.g0p_min,
                'mid_tested': self.mid_tested, 'mid_negative': self.mid_negative, 'mid_min': self.mid_min, 'mid_worst_row_and_t': self.mid_worst}


def shift_gt(gt, t, s=None, lower=None):
    """G-terms (lGpos, lGneg, lg0_lo, lg0_hi) scaled by e^-t (all four: G and G/mu together, i.e. G -> G e^-t at fixed mu);
    with s (lower tail): G further by e^-s at fixed G/mu (mu inside the tail)."""
    a, b, c, d = (np.asarray(x, dtype=float) - t for x in gt)
    if s is not None:
        a = np.where(lower, a - s, a); b = np.where(lower, b - s, b)
    return a, b, c, d


# ------------------------------------------------------------------ certificate (ii): wedge and mu tails (cq_core)
def make_ok(kind, stat, Dd=10.0):
    orig = CQ.ok_w if kind == 'wedge' else CQ.ok_ii
    def ok(X, *a, **k):
        res = orig(X, *a, **k)
        if kind == 'wedge':
            m = X[:, 8].astype(int); ext = (X[:, 6], X[:, 7], Dd); phi = Q.phi_bl
        else:
            ell_a, ell_b = a[0], a[1]; Dd_ = a[2] if len(a) > 2 else k.get('Dd', 10.0); m = X[:, 6].astype(int); ext = (ell_a, ell_b, Dd_); phi = Q.phi_b
        args = [X[:, j] for j in range(6)]
        corners = (0, 1) if (m == CQ.LOWER).any() else (1,)
        with np.errstate(all='ignore'):
            v0 = phi(*args, *ext, None)
            ok0p = np.ones(len(X), dtype=bool)
            for cn in corners:
                gt = CQ.g_terms(m, X[:, 0], X[:, 1], X[:, 2], X[:, 3], cn)
                v0p = phi(*args, *ext, shift_gt(gt, SHIFT0)); ok0p &= v0p > ACC
                stat.g0p(v0, v0p, res & (v0p > ACC))
            okall = res & ok0p
            rows = stat.pick(okall)
            if len(rows):
                Xs = X[rows]; ms = m[rows]; argss = [Xs[:, j] for j in range(6)]
                exts = (Xs[:, 6], Xs[:, 7], Dd) if kind == 'wedge' else ext
                for lo, hi in ((0.1, 700.0), (0.1, 40.0)):
                    t = stat.rng.uniform(lo, hi, len(rows)); s = stat.rng.uniform(0.0, 700.0, len(rows))
                    for cn in corners:
                        gt = CQ.g_terms(ms, Xs[:, 0], Xs[:, 1], Xs[:, 2], Xs[:, 3], cn)
                        v = phi(*argss, *exts, shift_gt(gt, t, s, ms == CQ.LOWER))
                        stat.mid(v, Xs, t)
        return okall
    return ok


def run_cq(kind, w=0, nw=1):
    KL.set_tau_min(20.0)
    src = 'cq_wedge.jsonl' if kind == 'wedge' else 'cq_tails_ii.jsonl'
    pub = [json.loads(l) for l in open(os.path.join(RES, src)) if l.strip()]
    fn = os.path.join(OUTD, 'ct_v11_wedge_G0plus_w%d.jsonl' % w if kind == 'wedge' else 'ct_v11_tails_ii_G0plus.jsonl')
    have = {round(json.loads(l)['ell_a'], 6) for l in open(fn) if l.strip()} if os.path.exists(fn) else set()
    orig_w, orig_ii = CQ.ok_w, CQ.ok_ii; maxc = int(os.environ.get('CT_MAXCOL', '1000000')); nd = 0
    for i, p in enumerate(pub):
        if i % nw != w or round(p['ell_a'], 6) in have: continue
        if nd >= maxc: break
        nd += 1
        st = Stat(int(round(p['ell_a'] * 1000)) + (7 if kind == 'wedge' else 3)); t0 = time.time()
        try:
            if kind == 'wedge':
                CQ.ok_w = make_ok('wedge', st, p['D']); d = CQ.certify_wedge(p['ell_a'], p['ell_b'], p['D'])
            else:
                CQ.ok_ii = make_ok('ii', st, p['D']); d = CQ.certify_tails_ii(p['ell_a'], p['ell_b'], p['D'])
        finally:
            CQ.ok_w, CQ.ok_ii = orig_w, orig_ii
        d.update({'published_cells': p['cells_evaluated'], 'published_skipped': p.get('skipped_boxes'), 'seconds': round(time.time() - t0, 1), **st.rec()})
        d.pop('stuck_cell', None) if d.get('closed') else None
        open(fn, 'a').write(json.dumps(d) + '\n'); print(json.dumps(d), flush=True)
        if not (d['closed'] and d['cells_evaluated'] == p['cells_evaluated']): print('STOP: a column changed', flush=True); break


# ------------------------------------------------------------------ certificate (iii): the 17 windows (cl_link through cf_windows)
class KProxy:
    """cl_link with phi_cells wrapped: the G = G_cell call also requires G -> 0+; middle G recorded on accepted cells."""
    def __init__(self, stat): self.stat = stat; self.last_v0 = None

    def __getattr__(self, a): return getattr(KL, a)

    def phi_cells(self, lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, Dd, lG):
        v = KL.phi_cells(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, Dd, lG)
        if np.isscalar(lG) and not np.isfinite(lG): self.last_v0 = v; return v
        v0p = KL.phi_cells(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, Dd, lG - SHIFT0)
        okall = (self.last_v0 > ACC) & (v > ACC) & (v0p > ACC); self.stat.g0p(self.last_v0, v0p, okall)
        rows = self.stat.pick(okall)
        if len(rows):
            A = [np.asarray(x)[rows] if np.ndim(x) else x for x in (lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b)]
            for lo, hi in ((0.1, 700.0), (0.1, 40.0)):
                t = self.stat.rng.uniform(lo, hi, len(rows))
                vm = KL.phi_cells(*A, ell_a, ell_b, Dd, lG[rows] - t)
                self.stat.mid(vm, np.stack([lr_a[rows], lr_b[rows], lb_a[rows], lb_b[rows]], axis=1), t)
        return np.minimum(v, v0p)


def run_windows():
    fn = os.path.join(OUTD, 'ct_v11_windows_G0plus.jsonl')
    pub = [json.loads(l) for l in open(os.path.join(RES, 'cl_windows.jsonl')) if l.strip()]
    W = [json.loads(l) for l in open(os.path.join(RES, 'cf_windows.jsonl')) if l.strip()]
    ref = {48: (-17.49, 3.92), 96: (-40.88, 6.39), 300: (-142.12, 12.46)}          # as cl_windows_run.py
    have = {json.loads(l)['index'] for l in open(fn) if l.strip()} if os.path.exists(fn) else set()
    saved = dict(KL.FAC); old = CF.K
    for p in pub:
        i = p['index']
        if i in have: continue
        w = W[i]; e = w['sig_exp']; x = 10 ** ref[e][0]; u = ref[e][1]
        ws = SM.window_general(x, u, w['eta_target'], w['u_block_over_u_band'] * u, w['block_width_decay_lengths'])
        S0 = ws['S0_lo'] + ws['S0_hi']; S2 = ws['S2_lo'] + ws['S2_hi']
        st = Stat(1000 + i); prox = KProxy(st); CF.K = prox; t0 = time.time()
        try:
            d = CF.certify_window(e * math.log(10), math.log(x * S2 / S0), math.log(x * S2), math.log(x), w['factors'])
        finally:
            CF.K = old; KL.FAC.clear(); KL.FAC.update(saved)
        rec = {'index': i, 'sig_exp': e, 'closed': d['closed'], 'cells_evaluated': d.get('cells_evaluated', d.get('cells')), 'published_cells': p['cells_evaluated'],
               'seconds': round(time.time() - t0, 1), **st.rec()}
        open(fn, 'a').write(json.dumps(rec) + '\n'); print(json.dumps(rec), flush=True)
        if not (rec['closed'] and rec['cells_evaluated'] == p['cells_evaluated']): print('STOP: a window changed', flush=True); break


# ------------------------------------------------------------------ cp_link (stage CP beta boxes; not used by Theorem H1 in v1.1, see the report)
def run_cp_control():
    fn = os.path.join(OUTD, 'ct_v11_cp_control_G0plus.jsonl')
    pub = [json.loads(l) for l in open(os.path.join(RES, 'cp_link_control_w1.jsonl')) if l.strip()]
    orig = CP.cell_ok
    out = []
    for p in pub:
        KL.set_tau_min(p['tau_min']); st = Stat(2000 + int(p['ell_a']))
        def ok(lmu_a, lmu_b, lr_a, lr_b, lB_a, lB_b, ell_a, ell_b, Dd):
            res = orig(lmu_a, lmu_b, lr_a, lr_b, lB_a, lB_b, ell_a, ell_b, Dd)
            lPi = np.where(lmu_a >= math.log(2.0), CP.CC.lPistar_vec(np.maximum(lmu_a, math.log(2.0))), lmu_b + 2 * math.log(KL.QB)) + 1e-6   # as cp_link.cell_ok
            lGc = lPi + 0.5 * (lr_b - np.log(2 + np.exp(lr_b)))
            v0 = CP.phi_cells_beta(lmu_a, lmu_b, lr_a, lr_b, lB_a, lB_b, ell_a, ell_b, Dd, -np.inf)
            v0p = CP.phi_cells_beta(lmu_a, lmu_b, lr_a, lr_b, lB_a, lB_b, ell_a, ell_b, Dd, lGc - SHIFT0)
            okall = res & (v0p > ACC); st.g0p(v0, v0p, okall)
            rows = st.pick(okall)
            if len(rows):
                A = [x[rows] for x in (lmu_a, lmu_b, lr_a, lr_b, lB_a, lB_b)]
                for lo, hi in ((0.1, 700.0), (0.1, 40.0)):
                    t = st.rng.uniform(lo, hi, len(rows))
                    st.mid(CP.phi_cells_beta(*A, ell_a, ell_b, Dd, lGc[rows] - t), np.stack(A, axis=1), t)
            return okall
        u = CP.certify(p['ell_a'], p['ell_b'], p['D'], p['lBmin'])      # the published code as it is (the stored record predates its last edit)
        CP.cell_ok = ok; t0 = time.time()
        try:
            d = CP.certify(p['ell_a'], p['ell_b'], p['D'], p['lBmin'])
        finally:
            CP.cell_ok = orig
        d.update({'stored_cells': p['cells_evaluated'], 'unwrapped_cells_now': u.get('cells_evaluated'), 'unwrapped_closed_now': u['closed'],
                  'seconds': round(time.time() - t0, 1), **st.rec()}); out.append(d); print(json.dumps(d), flush=True)
    open(fn, 'w').write(''.join(json.dumps(d) + '\n' for d in out))


# ------------------------------------------------------------------ CL: middle G on chosen columns (the G -> 0+ end of all 92 columns: stage CS)
def run_cl_mid(idx, fn=None):
    fn = fn or os.path.join(OUTD, 'ct_v11_cl_mid.jsonl')
    pub = [json.loads(l) for l in open(os.path.join(RES, 'cl_link_cells.jsonl')) if l.strip()]
    have = {json.loads(l)['index'] for l in open(fn) if l.strip()} if os.path.exists(fn) else set()
    orig = KL.cell_ok
    for i in idx:
        if i in have: continue
        p = pub[i]; KL.set_tau_min(p['tau_min']); st = Stat(3000 + i)
        def ok(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, Dd):
            res = orig(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, Dd)
            lPi = np.where(lmu_a >= math.log(2.0), KL.CC.lPistar_vec(np.maximum(lmu_a, math.log(2.0))), lmu_b + 2 * math.log(KL.QB)) + 1e-6   # as cl_link.cell_ok
            lGc = lPi + 0.5 * (lr_b - np.log(2 + np.exp(lr_b)))
            v0 = KL.phi_cells(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, Dd, -np.inf)
            v0p = KL.phi_cells(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, Dd, lGc - SHIFT0)
            okall = res & (v0p > ACC); st.g0p(v0, v0p, okall)
            rows = st.pick(okall)
            if len(rows):
                A = [x[rows] for x in (lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b)]
                for lo, hi in ((0.1, 700.0), (0.1, 40.0)):
                    t = st.rng.uniform(lo, hi, len(rows))
                    st.mid(KL.phi_cells(*A, ell_a, ell_b, Dd, lGc[rows] - t), np.stack(A, axis=1), t)
            return okall
        KL.cell_ok = ok; t0 = time.time()
        try:
            d = KL.certify(p['ell_a'], p['ell_b'], p['D'], p['lbmin_rel'])
        finally:
            KL.cell_ok = orig
        d.update({'index': i, 'published_cells': p['cells_evaluated'], 'seconds': round(time.time() - t0, 1), **st.rec()})
        open(fn, 'a').write(json.dumps(d) + '\n'); print(json.dumps(d), flush=True)
        if not (d['closed'] and d['cells_evaluated'] == p['cells_evaluated']): print('STOP: a column changed', flush=True); break


if __name__ == '__main__':
    os.makedirs(OUTD, exist_ok=True); m = sys.argv[1]
    if m == 'wedge': run_cq('wedge', *(int(v) for v in sys.argv[2:4]))
    elif m == 'tails_ii': run_cq('ii')
    elif m == 'windows': run_windows()
    elif m == 'cp_control': run_cp_control()
    elif m == 'cl_mid': run_cl_mid([int(v) for v in sys.argv[2].split(',')])
    elif m == 'cl':                                   # all 92 CL columns at the LIMIT G -> 0+ (stage CS evaluated G_cell e^-700 only), worker W of NW
        w, nw = int(sys.argv[2]), int(sys.argv[3])
        n = sum(1 for l in open(os.path.join(RES, 'cl_link_cells.jsonl')) if l.strip())
        fn = os.path.join(OUTD, 'ct_v11_cl_G0plus_w%d.jsonl' % w)
        have = {json.loads(l)['index'] for l in open(fn) if l.strip()} if os.path.exists(fn) else set()
        todo = [i for i in range(n) if i % nw == w and i not in have]
        if len(sys.argv) > 4: todo = todo[:int(sys.argv[4])]      # optional: at most K new columns in this run (bounded chunks)
        run_cl_mid(todo, fn)
    else: raise SystemExit(__doc__)
