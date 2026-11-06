"""Stage CQ 2.4: the covering check of cp_cover.py on the REPAIRED tables.  Python 3.11.  -> results/cq_cover.json
(ii): v1.0 b boxes (results/cl_link_cells.jsonl: mu in [e^-30, mu_hi e^60], r in [sigma_bot^3, 1/4], b in [1e4 sigma_top^2, b_max])
      + mu tails (results/cq_tails_ii.jsonl: mu in (0, e^-30] and [mu_hi, inf), same r and b)
      + wedge (results/cq_wedge.jsonl: all mu, same r, beta in [1e4, 1e4 e^{2(ell_b - ell_a)}], ell in the column)
      + Lemma BE (0 < b <= 1e4 sigma^2).  A column counts only if its result line says closed.
(i):  for each of its two claims (D = 0 on [55.26, 3062.6], D = 10 on [55.26, 6907.8]): v1.0 boxes (results/cl_main_cells.jsonl: mu in [1e-4, mu_hi] with
      b in [e^{-2100-3 ell_b}, 1], mu in [e^-30, 1e-4] and [mu_hi, mu_hi e^60] with b in [sigma_bot^3, 1]; r in [sigma_bot^3, 1/4]) + tails
      (results/cq_tails_i_w*.jsonl: mu tails with every b in (0, 1]; b below the v1.0 b ends for the v1.0 mu range).
Power: deliberately broken tables (a wedge column removed, the wedge beta top lowered by e, the upper mu tail removed, the v1.0 table alone; for (i)
the b -> 0 strip removed and a column removed) must be caught."""
import json, math, glob, sys, time
import numpy as np
import cp_cover as V

LN = math.log; B_MAX = 0.003 / 20


def closed_cols(fns, key=None):
    rows = []
    for fn in fns:
        try: rows += [json.loads(l) for l in open(fn) if l.strip()]
        except FileNotFoundError: pass
    if key: rows = [r for r in rows if key(r)]
    return sorted((r['ell_a'], r['ell_b']) for r in rows if r.get('closed'))


def in_ii(s, body, tails, wedge, wedge_top_shift=0.0, no_upper=False):
    ell, lx, lb, lmu, lr = s
    if lb <= LN(1e4) - 2 * ell: return 'BE'
    lB = lb + 2 * ell
    for cols, kind in ((body, 'body'), (tails, 'tails'), (wedge, 'wedge')):
        c = V.find_col(cols, ell)
        if c is None: continue
        ea, eb = c
        if not (-3 * eb <= lr <= LN(0.25)): continue
        mh = LN(max(4.0 * eb, 50.0))
        if kind == 'body':
            if -30.0 <= lmu <= mh + 60 and LN(1e4) - 2 * ea <= lb <= LN(B_MAX): return 'body'
        elif kind == 'tails':
            if (lmu <= -30.0 or (lmu >= mh and not no_upper)) and LN(1e4) - 2 * ea <= lb <= LN(B_MAX): return 'tails'
        else:
            if LN(1e4) <= lB <= LN(1e4) + 2 * (eb - ea) - wedge_top_shift and (not no_upper or lmu < mh): return 'wedge'
    return None


def cols_with(cols, ell):
    return [(a, b) for a, b in cols if a - 1e-9 <= ell <= b + 1e-9]


def in_i(s, body, tails, no_strip=False):
    ell, lb, lmu, lr = s
    for cols, kind in ((body, 'body'), (tails, 'tails')):
        for ea, eb in cols_with(cols, ell):
            if not (-3 * eb <= lr <= LN(0.25)) or lb > 0: continue
            mh = LN(max(4.0 * eb, 50.0)); strip_top = -2100.0 - 3 * eb
            if kind == 'body':      # v1.0: 16 mu boxes [1e-4, mu_hi] with b in [e^{-2100-3 ell_b}, 1]; mu tails [e^-30, 1e-4], [mu_hi, mu_hi e^60] with b in [sigma_bot^3, 1]
                if LN(1e-4) <= lmu <= mh and lb >= strip_top: return 'body'
                if (-30.0 <= lmu <= LN(1e-4) or mh <= lmu <= mh + 60) and lb >= -3 * eb: return 'body'
            else:
                if lmu <= -30.0 or lmu >= mh: return 'tails'
                if not no_strip and LN(1e-4) <= lmu <= mh and lb <= strip_top: return 'strip'
                if not no_strip and (-30.0 <= lmu <= LN(1e-4) or mh <= lmu <= mh + 60) and lb <= -3 * eb: return 'strip'
    return None


def sample_i(n, rng, ell_hi, edges):
    out = []
    for _ in range(n):
        e = (float(rng.choice(edges)) + rng.uniform(-1e-3, 1e-3)) if rng.uniform() < 0.3 else rng.uniform(55.26, ell_hi)
        e = min(max(e, 55.26), ell_hi)
        u = rng.uniform()
        lb = rng.uniform(-3 * e, 0) if u < 0.4 else (rng.uniform(-2100 - 3 * e - 5, -2100 - 3 * e + 5) if u < 0.55 else (rng.uniform(-1e5, -2100 - 3 * e) if u < 0.7 else rng.uniform(-2100 - 3 * e, -3 * e)))
        v = rng.uniform()
        lmu = rng.uniform(-4000, -30) if v < 0.2 else (rng.uniform(-31, -29) if v < 0.3 else (LN(max(4 * e, 50)) + rng.uniform(-1, 400) if v < 0.5 else rng.uniform(LN(1e-4), LN(max(4 * e, 50)))))
        w = rng.uniform()
        lr = -3 * e + rng.uniform(0, 1e-3) if w < 0.1 else (LN(0.25) - rng.uniform(0, 1e-3) if w < 0.2 else rng.uniform(-3 * e, LN(0.25)))
        out.append((e, lb, lmu, lr))
    return out


def run(n=100000, seed=1):
    t0 = time.time(); rng = np.random.default_rng(seed); res = {}
    body = V.columns('results/cl_link_cells.jsonl'); tails = closed_cols(['results/cq_tails_ii.jsonl']); wedge = closed_cols(['results/cq_wedge.jsonl'])
    edges = sorted(set([a for a, b in body] + [b for a, b in body]))
    S = V.sample(n, rng, edges=edges)
    def tally(name, f, SS=S):
        miss = [s for s in SS if f(s) is None]
        res[name] = {'states': len(SS), 'missed': len(miss), 'examples': [list(m) for m in miss[:3]]}; print(name, len(miss), flush=True)
    res['columns'] = {'ii_body': len(body), 'ii_tails_closed': len(tails), 'ii_wedge_closed': len(wedge)}
    tally('ii repaired (body + tails + wedge + BE)', lambda s: in_ii(s, body, tails, wedge))
    tally('POWER ii: v1.0 table alone', lambda s: in_ii(s, body, [], []))
    tally('POWER ii: one wedge column removed', lambda s: in_ii(s, body, tails, [c for c in wedge if not (455 < c[0] < 466)]))
    tally('POWER ii: wedge beta top lowered by e', lambda s: in_ii(s, body, tails, wedge, wedge_top_shift=1.0))
    tally('POWER ii: upper mu tail removed', lambda s: in_ii(s, body, tails, wedge, no_upper=True))
    # (i)
    tl_rows = []
    for fn in glob.glob('results/cq_tails_i_w*.jsonl'): tl_rows += [json.loads(l) for l in open(fn) if l.strip()]
    main = [json.loads(l) for l in open('results/cl_main_cells.jsonl') if l.strip()]
    for kind, hi in (('cc', 3062.6), ('cd', 6907.8)):
        # the D = 10 columns also prove the D = 0 claim (T_10 >= T_0): the cc columns start at 55.262, and [55.26, 55.262) is covered by cd
        kinds = ('cc', 'cd') if kind == 'cc' else ('cd',)
        bcols = sorted((r['ell_a'], r['ell_b']) for r in main if r['cert'] in kinds and r['closed'])
        tcols = sorted((r['ell_a'], r['ell_b']) for r in tl_rows if r['cert'] in kinds and r['closed'])
        res['columns']['i_%s_body' % kind] = len(bcols); res['columns']['i_%s_tails_closed' % kind] = len(tcols)
        ed = sorted(set([a for a, b in bcols] + [b for a, b in bcols])); SI = sample_i(n // 2, rng, hi, ed)
        tally('i %s repaired (body + mu tails + b strip)' % kind, lambda s: in_i(s, bcols, tcols), SI)
        tally('POWER i %s: v1.0 table alone' % kind, lambda s: in_i(s, bcols, []), SI)
        tally('POWER i %s: b -> 0 strip removed' % kind, lambda s: in_i(s, bcols, tcols, no_strip=True), SI)
        tally('POWER i %s: one column removed' % kind, lambda s: in_i(s, [c for c in bcols if not (455 < c[0] < 466)], tcols), SI)
        if kind == 'cc':
            only = sorted((r['ell_a'], r['ell_b']) for r in main if r['cert'] == 'cc' and r['closed'])
            tally('POWER i cc: the cc columns alone (start at 55.262)', lambda s: in_i(s, only, tcols), SI)
    res['seconds'] = round(time.time() - t0)
    return res


if __name__ == '__main__':
    r = run(int(sys.argv[1]) if len(sys.argv) > 1 else 100000); json.dump(r, open('results/cq_cover.json', 'w'), indent=1)
