"""Stage CP 2.2: does the union of boxes (+ tails + lemmas) cover the set of states the theorem talks about?  Random states, looked up in the
actual box tables of the code (initial grids and tail boxes; refined boxes are subsets of initial ones, and a column counts only if it closed).
Python 3.11.  -> results/cp_cover.json

THE STATES (Theorem H1 and the claims of the certificates; paper H v1.0 Sec. 2, 5, 6, 6a, 7):
  ell in [55.26, 960.8] (sigma = e^-ell), x in [5e11 sigma^2, 0.003], b = x/tau with tau >= 20 (b <= x/20 <= 1.5e-4), b > 0,
  mu = x S2/S0 >= x (every harmonic has n >= 1) and unbounded above, r = s^2/L^2 in [sigma^3, 1/4], A in [4 sigma^3, 1] (inside the cell bound), G in [0, G_max].
COVERS:
  'ii_v10'   certificate (ii) of v1.0 (cl_link.py, results/cl_link_cells.jsonl) + Lemma BE (0 < b <= 1e4 sigma^2)
  'ii_beta'  the beta rebuild (cp_link.initial_boxes) on the same 92 columns + Lemma BE
  'i'        certificate (i) (cd_cells.py; its claim is for every mu > 0, r in [sigma^3, 1/4], b in (0, 1], ell in [55.26, 6907.8])
  'iii'      the 17 windows (cf_windows.py): mu, x, P fixed, (r, b) in [sigma^3, 1/4] x [1e4 sigma^2, x/20] at the window's sigma
POWER: deliberately broken tables must be caught: 'ii_beta' with one column removed, 'ii_beta' with the beta edge at 2e4, 'ii_v10' itself (the wedge).
"""
import json, math, os, sys, time
import numpy as np

LN = math.log
B_MAX = 0.003 / 20


def columns(fn):
    rows = [json.loads(l) for l in open(fn) if l.strip()]
    return sorted([(r['ell_a'], r['ell_b']) for r in rows if r.get('closed')])


def find_col(cols, ell):
    for a, b in cols:
        if a - 1e-9 <= ell <= b + 1e-9: return a, b
    return None


def mu_range(ell_b, ext=False):
    mu_hi = max(4.0 * ell_b, 50.0)
    return (-30.0, LN(mu_hi) + 60) if not ext else (-math.inf, math.inf)


def in_ii_v10(s, cols, mu_ext=False):
    ell, lx, lb, lmu, lr = s
    if lb <= LN(1e4) - 2 * ell: return 'BE'
    c = find_col(cols, ell)
    if c is None: return None
    ea, eb = c; m0, m1 = mu_range(eb, mu_ext)
    if not (m0 <= lmu <= m1): return None
    if not (-3 * eb <= lr <= LN(0.25)): return None
    if not (LN(1e4) - 2 * ea <= lb <= LN(B_MAX)): return None
    return 'box'


def in_ii_beta(s, cols, lBmin=LN(1e4), mu_ext=False):
    ell, lx, lb, lmu, lr = s
    lB = lb + 2 * ell
    if lB <= LN(1e4): return 'BE'
    c = find_col(cols, ell)
    if c is None: return None
    ea, eb = c; m0, m1 = mu_range(eb, mu_ext)
    if not (m0 <= lmu <= m1): return None
    if not (-3 * eb <= lr <= LN(0.25)): return None
    if not (lBmin <= lB <= LN(B_MAX) + 2 * eb): return None
    return 'box'


def in_i(s, cols, mu_ext=False):
    ell, lx, lb, lmu, lr = s
    c = find_col(cols, ell)
    if c is None: return None
    ea, eb = c; m0, m1 = mu_range(eb, mu_ext)
    if not (m0 <= lmu <= m1): return None
    if not (-3 * eb <= lr <= LN(0.25)): return None
    if not (-700.0 * 3 - 3 * eb <= lb <= 0.0): return None
    return 'box'


def sample(n, rng, ell_lo=55.26, ell_hi=960.8, edges=()):
    out = []
    for _ in range(n):
        if edges and rng.uniform() < 0.3:
            e = float(rng.choice(edges)) + rng.uniform(-1e-3, 1e-3); e = min(max(e, ell_lo), ell_hi)
        else:
            e = rng.uniform(ell_lo, ell_hi)
        lx = rng.uniform(LN(5e11) - 2 * e, LN(0.003)) if rng.uniform() < 0.8 else LN(5e11) - 2 * e + rng.uniform(0, 1e-3)
        u = rng.uniform()
        if u < 0.3:   lb = LN(1e4) - 2 * e + rng.uniform(-3, 25)                       # near the b = 1e4 sigma^2 edge and the wedge
        elif u < 0.4: lb = min(lx - LN(20), LN(B_MAX)) - rng.uniform(0, 1e-3)            # near tau = 20 / b_max
        else:         lb = rng.uniform(LN(1e4) - 2 * e - 50, lx - LN(20))
        lb = min(lb, lx - LN(20))
        v = rng.uniform()
        if v < 0.15:  lmu = rng.uniform(lx, -30.0) if lx < -30 else lx                 # small mu, down to its true lower end x
        elif v < 0.25: lmu = LN(max(4 * e, 50)) + rng.uniform(55, 200)                   # large mu
        else:         lmu = rng.uniform(max(lx, LN(1e-4)), LN(max(4 * e, 50)))
        w = rng.uniform()
        lr = -3 * e + rng.uniform(0, 1e-3) if w < 0.1 else (LN(0.25) - rng.uniform(0, 1e-3) if w < 0.2 else rng.uniform(-3 * e, LN(0.25)))
        out.append((e, lx, lb, lmu, lr))
    return out


def run(n=100000, seed=1):
    rng = np.random.default_rng(seed); t0 = time.time()
    cols10 = columns('results/cl_link_cells.jsonl')
    edges = sorted(set([a for a, b in cols10] + [b for a, b in cols10]))
    S = sample(n, rng, edges=edges)
    res = {}
    def tally(name, f):
        miss = [s for s in S if f(s) is None]
        kinds = {'mu_below_e-30': sum(1 for s in miss if s[3] < -30), 'mu_above_tail': sum(1 for s in miss if s[3] > LN(max(4 * s[0] + 40, 50)) + 60),
                 'wedge': sum(1 for s in miss if LN(1e4) - 2 * s[0] < s[2] < LN(1e4) - 2 * find_col(cols10, s[0])[0]) if cols10 else 0}
        res[name] = {'states': len(S), 'missed': len(miss), 'by_kind': kinds, 'examples': [list(m) for m in miss[:3]]}
        print(name, res[name]['missed'], kinds, flush=True)
    tally('ii_v10 (published)', lambda s: in_ii_v10(s, cols10))
    tally('ii_v10 with mu tails extended to (0, inf)', lambda s: in_ii_v10(s, cols10, mu_ext=True))
    tally('ii_beta (92 columns as v1.0)', lambda s: in_ii_beta(s, cols10))
    tally('ii_beta with mu tails extended', lambda s: in_ii_beta(s, cols10, mu_ext=True))
    # power: broken beta tables
    cols_drop = [c for c in cols10 if not (455 < c[0] < 466)]
    tally('POWER ii_beta mu-ext, one column removed', lambda s: in_ii_beta(s, cols_drop, mu_ext=True))
    tally('POWER ii_beta mu-ext, beta edge 2e4', lambda s: in_ii_beta(s, cols10, lBmin=LN(2e4), mu_ext=True))
    # certificate (i): its own claim (b in (0, 1], any mu > 0) on its D = 10 columns
    try:
        ci = columns('results/cd_cells.jsonl' if os.path.exists('results/cd_cells.jsonl') else 'paperH/repo_cl/results/cd_cells.jsonl')
    except Exception:
        ci = []
    if ci:
        tally('i (cd_cells columns, claim b in (0,1], mu > 0)', lambda s: in_i(s, ci))
        tally('i with mu tails extended', lambda s: in_i(s, ci, mu_ext=True))
    res['seconds'] = round(time.time() - t0)
    return res


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    r = run(n); json.dump(r, open('results/cp_cover.json', 'w'), indent=1)
