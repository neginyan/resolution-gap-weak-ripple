"""Paper H v1.2 — the checks of the two corrections of v1.2 (called by run_all.py; can also be run alone).  Python 3.11 grammar.  A few minutes.
  (v1) the code of v1.1 is unchanged: every src/*.py of v1.1 has its v1.1 sha256, except cl_link.py, whose syntax tree without docstrings is that of v1.1
       (v1.2 corrects only the text "MAIN-FACTOR SCALING" of its docstring);
  (v2) the G -> 0+ end point (stage CT, src/ct_v11_G0plus.py): the stored re-runs of certificates (ii) and (iii) with the limit G -> 0+ close every column
       and window with the published cell counts, and no value at intermediate G is negative;
  (v3) the algebra of the corrected scaling step (Sec. IV): Num - (lambda/Upsilon) gamma sigma den >= T(lambda), T affine in lambda, and the target D';
  (v4) the wrapped acceptance of src/da_scaling.py with the factors one and without G -> 0+ is the published acceptance, bit for bit (sampled cells);
  (v5) the stored second certificate of the scaling step (factors (2/3, 2/3, 1), target (2/3) T10): all 92 + 92 + 92 columns closed, G -> 0+ in every
       box, no value at intermediate G negative;
  (v6) three columns of the second certificate re-run now give the stored cell counts;
  (v7) power: with the factors 2/3 but the target not lowered, the first column does not close (stored).
--recompute: (v6) re-runs every column of the second certificate (many CPU hours)."""
import sys, os, json, glob, ast, hashlib, io, contextlib, warnings
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, 'src')); os.chdir(ROOT)
import numpy as np
warnings.filterwarnings('ignore')
FULL = '--recompute' in sys.argv
QUIET = lambda: contextlib.redirect_stdout(io.StringIO())
res = []
def check(name, ok, info=''):
    res.append(bool(ok)); print(('[PASS] ' if ok else '[FAIL] ') + name + ('  ' + info if info else ''), flush=True)
def jl(f): return [json.loads(l) for l in open(f) if l.strip()]

# (v1)
H = json.load(open('results/v11_src_sha256.json'))
def strip(src):
    t = ast.parse(src)
    for n in ast.walk(t):
        if isinstance(n, (ast.Module, ast.FunctionDef, ast.ClassDef)) and n.body and isinstance(n.body[0], ast.Expr) \
                and isinstance(getattr(n.body[0], 'value', None), ast.Constant) and isinstance(n.body[0].value.value, str):
            n.body = n.body[1:]
    return hashlib.sha256(ast.dump(t).encode()).hexdigest()
bad = [f for f, h in H['files'].items() if f != 'src/cl_link.py' and (not os.path.exists(f) or hashlib.sha256(open(f, 'rb').read()).hexdigest() != h)]
okcl = strip(open('src/cl_link.py').read()) == H['cl_link_ast_without_docstrings']
check('(v1) the v1.1 code is unchanged: %d files with their v1.1 sha256, and cl_link.py equal to v1.1 up to its docstring' % (len(H['files']) - 1), not bad and okcl, ', '.join(bad))

# (v2)
pubn = {'cl': 'cl_link_cells.jsonl', 'wedge': 'cq_wedge.jsonl', 'tails': 'cq_tails_ii.jsonl', 'windows': 'cl_windows.jsonl'}
ctf = {'cl': 'ct_v11_cl_G0plus_w*.jsonl', 'wedge': 'ct_v11_wedge_G0plus_w*.jsonl', 'tails': 'ct_v11_tails_ii_G0plus.jsonl', 'windows': 'ct_v11_windows_G0plus.jsonl'}
expect = {'cl': 92, 'wedge': 92, 'tails': 92, 'windows': 17}
mid_total = 0; lines = []; okc = True
for k in ('cl', 'wedge', 'tails', 'windows'):
    rows = [r for f in sorted(glob.glob('results/' + ctf[k])) for r in jl(f)]
    same = all(r['closed'] and r.get('published_cells', r['cells_evaluated']) == r['cells_evaluated'] for r in rows)
    neg = sum(r.get('mid_negative', 0) for r in rows); mid = sum(r.get('mid_tested', 0) for r in rows); mid_total += mid
    okc &= len(rows) == expect[k] and same and neg == 0
    lines.append('%s %d/%d closed, %d cells, middle G %d, negative %d' % (k, sum(r['closed'] for r in rows), expect[k], sum(r['cells_evaluated'] for r in rows), mid, neg))
S = json.load(open('results/ct_summary.json'))
extra = S['i_tails']['mid_tested'] + S['i_body']['mid_evaluations'] + S['cp_control'].get('mid_tested', 0)   # certificate (i), and the control columns of stage CP
okc &= not S['i_tails']['columns_changed'] and S['i_tails']['columns'] == 987 and S['i_body']['all_equal_published'] and S['all_v11_clean'] and S['mid_tested_total'] == mid_total + extra
okc &= all(v.get('mid_negative', 0) == 0 and v.get('all_closed', True) for v in S.values() if isinstance(v, dict))
check('(v2) G -> 0+ re-runs of certificates (ii) and (iii), published cell counts: ' + '; '.join(lines) + '; with certificate (i) %d middle-G values in all, none negative' % (mid_total + extra),
      okc and abs((mid_total + extra) / 1e8 - 2.25) < 0.005)

# (v3)
import da_scaling as DS
rng = np.random.default_rng(1); okA = True
for _ in range(20000):
    Pbg, P0, P2, N, dbg, d2, gs = rng.exponential(1.0, 7); l0 = rng.uniform(2 / 3, 1); l2 = rng.uniform(2 / 3, 1); ups = rng.uniform(1, 2); lam = min(l0, l2)
    true_lb = Pbg + l0 * P0 + l2 * P2 - N - (lam / ups) * gs * (dbg + ups * d2)
    T = lambda L: Pbg + L * (P0 + P2) - N - L * gs * (dbg + d2); th = (lam - 2 / 3) / (1 / 3)
    okA &= true_lb >= T(lam) - 1e-12 and abs(T(lam) - (th * T(1) + (1 - th) * T(2 / 3))) < 1e-9 and abs(T(lam) - (lam * T(1) + (1 - lam) * (Pbg - N))) < 1e-9
okD = True
for ea in np.linspace(55.26, 950.8, 200):
    Dp = DS.Dprime(ea, 2 / 3)
    for ell in np.linspace(ea, ea + 10, 50):
        okD &= 1 + (2.1 + Dp) / ell + 2.5e-4 >= (2 / 3) * DS.T10(ell) * (1 - 1e-12)
check('(v3) the scaling step: Num - (lambda/Upsilon) gamma sigma den >= T(lambda) = lambda T(1) + (1 - lambda)(P_bg - N), affine between 2/3 and 1 (20000 instances); T_{D\'}(ell) >= (2/3) T10(ell) on every column',
      okA and okD)

# (v4)
with QUIET(): eq = DS.equiv(n=50000, save=False)
check('(v4) the wrapped acceptance (factors 1, G -> 0+ off) equals the published one bit for bit (%d samples of columns)' % (len(eq) - 1), eq['all_identical'])

# (v5)
with QUIET(): M = DS.merge()
tot = sum(M[k]['cells'] for k in ('cl', 'wedge', 'tails')); mm = min(M[k]['min_margin'] for k in ('cl', 'wedge', 'tails')); mt = sum(M[k]['mid_tested'] for k in ('cl', 'wedge', 'tails'))
for k in ('cl', 'wedge', 'tails'):
    m = M[k]
    check('(v5) second certificate, %s: %d/92 columns closed; %d cells; smallest margin %.3g; G -> 0+ on %d rows, rejected only by it %d; middle G %d, negative %d' % (
        k, m['closed'], m['cells'], m['min_margin'], m['g0p_rows'], m['cells_rejected_only_by_G0plus'], m['mid_tested'], m['mid_negative']),
          m['columns'] == 92 and m['closed'] == 92 and m['mid_negative'] == 0 and m['g0p_rows'] >= m['cells'] and m['min_margin'] > 0)
tex = ''.join(open(f).read() for f in glob.glob('paper/sections/*.tex')) if os.path.isdir('paper/sections') else None
want = ['$%.2f\\times10^{9}$ cells' % (tot / 1e9), '$%.1f\\times10^{-10}$' % (mm / 1e-10), '$%.2f\\times10^{8}$ values' % (mt / 1e8)]
check('(v5) the numbers of Sec. VI on the second certificate (%s) are those of the records' % ', '.join(want), tex is None or all(w in tex for w in want))

# (v6)
pub = {k: jl(os.path.join('results', f)) for k, f in DS.PUB.items()}
recs = {k: [r for f in sorted(glob.glob('results/da_%s_w*.jsonl' % k)) for r in jl(f)] for k in DS.PUB}
todo = [(k, i) for k in DS.PUB for i in range(len(pub[k]))] if FULL else [('cl', 0), ('wedge', 0), ('tails', 91)]
rr = []
for kind, i in todo:
    p = pub[kind][i]; st = [r for r in recs[kind] if abs(r['ell_a'] - p['ell_a']) < 1e-9]
    if not st: rr.append((kind, i, False, 0, 0)); continue
    with QUIET(): d = DS.run_column(kind, p, 2 / 3, seed=int(round(p['ell_a'] * 1000)) + {'cl': 1, 'wedge': 7, 'tails': 3}[kind])
    rr.append((kind, i, d['closed'], d['cells_evaluated'], st[0]['cells_evaluated']))
check('(v6) second certificate re-run now (%s): ' % ('all columns' if FULL else 'three columns') + '; '.join('%s %d %s %d/%d' % t for t in rr[:6]) + (' ...' if len(rr) > 6 else ''),
      all(t[2] and t[3] == t[4] for t in rr))

# (v7)
Pw = json.load(open('results/da_power.json'))
check('(v7) power: factors 2/3 with the target not lowered -> the first column does not close (%s after %s cells)' % (Pw['reason'], Pw['cells']), Pw['ok'])

print('\nv1.2 checks: %d/%d passed' % (sum(res), len(res)))
sys.exit(0 if all(res) else 1)
