"""Paper H v1.3 -- the checks of v1.3 (called by run_all.py; can also be run alone).  Python 3.11 grammar.  About ten minutes.
v1.3 changes four files of the v1.1/v1.2 code (results/v13_code_changes.json) and adds src/dv_m2.py, src/dv_interval.py, src/dv_section9.py.
  (w1) undoing the documented v1.3 changes gives back the v1.1 sha256 of the four files; no other change is made to them;
  (w2) verification/run_v12.py, unchanged, passes on a temporary copy of the repository in which the four files are undone (all v1.2 checks:
       the G -> 0+ re-runs, the second certificate of the scaling step, three columns re-run, power) -- so v1.3 does not touch Theorem H1;
  (w3) the constant M2 of Lemma R_theta (N16): src/dv_m2.py re-run now (hand bound for y <= 1e-20, interval arithmetic with the decimal.Inexact
       trap, power checks) equals its stored record; M2_const of co_remainders.py and co_rtheta.py is dv_m2.M2_bound at y rounded up;
  (w4) the range of y = (1 - b)^2 sigma^2/(2b): every state of co_remainders.json has y <= 1e-20 (hand bound); the states of co_rtheta.json have
       y in [8e-4, 0.31]; Lemma R_theta recomputed with the new M2 on a sample of the stored states (all of them with --recompute) equals the record;
  (w5) cn_lemmas.py uses sqrt2/(16e) (the value of Lemma CM-1) and its stored checks pass;
  (w6) co_rtheta.py: default argument 40, the record has 46 states (40 + 6 packets), Lemma R_theta below quadrature on all;
  (w7) the top end ell = 960.8 in binary64: the six stored runs of src/top_end.py (lambda = 1 and 2/3; Lemma R* boxes, wedge, mu tails) have ell_b =
       the next binary64 number above 960.8, are closed, and have the cell counts of the last column with the old end (published certificate, or second certificate for 2/3);
  (w8) every number printed in Sec. IX A is the same in v1.2 and v1.3 (results/v12_section9_values.json against results/ now);
  (w9) the text of v1.3 (version paragraph, App. D lemma on M2, the binary64 sentence, the setting of CM-1) is in paper/sections; the header of
       cm_local_bound.py no longer says b <= 1;
  (w10) power: undoing only some of the changes does not give the v1.1 sha256; the v1.2 grid value of M2 is below the v1.3 bound (the bound is not
       the grid value); a v1.3 text check with a sentence removed fails.
  (w11) Lemma BE at the upper end sigma = e^{-55.26} (src/dw_bedge_low.py, stage DW) and the new range in the text.
--recompute: (w4) on every stored state; run_v12.py with --recompute."""
import sys, os, json, glob, math, hashlib, shutil, subprocess, tempfile, inspect, io, contextlib, warnings, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, 'src')); os.chdir(ROOT)
warnings.filterwarnings('ignore')
FULL = '--recompute' in sys.argv
QUIET = lambda: contextlib.redirect_stdout(io.StringIO())
res = []
def check(name, ok, info=''):
    res.append(bool(ok)); print(('[PASS] ' if ok else '[FAIL] ') + name + ('  ' + info if info else ''), flush=True)
sha = lambda s: hashlib.sha256(s.encode() if isinstance(s, str) else s).hexdigest()

# (w1)
H = json.load(open('results/v11_src_sha256.json'))['files']; CH = json.load(open('results/v13_code_changes.json'))['changes']
FILES = sorted({c[0] for c in CH})
def undo(f, changes):
    s = open(f).read()
    for g, new, old in changes:
        if g == f:
            if s.count(new) != 1: return None
            s = s.replace(new, old)
    return s
und = {f: undo(f, CH) for f in FILES}
ok1 = all(und[f] is not None and sha(und[f]) == H[f] for f in FILES)
check('(w1) undoing the %d documented v1.3 changes gives the v1.1 sha256 of %s' % (len(CH), ', '.join(FILES)), ok1)

# (w2)
t0 = time.time(); tmp = tempfile.mkdtemp(prefix='h13_v12_')
try:
    for d in ('src', 'results', 'verification', 'paper'):
        if os.path.isdir(d): shutil.copytree(d, os.path.join(tmp, d), ignore=shutil.ignore_patterns('__pycache__', '*.pdf'))
    for f in FILES: open(os.path.join(tmp, f), 'w').write(und[f] or '')
    r = subprocess.run([sys.executable, '-W', 'ignore', os.path.join(tmp, 'verification', 'run_v12.py')] + (['--recompute'] if FULL else []),
                       capture_output=True, text=True, timeout=36000 if not FULL else 10 ** 6)
    print(''.join('    ' + l + '\n' for l in r.stdout.splitlines() if l.strip()), end='')
finally:
    shutil.rmtree(tmp, ignore_errors=True)
check('(w2) verification/run_v12.py (unchanged) passes with the four files undone: all v1.2 checks of Theorem H1', r.returncode == 0 and 'v1.2 checks:' in r.stdout, '%.0fs' % (time.time() - t0))

# (w3)
import dv_m2, co_remainders as COR, co_rtheta as COT
from ca_dec import D
t0 = time.time()
with QUIET(): m = dv_m2.run(save=False)
st = json.load(open('results/dv_m2.json'))
same = m['ok'] and st['ok'] and all(abs(m['(b) interval bound vs grid'][k]['bound'] - v['bound']) == 0 for k, v in st['(b) interval bound vs grid'].items())
yv = [1e-30, 1e-21, 0.01, 0.2]; okc = True
for y in yv:
    tau = D('2.5'); beta = (tau * D(repr(y))).sqrt(); b1, b2 = COR.M2_const(beta, tau), COT.M2_const(beta, tau)
    okc &= b1 == b2 == dv_m2.M2_bound(math.nextafter(float(beta * beta / tau) * (1 + 1e-12), math.inf)) and b1 >= dv_m2.grid(y)
check('(w3) M2 of Lemma R_theta: dv_m2.py re-run now equals its record (hand bound 2.737605 for y <= 1e-20, interval bound above, power checks %s); '
      'M2_const of co_remainders.py and co_rtheta.py is that bound' % ('all broken as required' if all(m['power'].values()) else 'NOT all'), same and okc, '%.0fs' % (time.time() - t0))

# (w4)
CO = json.load(open('results/co_remainders.json')); ys = []
for key, row in CO.items():
    ell = {'55.26': 55.26204223, '110.52': 110.52408446, '221.05': 221.04816893, '455.26': 455.26, '690.78': 690.77552790, '960.8': 960.8}[key]
    for name, v in row.items():
        if isinstance(v, dict) and 'state' in v:
            lb = v['state']['log_b']; ys.append((2 * math.log1p(-math.exp(lb)) - math.log(2) - lb - 2 * ell) / math.log(10))   # log10 y, y = (1-b)^2 sigma^2/(2b)
RT = json.load(open('results/co_rtheta.json')); yr = [(1 - s['b']) ** 2 * s['sigma'] ** 2 / (2 * s['b']) for s in RT['states']]
import ca_reduction as CA
idx = list(range(len(RT['states']))) if FULL else sorted({0, 7, 20, 40, 45, min(range(len(RT['states'])), key=lambda i: RT['states'][i]['ratio_theta'] if RT['states'][i]['lemma_R_period'] < 0 else 9)})
t0 = time.time(); okr = True; worst = 0
for i in idx:
    s = RT['states'][i]; wn = CA.window_numbers(s['q0'], s['phases'], s['x'])
    E0, E1, kap, Rth, E1_2 = COT.local_numbers(s['q0'], s['phases'], s['x'], s['tau'])
    lr = COT.lemma_R_theta(s['sigma'], s['x'], s['b'], s['L2'], s['s2'], wn['eta'], wn['Qb'], wn['N1'], wn['N2'], E0, E1, kap, Rth, E1_2)
    d = abs(float(lr['lower']) / s['lemma_R_theta'] - 1); worst = max(worst, d); okr &= d < 1e-9
check('(w4) y on the states of Sec. IX A: co_remainders y <= 10^%.2f (<= 1e-20, hand bound), co_rtheta y in [%.2g, %.3g]; Lemma R_theta recomputed with the '
      'new M2 on %d stored states equals the record (largest relative difference %.1e)' % (max(ys), min(yr), max(yr), len(idx), worst),
      max(ys) <= -20 and 8e-4 <= min(yr) and max(yr) <= 0.31 and okr, '%.0fs' % (time.time() - t0))

# (w5)
import cn_lemmas as CL
src = open('src/cn_lemmas.py').read(); P, mu, b, kap, K, L = 1e-3, 5.0, 0.3, 0.2, 4, 2
E0 = P / (2 * mu); bp = b / (L * K * kap); A = math.sqrt(L * K * kap) * math.exp(1 / (2 * L)) / (1 - 1 / K); Eh = min(E0 / math.sqrt(L * K * kap), 1 / 16)
manual = A * ((1 + bp / 2) * Eh * math.log(2 * CL.B0(Eh, bp) / Eh) + math.sqrt(2) / (16 * math.e) * bp) / E0
okcl = abs(CL.rhs2(P, mu, b, kap, K, L) / manual - 1) < 1e-15 and '0.0325 * bp' not in src and all(v.get('ok', True) for v in json.load(open('results/cn_lemmas.json')).values() if isinstance(v, dict))
check('(w5) cn_lemmas.py uses sqrt2/(16e) = %.10f (Lemma CM-1), not 0.0325; its stored checks pass' % (math.sqrt(2) / (16 * math.e)), okcl)

# (w6)
dflt = inspect.signature(COT.run).parameters['nstates'].default
check('(w6) co_rtheta.py: default argument %d; the record has %d states, all below quadrature' % (dflt, len(RT['states'])), dflt == 40 and len(RT['states']) == 46 and RT['all_ok'])

# (w7)
top = [json.loads(l) for f in ('results/top_end_lam1.jsonl', 'results/top_end_lam2_3.jsonl') for l in open(f) if l.strip()]
okt = len(top) == 6 and {(t['kind'], t['lam']) for t in top} == {(k, l) for k in ('tails', 'cl', 'wedge') for l in (1.0, 2 / 3)} \
      and all(t['ell_b'] == math.nextafter(960.8, math.inf) and t['closed'] and t.get('mid_negative', 0) == 0 for t in top)
last2 = {}
for k in ('cl', 'wedge', 'tails'):
    rows = [json.loads(l) for f in sorted(glob.glob('results/da_%s_w*.jsonl' % k)) for l in open(f) if l.strip()]
    r = max(rows, key=lambda r: r['ell_a']); last2[k] = (r['ell_b'], r['cells_evaluated'])
okt &= all((t['cells_minus_published'] == 0) if t['lam'] == 1.0 else (last2[t['kind']][0] == 960.8 and t['cells_evaluated'] == last2[t['kind']][1]) for t in top)
check('(w7) top end 960.8 in binary64: six runs with ell_b = %r, all closed, with the cell counts of the last column with the old end (published certificate for '
      'lambda = 1, second certificate for lambda = 2/3; %s cells)' % (math.nextafter(960.8, math.inf), '{:,}'.format(sum(t['cells_evaluated'] for t in top))), okt)

# (w8)
import dv_section9 as S9
old = json.load(open('results/v12_section9_values.json')); new = S9.values('results'); diffs = []
okp = set(old) == set(new)
for k in old:
    a, f = old[k]; c = new[k][0]
    okp &= S9.printed(a, f) == S9.printed(c, f); diffs.append(abs(c / a - 1) if a else abs(c))
check('(w8) every number printed in Sec. IX A is unchanged from v1.2 (%d values; largest relative change of the raw values %.1e)' % (len(old), max(diffs)), okp)

# (w9)
tex = ''.join(open(f).read() for f in sorted(glob.glob('paper/sections/*.tex')))
WANT = ['\\paragraph*{Version~1.3.}', 'This version (v1.3)', '\\label{app:details:M2}', '$M_2(y)\\le2.737605$', 'next binary64 number',
        'whatever $b=x/\\tau>0$', '($40$ random states and $6$ wave packets']
okx = all(w in tex for w in WANT) and '(t0, tau = x/b, b <= 1)' not in open('src/cm_local_bound.py').read()
check('(w9) the v1.3 text is in paper/sections (%d sentences) and the header of cm_local_bound.py is corrected' % len(WANT), okx)

# (w10)
part = [undo(f, CH[:-1]) for f in FILES]
p1 = any(p is None or sha(p) != H[f] for p, f in zip(part, FILES))
p2 = dv_m2.grid(0.2) * (1 + 1e-6) < dv_m2.M2_bound(0.2)
p3 = not all(w in tex.replace(WANT[3], '') for w in WANT)
check('(w10) power: undoing all but one change misses the v1.1 sha256 (%s); the v1.2 grid value is below the v1.3 bound at y = 0.2 (%s); a removed sentence is caught (%s)' % (p1, p2, p3), p1 and p2 and p3)

# (w11) Lemma BE at the upper end sigma = e^{-55.26} (stage DW)
import dw_bedge_low as BL
with QUIET(): bl = BL.run(save=False)
st = json.load(open('results/dw_bedge_low.json'))
same = all(bl[k] == st[k] for k in st if k != 'seconds')
txt = ('Let $\\sigma\\in[10^{-3000},e^{-55.26}]$' in tex and 'for $\\sigma\\le e^{-55.26}$ and $U_1\\le29.4$' in tex and '\\code{dw\\_bedge\\_low.py}' in tex
       and 'Let $\\sigma\\in[10^{-3000},10^{-24}]$' not in tex)
check('(w11) Lemma BE for sigma <= e^{-55.26}: dw_bedge_low.py re-run equals its record (rho_W <= 0.0801, kappa >= 1.1e-5, the sigma^2 term %.1e < 1e-40, '
      'tail ratio 10^%.3g, kappa sqrt(ell)/sigma = %.3g >= 8e19; power: %d wrong versions broken); the text states the new range' % (
      bl['(2) sigma^2 term of rho_W at sigma_up'], bl['(3) log10 of the tail ratio at sigma_up'], bl['(4) kappa sqrt(ell)/sigma at ell = 55.26 (without 1 - epsbar)'],
      sum(bl['power (each wrong version broken)'].values())), bl['ok'] and same and txt)

print('\nv1.3 checks: %d/%d passed' % (sum(res), len(res)))
sys.exit(0 if all(res) else 1)
