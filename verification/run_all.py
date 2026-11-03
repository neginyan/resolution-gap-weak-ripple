"""Paper H (v1.0 draft) — verification suite (acceptance rules of App. E, stage CL).  Python 3.11.
Default (< 30 min): grammar; every number of the manuscript re-derived from results/ (hp_numbers.py); recomputation of a SAMPLE of every certificate
(first and last column of each cell proof, first thin-block window, the tau cells, Lemma A4* and Lemma S random checks), compared with the stored results.
--recompute: recompute EVERY column and window (hours; the cell counts must equal the stored ones).
Usage: python verification/run_all.py [--recompute] [--tex DIR]   (DIR = the folder with sections/, default paper/)."""
import sys, os, json, math, glob, ast, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, 'src')); sys.path.insert(0, HERE); os.chdir(ROOT)
FULL = '--recompute' in sys.argv
TEX = sys.argv[sys.argv.index('--tex') + 1] if '--tex' in sys.argv else 'paper'
res = []
def check(name, ok, info=''):
    res.append(ok); print(('[PASS] ' if ok else '[FAIL] ') + name + ('  ' + info if info else ''), flush=True)
def jl(f): return [json.loads(l) for l in open('results/' + f) if l.strip()]

bad = []
for f in glob.glob('src/*.py') + glob.glob('verification/*.py'):
    try: ast.parse(open(f).read(), feature_version=(3, 11))
    except SyntaxError: bad.append(f)
check('Python 3.11 grammar (%d files)' % len(glob.glob('src/*.py') + glob.glob('verification/*.py')), not bad, ', '.join(bad))

import hp_numbers
if os.path.isdir(os.path.join(TEX, 'sections')):
    rows, miss = hp_numbers.check('.', TEX)
    check('numbers of the manuscript re-derived from results/: %d numbers, %d missing' % (len(rows), len(miss)), not miss)
else:
    check('numbers: skipped (no %s/sections)' % TEX, True)

import cc_cells, cd_cells, ce_link, cf_windows, ce_split_map as SM, ce_a4, cd_split
import cl_cells, cl_link
# Stage CL: the certificates are those recomputed with the rounding-covering acceptance rules (App. E): cl_cells.cell_ok for (i),
# cl_link.phi_cells for (ii) and (iii).  The older ce_link/cc_cells records are kept only for the comparison of cell counts.
cc_cells.cell_ok = cl_cells.cell_ok
cf_windows.K = cl_link
MAIN = jl('cl_main_cells.jsonl')
CC = sorted([r for r in MAIN if r['cert'] == 'cc'], key=lambda r: r['ell_a'])
D = sorted([r for r in MAIN if r['cert'] == 'cd'], key=lambda r: r['ell_a'])
L = sorted(jl('cl_link_cells.jsonl'), key=lambda r: r['ell_a'])
W = jl('cf_windows.jsonl'); WCL = jl('cl_windows.jsonl')
for name, rows in (('main inequality, D = 0 (cc_cells, rule CL)', CC), ('main inequality, D = 10 (cd_cells, rule CL)', D), ('full Lemma R* (cl_link, rule CL)', L)):
    contig = all(abs(a['ell_b'] - b['ell_a']) < 1e-6 for a, b in zip(rows, rows[1:]))
    check('%s: %d columns, all closed, contiguous [%.2f, %.1f]' % (name, len(rows), rows[0]['ell_a'], rows[-1]['ell_b']), all(r['closed'] for r in rows) and contig)
check('certificate (i), rule CL: no audit violation; smallest accepted (lpos - lneg)/(S delta) = %.3g >= 1' % min(r['min_slack_over_S_delta'] for r in MAIN),
      all(r['audit_violations'] == 0 for r in MAIN) and min(r['min_slack_over_S_delta'] for r in MAIN) >= 1)
def recompute_cc(r): return cc_cells.certify(r['ell_a'], r['ell_b'], r['delta'])
def recompute_cd(r): return cd_cells.certify(r['ell_a'], r['ell_b'], r['delta'])
def recompute_link(r):
    cl_link.set_tau_min(r['tau_min']); return cl_link.certify(r['ell_a'], r['ell_b'], r['D'], r['lbmin_rel'])
for name, rows, fn in (('cc_cells (rule CL)', CC, recompute_cc), ('cd_cells (rule CL)', D, recompute_cd), ('cl_link (rule CL)', L, recompute_link)):
    sample = rows if FULL else [rows[0], rows[-1]]
    t0 = time.time(); same = True
    for r in sample:
        d = fn(r); same &= d['closed'] and d.get('cells_evaluated', d.get('cells')) == r['cells_evaluated']
    check('%s recomputed (%s): same cell counts' % (name, 'all %d columns' % len(sample) if FULL else 'first and last column'), same, '%.0fs' % (time.time() - t0))
ref = {48: (-17.49, 3.92), 96: (-40.88, 6.39), 300: (-142.12, 12.46)}
same = True
for w in (W if FULL else W[:1]):
    e = w['sig_exp']; x = 10 ** ref[e][0]; u = ref[e][1]
    ws = SM.window_general(x, u, w['eta_target'], w['u_block_over_u_band'] * u, w['block_width_decay_lengths'])
    S0 = ws['S0_lo'] + ws['S0_hi']; S2 = ws['S2_lo'] + ws['S2_hi']
    d = cf_windows.certify_window(e * math.log(10), math.log(x * S2 / S0), math.log(x * S2), math.log(x), w['factors'])
    same &= d['closed'] and d.get('cells_evaluated', d.get('cells')) == WCL[W.index(w)]['cells_evaluated']
check('thin-block windows (rule CL): %d stored, all closed; recomputed %s with the same cell counts' % (len(WCL), 'all' if FULL else 'the first'), len(WCL) == 17 and all(w['closed'] for w in WCL) and same)
c = ce_a4.check_A4(ntrial=300 if FULL else 60, seed=9)
check('Lemma A4* on %d random periodic windows: %d violations' % (c['trials'], c['violations']), c['violations'] == 0)
r = cd_split.check_random(ntrial=400 if FULL else 120, seed=5)
check('Lemma S on %d random gapped windows: %d violations' % (r['trials'], r['violations']), r['violations'] == 0)
# ---- Lemma BE (stage CH), the narrow focus and Lemma A* (stages CI, CJ)
import subprocess
def runpy(script):
    t0 = time.time(); r = subprocess.run([sys.executable, '-W', 'ignore', 'src/' + script], capture_output=True, text=True, timeout=7200); return r.returncode == 0, time.time() - t0
ok, dt = runpy('ch_bedge.py'); B = json.load(open('results/ch_bedge.json'))
check('Lemma BE: formula of sigma^2 G, mean-square bound, constants (kappa = %.3e) and every remainder (< 10^%.3g) re-checked' % (B['BE_3_constants']['kappa'], B['BE_3_remainders_log10']['all_below_log10']),
      ok and B['all_ok'] and B['BE_3_constants']['kappa'] >= 1.1e-5, '%.0fs' % dt)
ok, dt = runpy('ci_rtheta.py'); t = json.load(open('results/ci_rtheta.json'))
check('local corrections (L1)-(L3) on %d random windows' % t['trials'], ok and t['ok'], '%.0fs' % dt)
ok, dt = runpy('ci_local.py'); lo = json.load(open('results/ci_local.json'))
check('local Z3\' on (P2) windows: max ratio %.3f (paper: <= 2.13)' % lo['worst']['rho_max'], ok and lo['worst']['rho_max'] <= 2.13, '%.0fs' % dt)
ok, dt = runpy('ci_envloc.py'); ev = json.load(open('results/ci_envloc.json'))
check('local envelope with brightness excess C = 10 stays above T_10 at %d values of ell' % len(ev), ok and all(v['local_C10'] >= v['T10'] for v in ev.values()), '%.0fs' % dt)
ok, dt = runpy('ci_nop2.py'); npj = json.load(open('results/cj_nop2.json'))
check('windows without (P2) exceed Pi*(mu): ratios %.2f to %.0f' % (npj['min_ratio'], npj['max_ratio']), ok and npj['min_ratio'] > 3.4, '%.0fs' % dt)
ok, dt = runpy('cl_tau_be.py'); TB = json.load(open('results/cl_tau_be.json'))
check('tau cells and Lemma BE against rounding: tau-factor margin %.1e = %.1e x the derived bound; BE in 60-digit decimal' % (TB['tau_cells']['margin'], TB['tau_cells']['margin_over_bound']), ok and TB['all_ok'], '%.0fs' % dt)
ok, dt = runpy('ck_details.py'); DT = json.load(open('results/ck_details.json'))
check('App. D line by line (Lemma S, S\', A4*, Z3\', BE) in exact fractions and 60-digit decimal', ok and DT['all_ok'], '%.0fs' % dt)
import ci_link
P = jl('ci_link_cells.jsonl'); sample = P if FULL else [p for p in P if p['ell_a'] == 55.26 and p['tau_min'] == 12.0]
same = True; t0 = time.time()
for p in sample:
    ci_link.set_tau_min(p['tau_min']); d = ci_link.certify(p['ell_a'], p['ell_b'], p['D'], p['lbmin_rel'], max_cells=p.get('cells_evaluated', p.get('cells', 0)) + 10)
    same &= d['closed'] == p['closed'] and (not p['closed'] or d['cells_evaluated'] == p['cells_evaluated'])
check('single columns with Lemma A* recomputed (%s): same outcome and cell counts' % ('all %d' % len(sample) if FULL else 'ell [55.26, 65.26], tau >= 12'), same, '%.0fs' % (time.time() - t0))
print('\n%d/%d checks passed' % (sum(res), len(res))); sys.exit(0 if all(res) else 1)
