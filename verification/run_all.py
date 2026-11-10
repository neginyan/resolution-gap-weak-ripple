"""Paper H v1.2 — verification suite (acceptance rules of App. E; the covering of the states, Sec. 6).  Python 3.11.
Default (< 30 min): grammar; every number of the manuscript re-derived from results/ (hp_numbers.py); recomputation of a SAMPLE of every certificate
(first and last column of each cell proof, first thin-block window, the tau cells, Lemma A4* and Lemma S random checks), compared with the stored results.
v1.1 adds: the random covering check (1e5 states for (ii), 5e4 for each claim of (i)), the generated tail/wedge functions, sample columns of the mu tails
(i) and (ii) and of the wedge, the checks of the new steps (reduced scale), and the stored narrow-focus results of Sec. IX A.
v1.2 adds verification/run_v12.py: the stored G -> 0+ re-runs, the second certificate of the scaling step (records, three columns re-run, power).
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
check('the v1.0 search for the local Z3\' (packets with small raw coefficients) gives the ratio %.3f quoted in Sec. IX A as the source of the v1.0 error' % lo['worst']['rho_max'], ok and lo['worst']['rho_max'] <= 2.13, '%.0fs' % dt)
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

# ======================================================================== v1.1: the covering (Sec. 6, App. E)
import subprocess, cq_core as Q, cq_cover as CV, cq_checks as CQC, cl_link as KL
KL.set_tau_min(20.0)
KL.FAC.update({'low_S0': 1.0, 'low_S2': 1.0, 'up_S2': 1.0})     # the window check above (cf_windows) leaves its Lemma S factors in cl_link.FAC; certificate (ii) uses 1
cur = open('src/cq_phi.py').read()
r = subprocess.run([sys.executable, 'src/cq_make.py'], capture_output=True, text=True)
check('src/cq_phi.py is exactly what src/cq_make.py generates from the v1.0 term table (asserted replacements)', r.returncode == 0 and open('src/cq_phi.py').read() == cur)
cols = [(r_['ell_a'], r_['ell_b']) for r_ in L]
T2 = jl('cq_tails_ii.jsonl'); T1 = jl('cq_tails_i_w0.jsonl'); WG = jl('cq_wedge.jsonl')
check('(ii) mu tails: %d of %d columns closed (%d boxes)' % (sum(t['closed'] for t in T2), len(cols), sum(t['cells_evaluated'] for t in T2)),
      sorted((t['ell_a'], t['ell_b']) for t in T2 if t['closed']) == cols)
for kind, rows in (('cc', CC), ('cd', D)):
    t = [x for x in T1 if x['cert'] == kind]
    check('(i) %s mu tails and b -> 0 strip: %d of %d columns closed (%d boxes)' % ('D = 0' if kind == 'cc' else 'D = 10', sum(x['closed'] for x in t), len(rows), sum(x['cells_evaluated'] for x in t)),
          sorted((x['ell_a'], x['ell_b']) for x in t if x['closed']) == [(r_['ell_a'], r_['ell_b']) for r_ in rows])
check('(ii) wedge, boxes in (mu, r, beta, ell): %d of %d columns closed (%d boxes, %d skipped as inside the v1.0 boxes)' % (sum(w['closed'] for w in WG), len(cols), sum(w['cells_evaluated'] for w in WG), sum(w['skipped_boxes'] for w in WG)),
      sorted((w['ell_a'], w['ell_b']) for w in WG if w['closed']) == cols)
t0 = time.time(); same = True
samples = [('tails ii', Q.certify_tails_ii, (T2[0]['ell_a'], T2[0]['ell_b']), T2[0]), ('tails ii', Q.certify_tails_ii, (T2[-1]['ell_a'], T2[-1]['ell_b']), T2[-1]),
           ('wedge', Q.certify_wedge, (WG[0]['ell_a'], WG[0]['ell_b']), WG[0])]
mid = [w for w in WG if abs(w['ell_a'] - 455.26) < 1e-6][0]; samples.append(('wedge', Q.certify_wedge, (mid['ell_a'], mid['ell_b']), mid))
for kind in ('cc', 'cd'):
    t = sorted([x for x in T1 if x['cert'] == kind], key=lambda x: x['ell_a'])
    for x in (t[0], t[-1]): samples.append(('tails i ' + kind, lambda a, b, x=x, kind=kind: Q.certify_tails_i(a, b, x['delta'], kind), (x['ell_a'], x['ell_b']), x))
if FULL:
    samples = [('tails ii', Q.certify_tails_ii, (t['ell_a'], t['ell_b']), t) for t in T2] + [('wedge', Q.certify_wedge, (w['ell_a'], w['ell_b']), w) for w in WG] + \
              [('tails i', lambda a, b, x=x: Q.certify_tails_i(a, b, x['delta'], x['cert']), (x['ell_a'], x['ell_b']), x) for x in T1]
for name, f, args, ref in samples:
    d_ = f(*args); same &= d_['closed'] and d_['cells_evaluated'] == ref['cells_evaluated']
check('mu tails and wedge recomputed (%s): same cell counts' % ('all columns' if FULL else '%d sample columns' % len(samples)), same, '%.0fs' % (time.time() - t0))
cv = CV.run(100000, seed=1)
check('covering (ii): v1.0 boxes + mu tails + wedge + Lemma BE, %d random states, %d missed' % (cv['ii repaired (body + tails + wedge + BE)']['states'], cv['ii repaired (body + tails + wedge + BE)']['missed']), cv['ii repaired (body + tails + wedge + BE)']['missed'] == 0)
for kind in ('cc', 'cd'):
    k = 'i %s repaired (body + mu tails + b strip)' % kind
    check('covering (i) %s: %d random states, %d missed' % ('D = 0' if kind == 'cc' else 'D = 10', cv[k]['states'], cv[k]['missed']), cv[k]['missed'] == 0)
pw = {k: v['missed'] for k, v in cv.items() if k.startswith('POWER')}
check('the covering check finds every deliberately damaged table (%d tables; v1.0 table for (ii): %d missed)' % (len(pw), pw['POWER ii: v1.0 table alone']), all(v > 0 for v in pw.values()) and len(pw) == 11)
cq = CQC.run_all(0.1)
check('checks of the new steps at reduced scale: identity with the v1.0 functions, components, box soundness, tau factors; wrong versions broken %s' % cq['wrong_broken'],
      cq['F1_identity']['phi_b_maxdiff'] == 0 and cq['F1_identity']['ok_i_dis'] == 0 and cq['K2']['violations'] == 0 and cq['F3_ii']['violations'] == 0 and cq['F3_i']['violations'] == 0
      and cq['F3_w']['violations'] == 0 and sum(cq['K6']['increases'].values()) == 0 and all(cq['wrong_broken'].values()))
CPc = json.load(open('results/cp_checks.json'))
check('beta-box components (J2)-(J5) and the four wrong versions of stage CP (stored): violations %s, broken %s' % (CPc['J']['violations'], {k: v['broken'] for k, v in CPc['wrong'].items()}),
      sum(CPc['J']['violations'].values()) == 0 and all(v['broken'] for v in CPc['wrong'].values()))
cp0 = jl('cp_link_cells.jsonl')[0]
check('control: with beta boxes but ell fixed per column, [55.26, 65.26] does not close (%s); width-1 columns do' % cp0['reason'], not cp0['closed'] and all(c['closed'] for c in jl('cp_link_control_w1.jsonl')))
# ======================================================================== v1.1: the narrow focus (Sec. IX A), stored results
CM = json.load(open('results/cm_counter.json'))
check("local Z3' with one constant is false: P_theta/Pi*(mu_theta) = %.1f at x = 0.003 and %.1f at x = 1e-6" % (CM['headline'][0]['rho'], CM['best']['rho']), CM['headline'][0]['rho'] > 20 and CM['best']['rho'] > 70)
EO = json.load(open('results/co_envmin_opt.json'))
check('local envelope with the proved Psi_2* above T_0 at all six ell (numerical): ' + ', '.join('%s: %.4f > %.4f' % (k, v['min'], v['T0']) for k, v in EO.items()), all(v['min'] > v['T0'] for v in EO.values()) and len(EO) == 6)
RT = json.load(open('results/co_rtheta.json'))
check('Lemma R_theta below direct quadrature on %d states' % len(RT['states']), RT['all_ok'] and len(RT['states']) == 46)
CN = json.load(open('results/cn_search.json'))
check('joint-bound search: candidates C1, C2 survive, the four wrong versions are broken (%s)' % CN['verdict'], CN['verdict'] == {'C1': 'survived', 'C2': 'survived', 'W1': 'BROKEN', 'W2': 'BROKEN', 'W3': 'BROKEN', 'W4': 'BROKEN'})
# ---- v1.2: the end point G -> 0+ and the second certificate of the scaling step (verification/run_v12.py)
t0 = time.time(); r = subprocess.run([sys.executable, '-W', 'ignore', 'verification/run_v12.py'] + (['--recompute'] if FULL else []), capture_output=True, text=True, timeout=36000 if not FULL else 10 ** 6)
print(r.stdout, end='', flush=True)
check('v1.2: the G -> 0+ re-runs and the second certificate of the scaling step (verification/run_v12.py)', r.returncode == 0, '%.0fs' % (time.time() - t0))
print('\n%d/%d checks passed' % (sum(res), len(res))); sys.exit(0 if all(res) else 1)
