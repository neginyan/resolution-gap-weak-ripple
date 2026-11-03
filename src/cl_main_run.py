"""Stage CL: certificate (i) (main reduced inequality; cc_cells D = 0, 301 columns, and cd_cells D = 10, 686 columns) recomputed with the
rounding-covering rule of cl_cells.py on exactly the stored columns.  Resumable; one file per worker.  CLI: python cl_main_run.py WORKER NWORKERS"""
import json, math, sys, time
import cc_cells as CC, cd_cells as CD, cl_cells as L
CC.cell_ok = L.cell_ok
w, nw = int(sys.argv[1]), int(sys.argv[2])
cc = sorted([r for r in (json.loads(l) for l in open('results/cc_cells.jsonl') if l.strip()) if r['ell_a'] < 3062.59], key=lambda r: r['ell_a'])
cd = sorted([json.loads(l) for l in open('results/cd_cells.jsonl') if l.strip()], key=lambda r: r['ell_a'])
jobs = [('cc', r) for r in cc] + [('cd', r) for r in cd]
mine = [j for i, j in enumerate(jobs) if i % nw == w]
FN = 'results/cl_main_cells_w%d.jsonl' % w
have = set()
try:
    for l in open(FN):
        if l.strip(): d = json.loads(l); have.add((d['cert'], round(d['ell_a'], 6)))
except FileNotFoundError: pass
for kind, r in mine:
    if (kind, round(r['ell_a'], 6)) in have: continue
    L.MIN_SLACK[0] = math.inf; t0 = time.time()
    d = (CC.certify if kind == 'cc' else CD.certify)(r['ell_a'], r['ell_b'], r['delta'])
    rec = {'cert': kind, 'ell_a': r['ell_a'], 'ell_b': r['ell_b'], 'delta': r['delta'], 'closed': d['closed'], 'cells_evaluated': d.get('cells_evaluated', d.get('cells')),
           'old_cells': r['cells_evaluated'], 'audit_violations': d.get('audit_violations'), 'audit_points': d.get('audit_points'), 'min_slack_over_S_delta': L.MIN_SLACK[0], 'seconds': round(time.time() - t0, 1)}
    open(FN, 'a').write(json.dumps(rec) + '\n')
    print(rec, flush=True)
    if not d['closed']: print('STOP: a column did not close', flush=True); break
