"""Stage CL: recompute certificate (ii) (and other column lists) with the rounding-covering rule of cl_link.py, column by column,
on exactly the columns of the stored certificate.  Resumable; one output file per worker.  Python 3.11.
CLI: python cl_run.py WORKER NWORKERS"""
import json, math, sys, time
import cl_link as K
w, nw = int(sys.argv[1]), int(sys.argv[2])
cols = sorted([json.loads(l) for l in open('results/ce_link_cells.jsonl') if l.strip()], key=lambda r: r['ell_a'])
cols = [c for c in cols if c.get('tau_min') == 20.0 and abs(c['lbmin_rel'] - math.log(1e4)) < 1e-9 and c['D'] == 10.0 and c['ell_b'] <= 960.81]
# interleave by cost: big columns alternate between workers
mine = [c for i, c in enumerate(cols) if i % nw == w]
FN = 'results/cl_link_cells_w%d.jsonl' % w
have = set()
try:
    for l in open(FN):
        if l.strip(): have.add(round(json.loads(l)['ell_a'], 6))
except FileNotFoundError: pass
K.set_tau_min(20.0)
for c in mine:
    if round(c['ell_a'], 6) in have: continue
    t0 = time.time(); d = K.certify(c['ell_a'], c['ell_b'], 10.0, c['lbmin_rel'])
    d.update({'ell_b': c['ell_b'], 'lbmin_rel': c['lbmin_rel'], 'tau_min': 20.0, 'seconds': round(time.time() - t0, 1), 'old_cells': c['cells_evaluated'], 'rule': 'CL rounding-covering, S = 10'})
    open(FN, 'a').write(json.dumps(d) + '\n')
    print('ell [%.2f, %.2f] %s cells %s (old %d) %.0fs' % (c['ell_a'], c['ell_b'], 'CLOSED' if d['closed'] else 'OPEN ' + d.get('reason', ''), d.get('cells_evaluated', d.get('cells')), c['cells_evaluated'], d['seconds']), flush=True)
    if not d['closed']: print('STOP: a column did not close', flush=True); break
