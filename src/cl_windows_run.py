"""Stage CL: certificate (iii) (17 thin-block windows) recomputed with the rounding-covering rule (cf_windows.certify_window evaluated with
cl_link.phi_cells instead of ce_link.phi_cells).  Resumable.  Python 3.11.  -> results/cl_windows.jsonl"""
import json, math, os, time
import cl_link, cf_windows, ce_split_map as SM
cf_windows.K = cl_link
FN = 'results/cl_windows.jsonl'
ref = {48: (-17.49, 3.92), 96: (-40.88, 6.39), 300: (-142.12, 12.46)}
W = [json.loads(l) for l in open('results/cf_windows.jsonl') if l.strip()]
have = set()
if os.path.exists(FN):
    have = {json.loads(l)['index'] for l in open(FN) if l.strip()}
for i, w in enumerate(W):
    if i in have: continue
    e = w['sig_exp']; x = 10 ** ref[e][0]; u = ref[e][1]
    ws = SM.window_general(x, u, w['eta_target'], w['u_block_over_u_band'] * u, w['block_width_decay_lengths'])
    S0 = ws['S0_lo'] + ws['S0_hi']; S2 = ws['S2_lo'] + ws['S2_hi']
    t0 = time.time(); d = cf_windows.certify_window(e * math.log(10), math.log(x * S2 / S0), math.log(x * S2), math.log(x), w['factors'])
    rec = {'index': i, 'sig_exp': e, 'closed': d['closed'], 'cells_evaluated': d.get('cells_evaluated', d.get('cells')), 'old_cells': w['cells_evaluated'], 'seconds': round(time.time() - t0, 1)}
    open(FN, 'a').write(json.dumps(rec) + '\n'); print(rec, flush=True)
    if not d['closed']: print('STOP: a window did not close', flush=True); break
