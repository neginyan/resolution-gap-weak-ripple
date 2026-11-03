"""Stage CF 2.3 (control of the search power): the counterexample search of ce_split_map.py applied to DELIBERATELY WRONG versions of Lemma S / S'
(one term removed).  A search with power must report violations for the wrong versions, and none for the true ones.  NUMERICAL.
-> results/cf_split_control.json.  Python 3.11."""
import json, math
import numpy as np
import ce_split_map as C

def factors_variant(eta_lo, eta_hi, f_lo, kappa, gap, drop):
    eta = eta_lo + eta_hi; f_hi = 1 - f_lo
    t = {'T6': eta_lo * f_lo, 'T3': 2 * eta_hi * math.sqrt(max(f_lo * f_hi, 0)), 'T1': eta_hi * f_hi, 'T2': kappa * f_hi,
         'T5': (0 if gap else eta_hi * f_lo), 'T4': (0 if gap else 2 * eta_lo * math.sqrt(max(f_lo * f_hi, 0)))}
    E = sum(v for k, v in t.items() if k != drop)
    B = (eta_lo * math.sqrt(f_lo) + kappa * math.sqrt(f_hi) + eta_hi) ** 2
    up = 1 + E + (0 if drop == 'B' else B) / (1 - eta)
    lo = 1 - E
    if drop != 'oldsandwich':
        lo = max(1 / (1 + eta), lo); up = min(1 / (1 - eta), up)
    return lo, up

def run(drop, gap_windows, form_gap, ntrial=150, seed=23, climb=30):
    rng = np.random.default_rng(seed); worst = 9.0; n = 0
    orig = C.factors_general
    C.factors_general = lambda a, b, c, d, g: factors_variant(a, b, c, d, g, drop)
    try:
        for it in range(ntrial):
            M_ = int(rng.integers(1, 300)); m = int(rng.integers(2 * M_ + 1, 2 * M_ + 600)) if gap_windows else int(rng.integers(M_ + 1, 2 * M_ + 1))
            Mh = min(m + int(rng.integers(1, 300)), 4000); start = int(rng.integers(1, max(2, M_)))
            nlo = np.arange(start, M_ + 1); nhi = np.arange(m, Mh + 1); coh = rng.random() < 0.5
            w = [nlo, np.ones(len(nlo)) if coh else rng.uniform(-1, 1, len(nlo)), np.zeros(len(nlo)) if coh else rng.uniform(0, 2 * np.pi, len(nlo)),
                 nhi, rng.choice([-1.0, 1.0], len(nhi)) * (1 if rng.random() < 0.5 else rng.random(len(nhi))), rng.uniform(0, 2 * np.pi, len(nhi)), rng.uniform(0.05, 0.95), rng.uniform(0, 1)]
            r = C.window_eval(*w)
            if r is None: continue
            best = C.margins(r, form_gap); ww = w
            for _ in range(climb):
                cand = [x.copy() if isinstance(x, np.ndarray) else x for x in ww]
                cand[1] = cand[1] * np.exp(rng.normal(0, 0.3, len(cand[1]))); cand[7] = float(np.clip(cand[7] + rng.normal(0, 0.1), 0, 1)); cand[6] = float(np.clip(cand[6] + rng.normal(0, 0.05), 0.02, 0.95))
                rc = C.window_eval(*cand)
                if rc is None: continue
                mc = C.margins(rc, form_gap)
                if mc < best: best, ww = mc, cand
            worst = min(worst, best); n += 1
            if worst < -1e-9 and it > 40: break
    finally:
        C.factors_general = orig
    return {'drop': drop, 'gap_windows': gap_windows, 'windows': n, 'min_margin': worst, 'violation_found': worst < -1e-9}

if __name__ == '__main__':
    import os
    FNL = 'results/cf_split_control.jsonl'
    done = set()
    if os.path.exists(FNL):
        for l in open(FNL):
            if l.strip(): d = json.loads(l); done.add((d['drop'], d['gap_windows']))
    out = [json.loads(l) for l in open(FNL) if l.strip()] if os.path.exists(FNL) else []
    for drop in ('none', 'T6', 'T3', 'T2', 'T1', 'B', 'oldsandwich', 'T4'):
        for gw in ((True, False) if drop not in ('T4',) else (False,)):
            if (drop, gw) in done: continue
            d = run(drop, gw, gw); d = {k: (v.item() if hasattr(v, 'item') else v) for k, v in d.items()}; out.append(d); print(d, flush=True)
            open(FNL, 'a').write(json.dumps(d) + '\n')
    json.dump(out, open('results/cf_split_control.json', 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
