"""Stage CJ: windows WITHOUT (P2) (no dominant fundamental) against Lemma Z3' -- period averages only.  Python 3.11.
Beats (two neighbouring harmonics) and combs of narrow spikes, raw f >= 0, blurred with K = x, scaled to sum q0 <= Qbar.
Z3' is a lemma under (P2); this records by how much windows outside (P2) exceed Pi*(mu).  [numerical]  -> results/cj_nop2.json"""
import json, math
import numpy as np
import ci_local as L
import cc_core as C

rows = []
T = L.T
fams = [('beats m=%d' % m, 1 + 0.5 * np.cos(m * T) + 0.5 * np.cos((m + 1) * T - 1.3), m) for m in (60, 200)]
fams += [('comb n=%d' % n, np.exp(-(1 - np.cos(n * T)) / (n * n * d * d)), n) for n, d in ((40, 0.02), (120, 0.005))]
for name, f, n in fams:
    for mu_t in (2.0, 6.0, 20.0):
        x = mu_t / n ** 2
        R, R1, S0, S2, qs, eps = L.ripple_from_raw(f.copy(), x)
        mu = x * S2 / S0; P = x * S2
        rows.append({'window': name, 'mu': mu, 'P_over_Pistar': P / math.exp(C.log_Pistar(mu))})
out = {'rows': rows, 'min_ratio': min(r['P_over_Pistar'] for r in rows), 'max_ratio': max(r['P_over_Pistar'] for r in rows)}
json.dump(out, open('results/cj_nop2.json', 'w'), indent=1)
print(out['min_ratio'], out['max_ratio'])
for r in rows: print(r)
