"""Stage CE 2.3: the windows of the map NOT covered by Corollary S (Lambda_S T_10 kappa_tau < 1): judged by hi (upper side, as always), and the split lower
bound minimised over the state (NUMERICAL, the same reduced formula as CD).  -> results/ce_split_uncovered.json.  Python 3.11."""
import json, math
import by_eval as Y, cd_split as S, ce_split_map as C
rows = []
for r in json.load(open('results/ce_split_map.json'))['map']:
    if r['covered_by_corollary_S']: continue
    e = r['sig_exp']; ref = {48: (-17.49, 3.92), 96: (-40.88, 6.39), 300: (-142.12, 12.46)}[e]; x = 10 ** ref[0]; u = ref[1]
    ws = C.window_general(x, u, r['eta_target'], r['u_block_over_u_band'] * u, r['block_width_decay_lengths'])
    lam, (l0, u0, l2, u2, f0, f2) = C.lam_S(ws, 10.0 ** (-e))
    S0 = ws['S0_lo'] + ws['S0_hi']; S2 = ws['S2_lo'] + ws['S2_hi']; eta = ws['eta_lo'] + ws['eta_hi']; ell = e * math.log(10); lsig = -ell
    hi = Y.minimise(lsig, x, S0, S2, eta, 'hi')['g'] * math.sqrt(ell)
    mid = Y.minimise(lsig, x, S0, S2, eta, 'mid')['g'] * math.sqrt(ell)
    split = S.minimise_split(lsig, x, S0, S2, {'low_S0': l0, 'up_S0': u0, 'low_S2': l2, 'up_S2': u2}) * math.sqrt(ell)
    rows.append({**r, 'gsl_hi': hi, 'gsl_mid': mid, 'gsl_split': split})
    print('1e-%d eta %.3f u_c/u %.2f width %g: hi %.4f mid %.4f split %.4f  (Corollary S product %.4f)' % (e, eta, r['u_block_over_u_band'], r['block_width_decay_lengths'], hi, mid, split, r['product_with_T10_kappa']), flush=True)
json.dump(rows, open('results/ce_split_uncovered.json', 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
