"""Stage CC: the counterexample to conjecture C (eta e^{mu/2} <= C), and how these windows fare (NUMERICAL, judged by hi as always).
Window: a Rudin-Shapiro band (the hardest band of BZ: carries mu) plus a WIDE coherent block of N_c ~ 1e8 low harmonics (all phases 0) with a
tiny raw amplitude A.  The block's blurred coefficients add up coherently (eta ~ N_c A e^{-u_c^2/2}, an l1 quantity) while its energy is
N_c A^2 (an l2 quantity), so mu is almost unchanged.  Proved properties of the window (cc_mid.profile): f >= 1 - c_RS(1 + sqrt(2W)) - A(1 + N_c) >= 0,
(P2) by a_1 = max amplitude, (P3) by scaling the blurred sum to Qbar = 1/4.  eta is the upper bound used on the hi side (band RS bound + block sum);
the true sup is at least (block sum) - (band RS bound), within 1e-3 relative of the value used.  -> results/cc_counter.json.  Python 3.11."""
import json, math
import cc_mid as M, by_eval as Y
rows = []
for e, lx, u in ((48, -17.49, 3.92), (96, -40.88, 6.39), (300, -142.12, 12.46)):
    ell = e * math.log(10); lsig = -e * math.log(10); x = 10 ** lx
    for et in (0.0, 0.03, 0.1, 0.2):
        pr = M.profile(x, u, et, Nc_dl=0.2, uc=0.3) if et > 0 else M.profile(x, u, 0.0)
        if pr is None: continue
        g = {s: Y.minimise(lsig, x, pr['S0'], pr['S2'], pr['eta'], s)['g'] * math.sqrt(ell) for s in ('hi', 'mid', 'lo')}
        mb = max(3, int(round(u / math.sqrt(x)))); pb = Y.profile_rs(x, mb)
        eta_lower = max(0.0, pr['eta'] - 2 * pb['eta'])
        rows.append({'sig_exp': e, 'log10_x': lx, 'u_m': u, 'eta_target': et, 'eta_upper': pr['eta'], 'eta_lower': eta_lower, 'mu': pr['mu'],
                     'log10_eta_e_mu2_lower': (math.log(eta_lower) + pr['mu'] / 2) / math.log(10) if eta_lower > 0 else None,
                     'A_block': pr['A_block'], 'N_block': pr['Nc'], 'gsl_hi': g['hi'], 'gsl_mid': g['mid'], 'gsl_lo': g['lo']})
        print(rows[-1], flush=True)
json.dump(rows, open('results/cc_counter.json', 'w'), indent=1)
