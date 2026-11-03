"""Stage CC section 1.2 (NUMERICAL, judged by hi): the 'intermediate' windows, eta between 1/ell and Qbar.
A Rudin-Shapiro band (carries mu, eta ~ 0) plus a coherent flat block (all phases 0) of N_c harmonics just below the band start, whose blurred
peak sets eta.  Non-negativity (proved bound): f >= 1 - c_RS (1 + sqrt(2W)) - A_c (1 + N_c) >= 0 with the block's raw amplitude A_c (a flat cosine sum
of N_c terms is >= -N_c, crude); (P2): the fundamental equals the largest raw amplitude; (P3): the total blurred sum is scaled to Qbar = 1/4.
eta = (RS bound of the band) + (block's blurred coefficient sum) (an upper bound of sup|R|, used on the 'hi' side as in BY).
For sigma = 1e-24 .. 1e-300 and target eta in {1/ell, 0.003, 0.01, 0.03, 0.05, 0.1, 0.2}: minimise hi over the state and over (x, u_m) near the band
optimum.  -> results/cc_mid.jsonl (resumable).  Python 3.11."""
import json, math, os
from scipy.optimize import minimize
import by_eval as Y

FN = 'results/cc_mid.jsonl'; QB = 0.25


def profile(x, u, eta_t, Nc_dl=0.2, uc=None):
    m = max(3, int(round(u / math.sqrt(x)))); p = Y.profile_rs(x, m)          # band, Qbar-limited, eta_band = p['eta']
    K = x
    if uc is None: Nc = max(1, int(round(Nc_dl / (m * K)))); mc = max(2, m - Nc)
    else:
        mc = max(2, int(round(uc / math.sqrt(x)))); Nc = max(1, int(round(Nc_dl / (mc * K)))); Nc = min(Nc, m - mc)
        if Nc < 1: return None
    # block blurred coefficients q_n = A e^{-n^2 K/2}, n = mc .. m-1; sum_block = A * s_b
    s_b = Y.gsum(0, mc, K / 2, mc + Nc - 1); S0b = Y.gsum(0, mc, K, mc + Nc - 1); S2b = Y.gsum(2, mc, K, mc + Nc - 1)
    need = max(0.0, eta_t - p['eta'])
    A = need / s_b if s_b > 0 else 0.0
    # non-negativity of the raw profile with the band's raw amplitude c: 1 - c(1 + sqrt(2W)) - A (1 + Nc) >= 0
    if p['c'] * (1 + math.sqrt(2 * p['W'])) + A * (1 + Nc) > 1: return None
    a1 = max(p['c'], A)                                                         # (P2)
    e1 = math.exp(-K / 2)
    S0 = p['S0'] - (p['c'] * e1) ** 2 + (a1 * e1) ** 2 + A * A * S0b
    S2 = p['S2'] - (p['c'] * e1) ** 2 + (a1 * e1) ** 2 + A * A * S2b
    Q = p['Q'] - p['c'] * e1 + a1 * e1 + A * s_b
    eta = p['eta'] - p['c'] * e1 + a1 * e1 + A * s_b
    if Q > QB:
        s = QB / Q; S0 *= s * s; S2 *= s * s; eta *= s; Q = QB
    return {'S0': S0, 'S2': S2, 'eta': eta, 'mu': x * S2 / S0, 'A_block': A, 'Nc': Nc}


def hardest(e, eta_t, lx0, u0):
    lsig = -e * math.log(10); ell = e * math.log(10)
    def f(z):
        x = 10.0 ** z[0]
        if z[1] < 0.5 or z[0] > math.log10(0.003) or not (0.05 <= z[2] < z[1]) or not (-4 <= z[3] <= 1): return 1e9
        pr = profile(x, z[1], eta_t, Nc_dl=10 ** z[3], uc=z[2])
        if pr is None or pr['eta'] >= 0.9: return 1e9
        r = Y.minimise(lsig, x, pr['S0'], pr['S2'], pr['eta'], 'hi')
        return r['g'] * math.sqrt(ell) if r else 1e9
    best = None
    for dz in ((0, 0), (e * 0.03, -0.5), (-e * 0.03, 0.5)):
        for uc0 in (0.3, 1.0, 2.0):
            for ln0 in (-2.0, -0.7):
                z0 = [lx0 + dz[0], u0 + dz[1], min(uc0, u0 + dz[1] - 0.2), ln0]
                if f(z0) >= 1e9: continue
                r = minimize(f, z0, method='Nelder-Mead', options={'xatol': 1e-3, 'fatol': 1e-6, 'maxiter': 300})
                if best is None or r.fun < best.fun: best = r
    if best is None:
        return {'sig_exp': e, 'eta_target': eta_t, 'judged_by': 'hi', 'gsl_hi': None, 'feasible': False}
    pr = profile(10.0 ** best.x[0], best.x[1], eta_t, Nc_dl=10 ** best.x[3], uc=best.x[2])
    return {'sig_exp': e, 'eta_target': eta_t, 'judged_by': 'hi', 'gsl_hi': best.fun, 'log10_x': best.x[0], 'u_m': best.x[1],
            'u_block': best.x[2], 'block_width_decay_lengths': 10 ** best.x[3], 'feasible': True, 'eta': pr['eta'] if pr else None, 'mu': pr['mu'] if pr else None, 'H1_factor_(1-eta)/(1+eta)': (1 - pr['eta']) / (1 + pr['eta']) if pr else None}


if __name__ == '__main__':
    import sys
    have = set()
    if os.path.exists(FN):
        for l in open(FN):
            if l.strip(): d = json.loads(l); have.add((d['sig_exp'], d['eta_target']))
    ref = {json.loads(l)['sig_exp']: json.loads(l) for l in open('results/bz_asym.jsonl') if l.strip()}
    for e in [int(a) for a in sys.argv[1:]] or [24, 48, 96, 300]:
        ell = e * math.log(10)
        for et in (1 / ell, 0.003, 0.01, 0.03, 0.05, 0.1, 0.2):
            if (e, et) in have: continue
            d = hardest(e, et, ref[e]['log10_x'], ref[e]['u_m']); d['band_reference_hi'] = ref[e]['gsl_hi']
            open(FN, 'a').write(json.dumps(d) + '\n')
            if d['gsl_hi'] is None: print('sigma 1e-%d eta_target %.4f: no feasible profile' % (e, et), flush=True); continue
            print('sigma 1e-%d eta_target %.4f: hi g√ℓ %.4f (band %.4f)  eta %.4f mu %.2f u_block %.2f' % (e, et, d['gsl_hi'], d['band_reference_hi'], d['eta'] or -1, d['mu'] or -1, d['u_block']), flush=True)
