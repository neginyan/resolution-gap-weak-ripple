"""Stage CB section 2.1: the bridge to tau < 80.
(A) NUMERICAL (reduced formula, judged by hi): the hardest Rudin-Shapiro band window when the state is FORCED into tau = x/b <= tau_c (b >= x/tau_c),
    by the same continuous optimisation as bz_asym.py but with b in [x/tau_c, min(1, u_M^2/80)] (the band harmonics stay separated: n^2 tau >= u_m^2/b
    >= 80).  The reduced formula is Lemma R without the theta_tau factors; for the band part these are e^{-m^2 tau/2}-small (Lemma H, cb_split.py), so
    this is the value of the band part in the narrow-envelope regime.  -> results/cb_tau_scan.jsonl (one line per (sigma, tau_c), resumable).
(B) COMPUTER-ASSISTED PROOF (1-D cells in tau): for tau in [tau_c, 80] the tau-dependent factors of Lemma R / Z1R,
        F_tau = sqrt((1 - eps_1)(1 - nu_0)/(1 + z)) / ((1 + nu_0)/(1 - z)) * (1 - nu_0)/(1 + z) ... (the product of every factor in T_B that depends on tau),
    are bounded below on cells [tau_a, tau_b] by their values at tau_a with beta at its largest value (every nu_k(tau, beta) = 2 sum n^k e^{-n^2 tau/2 + n beta}
    decreases in tau and increases in beta; beta = alpha k' V <= sigma tau/sqrt(2x) <= tau/sqrt(2 T_s) with x > T_s sigma^2, T_s = 4 (v1.5 split point);
    here we use beta <= tau * 1e-6, valid for x >= 5e11 sigma^2, and report the cell-wise minimum.  In-cell audit: 50 interior points per cell.
    -> results/cb_tau_cells.json.  Python 3.11."""
import json, math, os, sys
import numpy as np
from scipy.optimize import minimize
import by_eval as Y

FN = 'results/cb_tau_scan.jsonl'


def minimise_b(lsig, x, S0, S2, eta, lbmin, lbmax, side='hi'):
    lL2max = -lsig
    def unpack(z):
        lL2 = min(z[0], lL2max); ls2 = min(max(z[1], 2 * lsig), lL2 - math.log(4)); lb = min(max(z[2], lbmin), lbmax)
        return lL2, ls2, lb
    f = lambda z: sum(Y.branches_log(lsig, x, S0, S2, eta, *unpack(z), side))
    best = None
    for dL in (0.5, 0.0, -2.0, -6.0):
        for ds in (0.5, -1.0):
            for lb0 in (lbmin, lbmin + 2.0, 0.5 * (lbmin + lbmax), lbmax):
                r = minimize(f, [0.5 * lsig + dL, lsig + ds, lb0], method='Nelder-Mead', options={'xatol': 1e-9, 'fatol': 1e-14, 'maxiter': 3000})
                if best is None or r.fun < best.fun: best = r
    return best.fun, unpack(best.x)


def hardest_forced(sig_exp, tau_c, lx0, u0):
    lsig = -sig_exp * math.log(10); ell = sig_exp * math.log(10)
    def inner(z):
        lx, u = z; x = 10.0 ** lx
        if u < 0.5 or lx > math.log10(0.003): return 1e9, None        # x >= K_0 = 0.003 is covered by v1.5 Theorem F11
        if math.log(x / tau_c) >= 0: return 1e9, None                  # b = sigma^2/V <= 1: no state with tau <= tau_c
        m = max(2, int(round(u / math.sqrt(x)))); p = Y.profile_rs(x, m)
        if p['eta'] >= 1: return 1e9, None
        uM = (m + p['W']) * math.sqrt(x)
        lbmin = math.log(x / tau_c); lbmax = math.log(min(1.0, uM * uM / 80))        # tau <= tau_c, band separated
        if lbmax <= lbmin: return 1e9, None
        g, st = minimise_b(lsig, x, p['S0'], p['S2'], p['eta'], lbmin, lbmax)
        return g * math.sqrt(ell), (st, p)
    best = None
    for dz in ((0, 0), (sig_exp * 0.1, -0.6), (-sig_exp * 0.1, 0.6), (sig_exp * 0.25, 0.0)):
        z0 = [lx0 + dz[0], u0 + dz[1]]
        r = minimize(lambda z: inner(z)[0], z0, method='Nelder-Mead', options={'xatol': 1e-3, 'fatol': 1e-6, 'maxiter': 200,
                     'initial_simplex': [z0, [z0[0] + max(1, 0.05 * sig_exp), z0[1]], [z0[0], z0[1] + 0.5]]})
        if best is None or r.fun < best.fun: best = r
    g, (st, p) = inner(best.x)
    lx, u = best.x
    return {'sig_exp': sig_exp, 'tau_c': tau_c, 'judged_by': 'hi', 'gsl_hi_forced_tau_le_tau_c': g, 'log10_x': lx, 'u_m': u,
            'log10_tau': (lx * math.log(10) - st[2]) / math.log(10), 'b_over_sigma_log10': (st[2] - lsig) / math.log(10), 'mu': 10 ** lx * p['S2'] / p['S0']}


def tau_cells(tau_c_list=(8.0, 12.0, 20.0, 40.0), beta_rel=1e-6, ncell=400, Qb=0.25):
    out = {}
    def factors(tau, beta):
        n = np.arange(1, 80); en = np.exp(-n * n * tau / 2 + n * beta)
        nu0 = 2 * en.sum(); z = Qb * nu0 / 2; eps1 = nu0 * (1 + 2 * Qb) / (1 - z)
        k1 = 1 - eps1; k2 = (1 - nu0) / (1 + z); c0 = (1 + nu0) / (1 - z)
        return math.sqrt(k1 * k2) / c0, nu0
    for tc in tau_c_list:
        edges = np.geomspace(tc, 80.0, ncell + 1); worst = 9.0; audit_bad = 0
        for a, b in zip(edges[:-1], edges[1:]):
            Fcell, nu0 = factors(a, beta_rel * b)                      # every factor is monotone: worst at the left end and the largest beta
            worst = min(worst, Fcell)
            for t in np.linspace(a, b, 50):                             # in-cell audit
                if factors(t, beta_rel * t)[0] < Fcell * (1 - 1e-12): audit_bad += 1
        out['tau_c=%g' % tc] = {'min_tau_factor_F': worst, 'one_minus_F': 1 - worst, 'nu0_at_tau_c': factors(tc, beta_rel * tc)[1], 'cells': ncell,
                                'audit_points': ncell * 50, 'audit_violations': audit_bad}
        print('tau in [%g, 80]: tau-factor of Z1R >= %.8f  (1 - F = %.2e), nu0(tau_c) = %.2e, audit violations %d' % (tc, worst, 1 - worst, out['tau_c=%g' % tc]['nu0_at_tau_c'], audit_bad))
    json.dump(out, open('results/cb_tau_cells.json', 'w'), indent=1)


if __name__ == '__main__':
    if sys.argv[1:2] == ['cells']:
        tau_cells(); sys.exit()
    have = set()
    if os.path.exists(FN):
        for l in open(FN):
            if l.strip(): d = json.loads(l); have.add((d['sig_exp'], d['tau_c']))
    ref = {json.loads(l)['sig_exp']: json.loads(l) for l in open('results/bz_asym.jsonl') if l.strip()}
    for e in [int(a) for a in sys.argv[1:]] or [24, 48, 96, 150, 300]:
        for tc in (80.0, 1.0, 1e-3):
            if (e, tc) in have: continue
            d = hardest_forced(e, tc, ref[e]['log10_x'], ref[e]['u_m'])
            d['separated_reference_hi'] = ref[e]['gsl_hi']
            open(FN, 'a').write(json.dumps(d) + '\n')
            print('sigma 1e-%d tau <= %g: hi g√ℓ %.4f (separated %.4f)  log10 x %.2f u %.2f log10 tau %.2f log10 b/sigma %.1f mu %.1f' % (
                e, tc, d['gsl_hi_forced_tau_le_tau_c'], d['separated_reference_hi'], d['log10_x'], d['u_m'], d['log10_tau'], d['b_over_sigma_log10'], d['mu']), flush=True)
