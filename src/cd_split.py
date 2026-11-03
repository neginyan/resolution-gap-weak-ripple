"""Stage CD 1.1: LEMMA S (the split sandwich) and its values on the CC counterexample windows.  Python 3.11.

LEMMA S [proof].  Let R(t) = R_lo(t) + R_hi(t) be a real 2pi-periodic function with no constant term, R_lo having frequencies in [1, M] and R_hi
frequencies >= m, with the GAP 2M < m.  Let psi(w) in [0, 1], h = 1 + r, r = psi(w) R(t), eta = sup|R| < 1, and let Y(t, w) = phi(w) y(t) with
y = R or y = R' (so y = y_lo + y_hi with the same frequency split).  Write <.> for the t-period average integrated against any positive measure
rho(dw) (product structure; this is the form of J_per and W_per in (R1) and of the main F-term of (R4) in Lemma R), ||.|| for the L2 norm in it.
Put f_lo = ||Y_lo||^2/||Y||^2, f_hi = 1 - f_lo (Y_lo, Y_hi are orthogonal: disjoint frequencies), eta_lo = sup|R_lo|, eta_hi = sup|R_hi|,
kappa = min(eta_lo, sup|y_hi| ||R_lo||_t / ||y_hi||_t)   (||.||_t: L2 norm of the period average; ||R_lo||_t^2 = S0_lo/2), and
    E = eta_lo f_lo + 2 eta_hi sqrt(f_lo f_hi) + (eta_hi + kappa) f_hi,     B = (eta_lo sqrt(f_lo) + kappa sqrt(f_hi) + eta_hi)^2.
Then
    ||Y||^2 max(1/(1 + eta), 1 - E)  <=  <Y^2/h>  <=  ||Y||^2 min(1/(1 - eta), 1 + E + B/(1 - eta)).
Proof.  1/(1 + r) = 1 - r + r^2/(1 + r), and 0 <= r^2/(1 + r) <= r^2/(1 - eta).  So <Y^2/h> lies between ||Y||^2 - <Y^2 r> and
||Y||^2 - <Y^2 r> + <Y^2 r^2>/(1 - eta); the old sandwich gives the other halves of max/min.  Expand Y^2 r = (Y_lo^2 + 2 Y_lo Y_hi + Y_hi^2)(r_lo + r_hi):
  <Y_lo^2 r_hi> = int phi^2 psi drho * <y_lo^2 R_hi>_t = 0, because y_lo^2 has frequencies <= 2M < m (gap);
  <Y_lo Y_hi r_lo> = int phi^2 psi drho * <(y_lo R_lo) y_hi>_t = 0, the same;
  |<Y_lo^2 r_lo>| <= eta_lo ||Y_lo||^2;  |2 <Y_lo Y_hi r_hi>| <= 2 eta_hi ||Y_lo|| ||Y_hi|| (Cauchy-Schwarz);  |<Y_hi^2 r_hi>| <= eta_hi ||Y_hi||^2;
  |<Y_hi^2 r_lo>| = int phi^2 psi drho |<y_hi^2 R_lo>_t| <= ||Y_hi||^2 min(eta_lo, sup|y_hi| ||R_lo||_t/||y_hi||_t)   (0 <= psi <= 1; Cauchy-Schwarz in t:
     <y_hi^2 |R_lo|>_t <= ||y_hi^2||_t ||R_lo||_t <= sup|y_hi| ||y_hi||_t ||R_lo||_t).
  Hence |<Y^2 r>| <= E ||Y||^2.  For r^2: ||Y r|| <= ||Y r_lo|| + ||Y r_hi|| (Minkowski), ||Y r_hi|| <= eta_hi ||Y||, and
  ||Y r_lo|| <= ||Y_lo r_lo|| + ||Y_hi r_lo|| <= eta_lo ||Y_lo|| + kappa ||Y_hi|| (psi^2 <= 1, <y_hi^2 R_lo^2>_t <= sup|y_hi|^2 ||R_lo||_t^2).  So <Y^2 r^2> <= B ||Y||^2.  QED.
Meaning: the low block enters the S2 channel (y = R', f_lo = S2_lo/S2 tiny) only through eta_hi + kappa, kappa ~ crest(R'_hi) sqrt(S0_lo/2), its ENERGY,
not its sup eta_lo.  In the S0 channel it enters with the weight f_lo = S0_lo/S0.

Sup bounds used for the Rudin-Shapiro band (m^2 K >= 1, so d_k = (m+k) e^{-(m+k)^2 K/2} decreases) [proof, Abel summation as in by_eval.py]:
  sup|R'_hi| <= c C_RS e^{-m^2 K/2} [m + m sqrt(pi/(mK))/2 + (sqrt(pi)/4)(mK)^{-3/2}]   (d_k <= (m+k) e^{-m^2K/2} e^{-mkK}; d(k)/(2 sqrt k) decreasing).
Parts: (1) numerical check of Lemma S on random gapped windows (both sides: low/high frequency, many/few harmonics, coherent/random phases);
(2) the CC counterexample windows: old lo / mid / hi and the new split lower bound 'split'.  -> results/cd_split.json.
"""
import json, math, sys
import numpy as np
import by_eval as Y, cc_mid as M

C_RS = 2 + math.sqrt(2)


def lemma_S_factors(eta_lo, eta_hi, f_lo, kappa):
    eta = eta_lo + eta_hi
    f_hi = 1 - f_lo
    E = eta_lo * f_lo + 2 * eta_hi * math.sqrt(max(f_lo * f_hi, 0)) + (eta_hi + kappa) * f_hi
    B = (eta_lo * math.sqrt(f_lo) + kappa * math.sqrt(f_hi) + eta_hi) ** 2
    return max(1 / (1 + eta), 1 - E), min(1 / (1 - eta), 1 + E + B / (1 - eta))


# ---------------- (1) numerical check on random gapped windows ----------------
def check_random(ntrial=400, seed=1):
    rng = np.random.default_rng(seed); worst = []
    N = 4096; t = np.arange(N) * 2 * np.pi / N
    for it in range(ntrial):
        M_ = int(rng.integers(1, 300)); m = int(rng.integers(2 * M_ + 1, 2 * M_ + 400)); Mh = m + int(rng.integers(1, 300))
        if Mh >= N // 2 - 2: Mh = N // 2 - 3
        nlo = np.arange(1, M_ + 1); nhi = np.arange(m, Mh + 1)
        kind = it % 4
        alo = rng.random(len(nlo)) if kind in (0, 1) else np.ones(len(nlo))                 # coherent low blocks (kind 2, 3)
        ahi = rng.choice([-1.0, 1.0], len(nhi)) * (rng.random(len(nhi)) if kind in (0, 2) else 1.0)
        plo = rng.uniform(0, 2 * np.pi, len(nlo)) if kind in (0, 1) else np.zeros(len(nlo))
        phi_ = rng.uniform(0, 2 * np.pi, len(nhi))
        Rlo = (alo[:, None] * np.cos(nlo[:, None] * t + plo[:, None])).sum(0); Rhi = (ahi[:, None] * np.cos(nhi[:, None] * t + phi_[:, None])).sum(0)
        dRlo = (-alo[:, None] * nlo[:, None] * np.sin(nlo[:, None] * t + plo[:, None])).sum(0); dRhi = (-ahi[:, None] * nhi[:, None] * np.sin(nhi[:, None] * t + phi_[:, None])).sum(0)
        # target total eta in (0, 0.9), random share of the low block
        target = rng.uniform(0.02, 0.9); share = rng.uniform(0, 1)
        sl = share * target / max(np.max(np.abs(Rlo)), 1e-300); Rlo *= sl; dRlo *= sl
        sh = (1 - share) * target / max(np.max(np.abs(Rhi)), 1e-300); Rhi *= sh; dRhi *= sh
        R = Rlo + Rhi
        eta_lo = np.max(np.abs(Rlo)); eta_hi = np.max(np.abs(Rhi))      # sup on a grid 8x finer than the top frequency: measured
        for psi in (1.0, 0.6, 0.1):
            for (ylo, yhi) in ((Rlo, Rhi), (dRlo, dRhi)):
                y = ylo + yhi; h = 1 + psi * R
                exact = np.mean(y * y / h); nY = np.mean(y * y); f_lo = np.mean(ylo * ylo) / nY
                kappa = min(eta_lo, np.max(np.abs(yhi)) * math.sqrt(np.mean(Rlo * Rlo)) / math.sqrt(np.mean(yhi * yhi)))
                lo, up = lemma_S_factors(eta_lo, eta_hi, f_lo, kappa)
                worst.append((exact / nY - lo, up - exact / nY))
    w = np.array(worst)
    return {'trials': ntrial, 'checks': len(w), 'min_margin_lower': float(w[:, 0].min()), 'min_margin_upper': float(w[:, 1].min()),
            'violations': int(np.sum(w < -1e-12)), 'note': 'sup measured on a 4096-point grid (top frequency <= 2047 rendered exactly by the grid average)'}


# ---------------- (2) counterexample windows ----------------
def window_split(x, u, eta_t, Nc_dl=0.2, uc=0.3):
    """the cc_mid.profile window, split into lo = fundamental + block, hi = Rudin-Shapiro band, with the sup bounds."""
    K = x; m = max(3, int(round(u / math.sqrt(x)))); p = Y.profile_rs(x, m); c = p['c']; W = p['W']
    mc = max(2, int(round(uc / math.sqrt(x)))); Nc = max(1, int(round(Nc_dl / (mc * K)))); Nc = min(Nc, m - mc)
    top = mc + Nc - 1
    assert 2 * top < m, 'no spectral gap'
    s_b = Y.gsum(0, mc, K / 2, top); S0b = Y.gsum(0, mc, K, top); S2b = Y.gsum(2, mc, K, top)
    A = max(0.0, eta_t - p['eta']) / s_b if eta_t > 0 else 0.0
    assert c * (1 + math.sqrt(2 * W)) + A * (1 + Nc) <= 1
    a1 = max(c, A); e1 = math.exp(-K / 2)
    S0_hi = c * c * Y.gsum(0, m, K, m + W - 1); S2_hi = c * c * Y.gsum(2, m, K, m + W - 1)
    S0_lo = (a1 * e1) ** 2 + A * A * S0b; S2_lo = (a1 * e1) ** 2 + A * A * S2b
    Q = a1 * e1 + A * s_b + c * Y.gsum(0, m, K / 2, m + W - 1)
    eta_hi = c * C_RS * (math.exp(-m * m * K / 2) + math.exp(-m * m * K / 2) * 0.5 * math.sqrt(math.pi / (m * K)))
    eta_lo = a1 * e1 + A * s_b
    d1_hi = c * C_RS * math.exp(-m * m * K / 2) * (m + 0.5 * m * math.sqrt(math.pi / (m * K)) + math.sqrt(math.pi) / 4 * (m * K) ** -1.5)
    s = min(1.0, M.QB / Q)
    S0_lo, S2_lo, S0_hi, S2_hi = (v * s * s for v in (S0_lo, S2_lo, S0_hi, S2_hi)); eta_lo, eta_hi, d1_hi = (v * s for v in (eta_lo, eta_hi, d1_hi))
    pr = M.profile(x, u, eta_t, Nc_dl=Nc_dl, uc=uc) if eta_t > 0 else M.profile(x, u, 0.0)
    return {'m': m, 'mc': mc, 'Nc': Nc, 'gap_ratio_m_over_2top': m / (2 * top), 'S0_lo': S0_lo, 'S2_lo': S2_lo, 'S0_hi': S0_hi, 'S2_hi': S2_hi,
            'eta_lo': eta_lo, 'eta_hi': eta_hi, 'sup_dR_hi': d1_hi, 'profile_check': pr}


def factors(ws):
    nRlo = math.sqrt(ws['S0_lo'] / 2)
    k0 = min(ws['eta_lo'], ws['eta_hi'] * nRlo / math.sqrt(ws['S0_hi'] / 2))
    k2 = min(ws['eta_lo'], ws['sup_dR_hi'] * nRlo / math.sqrt(ws['S2_hi'] / 2))
    f0 = ws['S0_lo'] / (ws['S0_lo'] + ws['S0_hi']); f2 = ws['S2_lo'] / (ws['S2_lo'] + ws['S2_hi'])
    l0, u0 = lemma_S_factors(ws['eta_lo'], ws['eta_hi'], f0, k0); l2, u2 = lemma_S_factors(ws['eta_lo'], ws['eta_hi'], f2, k2)
    return {'f_lo_S0': f0, 'f_lo_S2': f2, 'kappa_S0': k0, 'kappa_S2': k2, 'low_S0': l0, 'up_S0': u0, 'low_S2': l2, 'up_S2': u2}


def branches_split(lsig, x, S0, S2, F, lL2, ls2, lb):
    """by_eval.branches_log with separate factors per channel: lower factors in the numerators, upper factors in den (a lower bound of g)."""
    LN2, lae = Y.LN2, Y.lae
    lse2 = -lae(-lL2, LN2 - ls2); lgs = (lse2 - lL2) / 2; lV = 2 * lsig - lb
    Jww_n = lae(2 * lsig - lL2, 2 * lsig + math.log(S0 * F['low_S0'] / 2) + lgs + lse2 - 2 * ls2)
    Jww_d = lae(2 * lsig - lL2, 2 * lsig + math.log(S0 * F['up_S0'] / 2) + lgs + lse2 - 2 * ls2)
    Jvv_d = lae(lb, math.log(x * S2 * F['up_S2'] / 2) + lgs)
    lden = lae(Jvv_d, Jww_d)
    lN = math.log(0.5) + lae(lb + math.log(0.25) + lae(math.log(3) + lV, lL2), math.log(x * S2 * F['low_S2'] / 8) + lgs + lae(lV, lse2))
    return 2 * math.exp(Jww_n - lsig - lden) + math.exp(lN - lsig - lden)


def minimise_split(lsig, x, S0, S2, F):
    from scipy.optimize import minimize
    lL2max = -lsig; lbmax = math.log(x / 80); lbmin = 2 * lsig - 6 * math.log(10)
    def unpack(z):
        lL2 = min(z[0], lL2max); ls2 = min(max(z[1], 2 * lsig), lL2 - math.log(4)); lb = min(max(z[2], lbmin), lbmax); return lL2, ls2, lb
    best = None
    for dL in (0.5, 0.0, -2.0):
        for ds in (0.5, -1.0):
            for db in (3.0, 6.0):
                z0 = [0.5 * lsig + dL, lsig + ds, min(lsig + db, lbmax)]
                r = minimize(lambda z: branches_split(lsig, x, S0, S2, F, *unpack(z)), z0, method='Nelder-Mead', options={'xatol': 1e-9, 'fatol': 1e-14, 'maxiter': 2500})
                if best is None or r.fun < best.fun: best = r
    return best.fun


def counter_windows():
    rows = []
    for e, lx, u in ((48, -17.49, 3.92), (96, -40.88, 6.39), (300, -142.12, 12.46)):
        ell = e * math.log(10); lsig = -ell; x = 10 ** lx
        for et in (0.03, 0.1, 0.2):
            ws = window_split(x, u, et); pr = ws.pop('profile_check')
            S0 = ws['S0_lo'] + ws['S0_hi']; S2 = ws['S2_lo'] + ws['S2_hi']
            F = factors(ws)
            g_split = minimise_split(lsig, x, S0, S2, F) * math.sqrt(ell)
            old = {s: Y.minimise(lsig, x, pr['S0'], pr['S2'], pr['eta'], s)['g'] * math.sqrt(ell) for s in ('lo', 'mid', 'hi')}
            row = {'sig_exp': e, 'log10_x': lx, 'u_m': u, 'eta_target': et, 'eta_upper': ws['eta_lo'] + ws['eta_hi'], 'mu': x * S2 / S0,
                   'consistency_S0_rel': S0 / pr['S0'] - 1, 'consistency_S2_rel': S2 / pr['S2'] - 1, **ws, **F,
                   'gsl_old_lo': old['lo'], 'gsl_mid': old['mid'], 'gsl_hi': old['hi'], 'gsl_split': g_split}
            # Corollary S (scaling argument, see the report): g >= Lambda_S * (main), Lambda_S = (low_S2/up_S2) min(1, sqrt(low_S0*/low_S2)),
            # low_S0* = low_S0 - sigma up_S0 (the gamma sigma part of the S0 term of den); H1 with the raised target T_D needs Lambda_S T_D kappa_tau >= 1
            low0s = F['low_S0'] - 10.0 ** (-e) * F['up_S0']
            LamS = F['low_S2'] / F['up_S2'] * min(1.0, math.sqrt(low0s / F['low_S2']))
            TD = 1 + 2.1 / ell + 2.5e-4 + 10 / ell
            row.update({'Lambda_S': LamS, 'Lambda_old_(1-eta)/(1+eta)': (1 - row['eta_upper']) / (1 + row['eta_upper']), 'T_D10': TD,
                        'H1_split_product_Lambda_S_T_D_kappa_tau': LamS * TD * 0.99978})
            rows.append(row)
            print('Lambda_S %.5f  product %.5f' % (row['Lambda_S'], row['H1_split_product_Lambda_S_T_D_kappa_tau'])); print('1e-%d eta %.3f: split %.4f  (old lo %.4f mid %.4f hi %.4f)  f_lo S0 %.2e S2 %.2e kappa S0 %.2e S2 %.2e' % (
                e, row['eta_upper'], g_split, old['lo'], old['mid'], old['hi'], F['f_lo_S0'], F['f_lo_S2'], F['kappa_S0'], F['kappa_S2']), flush=True)
    return rows


if __name__ == '__main__':
    out = {'random_check': check_random()}
    print(out['random_check'], flush=True)
    out['counter_windows'] = counter_windows()
    json.dump(out, open('results/cd_split.json', 'w'), indent=1)
