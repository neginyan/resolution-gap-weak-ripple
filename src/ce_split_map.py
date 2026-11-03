"""Stage CE 2.3: counterexample search for Lemma S, the gap-free Lemma S', and the map of windows without a spectral gap.  Python 3.11.

LEMMA S' (no gap) [proof].  Same setting as Lemma S (cd_split.py) but WITHOUT 2M < m.  In the proof of Lemma S only the two terms
<Y_lo^2 r_hi> and <Y_lo Y_hi r_lo> used the gap; without it  |<Y_lo^2 r_hi>| <= eta_hi ||Y_lo||^2  and  |<Y_lo Y_hi r_lo>| <= ||Y_lo r_lo|| ||Y_hi||
<= eta_lo ||Y_lo|| ||Y_hi||  (Cauchy-Schwarz; 0 <= psi <= 1), so Lemma S holds with
    E' = E + eta_hi f_lo + 2 eta_lo sqrt(f_lo f_hi)          (B unchanged: its proof used Minkowski only).
With a gap, E' may be replaced by E.  [The overlap could be sharpened to ||(y_lo R_lo)_{>= m}|| instead of eta_lo ||y_lo||; not used here.]

(1) counterexample search for Lemma S (gap) and S' (no gap): random windows of both kinds (low blocks of 1..300 harmonics placed anywhere below m/2 for S, anywhere
    below m for S'; high blocks of 1..300 harmonics; coherent and random phases; total sup eta up to 0.95), then a hill climb on the worst ratio
    (exact <Y^2/h>/||Y||^2 on a 8192-point grid against the two bounds).  A NEGATIVE CONTROL: Lemma S's E (gap form) applied to gap-free windows must fail
    somewhere, which shows the search can find violations.
(2) the map: CC-type windows (Rudin-Shapiro band at u_m + a coherent block) with the block moved from u_c = 0.3 up to the band and its width varied
    (both sides: narrow/wide, low/high); for each, Lambda_S from Lemma S (gap) or S' (no gap), and the H1 product Lambda_S * T_10 * kappa_tau.
-> results/ce_split_map.json.
"""
import json, math, sys
import numpy as np
import by_eval as Y, cc_mid as M, cd_split as S

C_RS = 2 + math.sqrt(2)


def factors_general(eta_lo, eta_hi, f_lo, kappa, gap):
    eta = eta_lo + eta_hi; f_hi = 1 - f_lo
    E = eta_lo * f_lo + 2 * eta_hi * math.sqrt(max(f_lo * f_hi, 0)) + (eta_hi + kappa) * f_hi
    if not gap: E += eta_hi * f_lo + 2 * eta_lo * math.sqrt(max(f_lo * f_hi, 0))
    B = (eta_lo * math.sqrt(f_lo) + kappa * math.sqrt(f_hi) + eta_hi) ** 2
    return max(1 / (1 + eta), 1 - E), min(1 / (1 - eta), 1 + E + B / (1 - eta))


# ---------------- (1) search ----------------
def window_eval(nlo, alo, plo, nhi, ahi, phi, target, share, N=8192):
    t = np.arange(N) * 2 * np.pi / N
    def build(ns, a, p):
        if len(ns) == 0: return np.zeros(N), np.zeros(N)
        return (a[:, None] * np.cos(ns[:, None] * t + p[:, None])).sum(0), (-a[:, None] * ns[:, None] * np.sin(ns[:, None] * t + p[:, None])).sum(0)
    Rlo, dRlo = build(nlo, alo, plo); Rhi, dRhi = build(nhi, ahi, phi)
    slo = share * target / max(np.max(np.abs(Rlo)), 1e-300); shi = (1 - share) * target / max(np.max(np.abs(Rhi)), 1e-300)
    Rlo *= slo; dRlo *= slo; Rhi *= shi; dRhi *= shi
    R = Rlo + Rhi; elo = float(np.max(np.abs(Rlo))); ehi = float(np.max(np.abs(Rhi)))
    if elo + ehi >= 0.97: return None
    worst = []
    for psi in (1.0, 0.5):
        h = 1 + psi * R
        for ylo, yhi in ((Rlo, Rhi), (dRlo, dRhi)):
            y = ylo + yhi; nY = np.mean(y * y)
            if nY <= 0: continue
            ex = np.mean(y * y / h) / nY; f_lo = np.mean(ylo * ylo) / nY
            kap = min(elo, float(np.max(np.abs(yhi))) * math.sqrt(np.mean(Rlo * Rlo)) / max(math.sqrt(np.mean(yhi * yhi)), 1e-300))
            worst.append((ex, f_lo, kap))
    return elo, ehi, worst


def margins(res, gap_form):
    elo, ehi, worst = res; m = 9.0
    for ex, f_lo, kap in worst:
        lo, up = factors_general(elo, ehi, f_lo, kap, gap_form)
        m = min(m, ex - lo, up - ex)
    return m


def search(ntrial=600, seed=17, climb=40):
    rng = np.random.default_rng(seed); out = {'S_gap_windows': [9.0, 0], 'Sprime_nogap_windows': [9.0, 0], 'negative_control_S_on_nogap': [9.0, 0]}
    def rand_window(gap):
        M_ = int(rng.integers(1, 300)); m = int(rng.integers(2 * M_ + 1, 2 * M_ + 600)) if gap else int(rng.integers(M_ + 1, 2 * M_ + 1))
        Mh = min(m + int(rng.integers(1, 300)), 4000)
        start = int(rng.integers(1, max(2, M_)))
        nlo = np.arange(start, M_ + 1); nhi = np.arange(m, Mh + 1)
        coh = rng.random() < 0.5
        alo = np.ones(len(nlo)) if coh else rng.uniform(-1, 1, len(nlo)); plo = np.zeros(len(nlo)) if coh else rng.uniform(0, 2 * np.pi, len(nlo))
        ahi = rng.choice([-1.0, 1.0], len(nhi)) * (1 if rng.random() < 0.5 else rng.random(len(nhi))); phi = rng.uniform(0, 2 * np.pi, len(nhi))
        return [nlo, alo, plo, nhi, ahi, phi, rng.uniform(0.05, 0.95), rng.uniform(0, 1)]
    for gap in (True, False):
        for it in range(ntrial // 2):
            w = rand_window(gap); r = window_eval(*w)
            if r is None: continue
            keys = (['S_gap_windows'] if gap else ['Sprime_nogap_windows', 'negative_control_S_on_nogap'])
            for k in keys:
                form = gap if k != 'negative_control_S_on_nogap' else True
                m = margins(r, form)
                # hill climb on amplitudes of the low block and the split (adversarial)
                best = m; ww = [x.copy() if isinstance(x, np.ndarray) else x for x in w]
                for _ in range(climb if m < 0.05 else 5):
                    cand = [x.copy() if isinstance(x, np.ndarray) else x for x in ww]
                    cand[1] = cand[1] * np.exp(rng.normal(0, 0.3, len(cand[1]))); cand[7] = float(np.clip(cand[7] + rng.normal(0, 0.1), 0, 1)); cand[6] = float(np.clip(cand[6] + rng.normal(0, 0.05), 0.02, 0.95))
                    rc = window_eval(*cand)
                    if rc is None: continue
                    mc = margins(rc, form)
                    if mc < best: best, ww = mc, cand
                out[k][0] = min(out[k][0], best); out[k][1] += 1
    return {k: {'min_margin': v[0], 'windows': v[1], 'violations_found': v[0] < -1e-12} for k, v in out.items()}


# ---------------- (2) the map ----------------
def window_general(x, u, eta_t, uc, Nc_dl):
    K = x; m = max(3, int(round(u / math.sqrt(x)))); p = Y.profile_rs(x, m); c = p['c']; W = p['W']
    mc = max(2, int(round(uc / math.sqrt(x)))); Nc = max(1, int(round(Nc_dl / (mc * K)))); Nc = min(Nc, m - mc)
    if Nc < 1: return None
    top = mc + Nc - 1; gap = 2 * top < m
    s_b = Y.gsum(0, mc, K / 2, top); S0b = Y.gsum(0, mc, K, top); S2b = Y.gsum(2, mc, K, top)
    A = max(0.0, eta_t - p['eta']) / s_b if s_b > 0 else 0.0
    if c * (1 + math.sqrt(2 * W)) + A * (1 + Nc) > 1: return None
    a1 = max(c, A); e1 = math.exp(-K / 2)
    S0_hi = c * c * Y.gsum(0, m, K, m + W - 1); S2_hi = c * c * Y.gsum(2, m, K, m + W - 1)
    S0_lo = (a1 * e1) ** 2 + A * A * S0b; S2_lo = (a1 * e1) ** 2 + A * A * S2b
    Q = a1 * e1 + A * s_b + c * Y.gsum(0, m, K / 2, m + W - 1)
    eta_hi = c * C_RS * (math.exp(-m * m * K / 2) * (1 + 0.5 * math.sqrt(math.pi / (m * K))))
    eta_lo = a1 * e1 + A * s_b
    d1_hi = c * C_RS * math.exp(-m * m * K / 2) * (m + 0.5 * m * math.sqrt(math.pi / (m * K)) + math.sqrt(math.pi) / 4 * (m * K) ** -1.5)
    sc = min(1.0, M.QB / Q)
    S0_lo, S2_lo, S0_hi, S2_hi = (v * sc * sc for v in (S0_lo, S2_lo, S0_hi, S2_hi)); eta_lo, eta_hi, d1_hi = (v * sc for v in (eta_lo, eta_hi, d1_hi))
    return {'gap': gap, 'S0_lo': S0_lo, 'S2_lo': S2_lo, 'S0_hi': S0_hi, 'S2_hi': S2_hi, 'eta_lo': eta_lo, 'eta_hi': eta_hi, 'sup_dR_hi': d1_hi, 'mu': x * (S2_lo + S2_hi) / (S0_lo + S0_hi)}


def lam_S(ws, sig):
    nRlo = math.sqrt(ws['S0_lo'] / 2)
    k0 = min(ws['eta_lo'], ws['eta_hi'] * nRlo / math.sqrt(ws['S0_hi'] / 2)); k2 = min(ws['eta_lo'], ws['sup_dR_hi'] * nRlo / math.sqrt(ws['S2_hi'] / 2))
    f0 = ws['S0_lo'] / (ws['S0_lo'] + ws['S0_hi']); f2 = ws['S2_lo'] / (ws['S2_lo'] + ws['S2_hi'])
    l0, u0 = factors_general(ws['eta_lo'], ws['eta_hi'], f0, k0, ws['gap']); l2, u2 = factors_general(ws['eta_lo'], ws['eta_hi'], f2, k2, ws['gap'])
    low0s = l0 - sig * u0
    return l2 / u2 * min(1.0, math.sqrt(max(low0s, 0) / l2)), (l0, u0, l2, u2, f0, f2)


def the_map():
    rows = []
    for e, lx, u in ((48, -17.49, 3.92), (96, -40.88, 6.39), (300, -142.12, 12.46)):
        ell = e * math.log(10); x = 10 ** lx; TD = 1 + 12.1 / ell + 2.5e-4
        for eta_t in (0.03, 0.1, 0.2):
            for uc_frac in (0.03, 0.1, 0.2, 0.35, 0.45, 0.6, 0.8, 0.95):
                for Nc_dl in (0.005, 0.05, 0.2, 1.0, 5.0):
                    ws = window_general(x, u, eta_t, uc_frac * u, Nc_dl)
                    if ws is None: continue
                    lam, f = lam_S(ws, 10.0 ** (-e))
                    rows.append({'sig_exp': e, 'eta_target': eta_t, 'eta': ws['eta_lo'] + ws['eta_hi'], 'u_block_over_u_band': uc_frac, 'block_width_decay_lengths': Nc_dl,
                                 'gap': ws['gap'], 'mu': ws['mu'], 'Lambda_S': lam, 'old_(1-eta)/(1+eta)': (1 - ws['eta_lo'] - ws['eta_hi']) / (1 + ws['eta_lo'] + ws['eta_hi']),
                                 'product_with_T10_kappa': lam * TD * 0.99978, 'covered_by_corollary_S': lam * TD * 0.99978 >= 1,
                                 'covered_by_old_eta_class': (ws['eta_lo'] + ws['eta_hi']) <= 1 / ell, 'f_lo_S0': f[4], 'f_lo_S2': f[5]})
    return rows


if __name__ == '__main__':
    res = {'search': search()}
    print(res['search'], flush=True)
    res['map'] = the_map()
    n = len(res['map']); cov = sum(r['covered_by_corollary_S'] for r in res['map'])
    res['map_summary'] = {'windows': n, 'covered': cov, 'with_gap': sum(r['gap'] for r in res['map']),
                          'covered_with_gap': sum(r['covered_by_corollary_S'] for r in res['map'] if r['gap']),
                          'covered_without_gap': sum(r['covered_by_corollary_S'] for r in res['map'] if not r['gap'])}
    print(res['map_summary'], flush=True)
    json.dump(res, open('results/ce_split_map.json', 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
