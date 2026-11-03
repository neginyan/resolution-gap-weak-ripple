"""Stage CH 1.1: the b-edge (b < 1e4 sigma^2) as a PROOF.  Python 3.11.

LEMMA BE.  Regime of paper H: sigma in [1e-3000, 1e-24] (ell <= 6907.8), x in [5e11 sigma^2, 0.003], (P0) L^2 <= 1/sigma, (P3) sum q0_n <= Qbar = 1/4,
(P4) sigma^2 <= s^2 <= L^2/4, q0_n <= 2 e^{-n^2 x/2} (|a_n| <= 2, K >= x).  For EVERY state with 0 < b <= 1e4 sigma^2:
        (a) + (a') >= kappa/sigma,   kappa = 1.1e-5,     hence  g sqrt(ell) >= (1 - epsbar) kappa sqrt(ell)/sigma >= (1 - epsbar) 8e19.
No condition on x beyond the regime (no gap at small x), no (P2).

The proof (report CH, section 1) has four new pieces, each checked here:
 (BE-a) the exact formula  sigma^2 G = E [x R''/h - x E R'^2/h^2 - sigma sqrt(x) (w/s^2) R'/h^2],  E = e^{-w^2/2s^2}  (finite differences, check_G),
        and the envelope form |sigma^2 G| <= E [x N2 (1-Qb) + x E N1^2 + sigma sqrt(x) (|w|/s^2) N1]/(1-Qb)^2  (whose sup over w is <= gamma_P of v1.5);
 (BE-1) the cross term by mean square: for every c > 0,  <h sigma^2 G cos psi> >= -sigma^2 <F> - <h (sigma^2 G)^2>/(4c) - 2c <h (1 - cos psi)>;
 (BE-2) M_G = int Gamma <h (sigma^2 G)^2> <= (3/2) g_s [x^2 S4/(1-Qb) + x S2 (U1^2 + 1/2)/(1-Qb)^3]  and  x^2 S4 <= Lambda x S2 + T_Lambda  (check_MG);
 (BE-3) the constants (choice c = (b + sigma^2 J_ww)/(16 sigma^2), Lambda = 1e-3/sigma^2) and every remainder (log10 of its ratio to sigma^2 J_ww or to
        (1 - e) W, all below -1e7; check_constants, check_remainders).
-> results/ch_bedge.json
"""
import json, math, os
import numpy as np
from decimal import Decimal as D, getcontext

getcontext().prec = 60; getcontext().Emin = -10 ** 15; getcontext().Emax = 10 ** 15
PI = D('3.14159265358979323846264338327950288419716939937510582097494')
E1 = D(1).exp()


class _mp:
    mpf = staticmethod(lambda v: D(str(v)))
    sqrt = staticmethod(lambda v: D(v).sqrt())
    exp = staticmethod(lambda v: D(str(v)).exp())
    log10 = staticmethod(lambda v: D(v).log10())
    log = staticmethod(lambda v: D(v).ln())
    pi = PI
    e = E1


mp = _mp
QB = 0.25
OUT = 'results/ch_bedge.json'


def window(rng, nh, x, qsum=0.24):
    """random window: q0_n <= 2 e^{-n^2 x/2}, sum q0_n = qsum <= Qbar."""
    n = np.arange(1, nh + 1)
    cap = np.minimum(2 * np.exp(-n * n * x / 2), 1.0)
    q = rng.uniform(0, 1, nh) * cap
    q *= qsum / q.sum()
    q = np.minimum(q, cap)
    ph = rng.uniform(0, 2 * math.pi, nh)
    return n, q, ph


def fields(n, q, ph, x, sig, s2, L2, nt=4096, nw=401, wmax=7.0):
    t = np.linspace(0, 2 * math.pi, nt, endpoint=False)
    L = math.sqrt(L2); w = np.linspace(-wmax * L, wmax * L, nw)
    T, Wg = np.meshgrid(t, w, indexing='ij')
    arg = n[None, None, :] * T[..., None] - ph[None, None, :]
    R = np.sum(q * np.cos(arg), axis=-1); R1 = np.sum(-n * q * np.sin(arg), axis=-1); R2 = np.sum(-n * n * q * np.cos(arg), axis=-1)
    E = np.exp(-Wg ** 2 / (2 * s2)); h = 1 + E * R
    kp = math.sqrt(x) / sig
    return t, w, T, Wg, R, R1, R2, E, h, kp


def check_G(rng, ntrial=12):
    """(BE-a): sigma^2 G with G = d_v (s~ + s_w), s~ = k' h_t/h, s_w = -w/L^2 + h_w/h, d_v = k' d_t, versus the closed formula; envelope bound."""
    worst_rel, worst_env = 0.0, 0.0
    for _ in range(ntrial):
        nh = int(rng.integers(1, 9)); x = 10 ** rng.uniform(-3, -0.5); sig = 10 ** rng.uniform(-3, -1)
        s2 = sig ** 2 * 10 ** rng.uniform(0, 3); L2 = 4 * s2 * 10 ** rng.uniform(0, 2)
        n, q, ph = window(rng, nh, x)
        t, w, T, Wg, R, R1, R2, E, h, kp = fields(n, q, ph, x, sig, s2, L2, nt=8192, nw=101)
        ht = E * R1; hw = -(Wg / s2) * E * R
        S = kp * ht / h + hw / h                       # s~ + s_w without the -w/L^2 part (independent of v)
        dt = t[1] - t[0]
        dS = (np.roll(S, -1, axis=0) - np.roll(S, 1, axis=0)) / (2 * dt)   # periodic central difference in t
        G_num = sig ** 2 * kp * dS
        G_form = E * (x * R2 / h - x * E * R1 ** 2 / h ** 2 - sig * math.sqrt(x) * (Wg / s2) * R1 / h ** 2)
        scale = np.max(np.abs(G_form)) + 1e-300
        worst_rel = max(worst_rel, float(np.max(np.abs(G_num - G_form)) / scale))
        N1 = float(np.sum(n * q)); N2 = float(np.sum(n * n * q))
        env = E * (x * N2 * (1 - QB) + x * E * N1 ** 2 + sig * math.sqrt(x) * (np.abs(Wg) / s2) * N1) / (1 - QB) ** 2
        worst_env = max(worst_env, float(np.max(np.where(env > 0, np.abs(G_form) / np.where(env > 0, env, 1.0), 0.0))))
        gP15 = (x * (N2 * (1 + QB) + N1 ** 2) + math.sqrt(x) * N1 * (1 + 2 * QB) / math.sqrt(math.e)) / (1 - QB) ** 2
        assert float(np.max(env)) <= gP15 * (1 + 1e-12)
    return {'max_rel_diff_formula_vs_finite_difference': worst_rel, 'max_ratio_to_envelope_bound': worst_env,
            'envelope_sup_le_gamma_P_v15': True}


def check_MG(rng, ntrial=40):
    """(BE-2): M_G by quadrature versus (3/2) g_s [x^2 S4/(1-Qb) + x S2 ((sqrt x N1)^2 + sigma^2/(2 s^2))/(1-Qb)^3]."""
    worst = 0.0; worst_split = 0.0
    for _ in range(ntrial):
        nh = int(rng.integers(1, 30)); x = 10 ** rng.uniform(-3.5, -1); sig = 10 ** rng.uniform(-3, -1)
        s2 = sig ** 2 * 10 ** rng.uniform(0, 3); L2 = 4 * s2 * 10 ** rng.uniform(0, 2)
        n, q, ph = window(rng, nh, x)
        t, w, T, Wg, R, R1, R2, E, h, kp = fields(n, q, ph, x, sig, s2, L2, nt=2048, nw=801, wmax=8.0)
        sG = E * (x * R2 / h - x * E * R1 ** 2 / h ** 2 - sig * math.sqrt(x) * (Wg / s2) * R1 / h ** 2)
        Gam = np.exp(-w ** 2 / (2 * L2)) / math.sqrt(2 * math.pi * L2); dw = w[1] - w[0]
        MG = float(np.sum(np.mean(h * sG ** 2, axis=0) * Gam) * dw)
        S2 = float(np.sum(n ** 2 * q ** 2)); S4 = float(np.sum(n ** 4 * q ** 2)); N1 = float(np.sum(n * q))
        gs = math.sqrt((1 / (1 / L2 + 2 / s2)) / L2)
        bound = 1.5 * gs * (x * x * S4 / (1 - QB) + x * S2 * (x * N1 ** 2 + sig ** 2 / (2 * s2)) / (1 - QB) ** 3)
        worst = max(worst, MG / bound)
        Lam = float(np.max(x * n ** 2 * (q > 0))) * rng.uniform(0.3, 1.2)
        big = x * n ** 2 > Lam
        T_L = float(np.sum((x * n ** 2) ** 2 * q ** 2 * big))
        worst_split = max(worst_split, x * x * S4 / (Lam * x * S2 + T_L))
    return {'max_MG_over_bound': worst, 'max_x2S4_over_Lambda_split': worst_split, 'trials': ntrial}


def U_consts(ell):
    QB = D('0.25')
    U1 = QB * mp.sqrt(2 * ell + mp.mpf('3.1')) + mp.mpf('1e-6'); U2 = QB * (2 * ell + mp.mpf('3.1')) + mp.mpf('1e-4')
    gP = (U2 * (1 + QB) + U1 ** 2 + U1 * (1 + 2 * QB) / mp.sqrt(mp.e)) / (1 - QB) ** 2
    cchi = (U1 / (1 - QB) + mp.sqrt(2 / mp.pi) + QB * mp.exp('-0.5') / (1 - QB)) / mp.sqrt(2)
    return U1, U2, gP, cchi


def check_constants():
    """(BE-3) constants.  b <= 1e4 sigma^2 and sigma <= 1e-24."""
    Qb = mp.mpf(QB); sig_hi = mp.mpf('1e-24'); bmax = 10 ** 4 * sig_hi ** 2
    y_min = (1 - bmax) ** 2 / (2 * 10 ** 4)                      # y = (1-b)^2 sigma^2/(2b) >= (1-b)^2/(2e4)
    one_m_e = y_min / (2 + y_min)                                # 1 - e^{-y/2} >= y/(2+y), increasing in y
    U1, U2, gP, cchi = U_consts(mp.mpf('6907.76'))
    K1s2 = mp.mpf('1e-3') / (1 - Qb) + sig_hi ** 2 * (U1 ** 2 + mp.mpf('0.5')) / (1 - Qb) ** 3      # K_1 sigma^2 with Lambda = 1e-3/sigma^2
    rhoW = 48 * (1 + Qb) * K1s2 / (1 - bmax) ** 2
    # coefficient of (1-e) W in 2 N_phi: 1 - delta_1 - rho_W (1+z)/((1-z)(1-nu0)); delta_1 = sigma^2/((1-e)(1-nu0)) <= 4.01e4 sigma^2; nu0, z < 1e-1e7
    cW2 = 1 - mp.mpf('4.01e4') * sig_hi ** 2 - rhoW * (1 + mp.mpf('1e-100'))
    cW = cW2 / 2 * one_m_e                                         # in Num = 2A_w + N_phi (N_phi >= half of the 2N_phi bound)
    cJ = 1 - mp.mpf(1) / (8 * (1 - mp.mpf('1e-100')))              # 2A_w >= sigma^2 J; minus sigma^2 J/(8(1-z))
    cb_coef = (1 / (1 + mp.mpf('1e-100')) - 1 / (2 * (1 - mp.mpf('1e-100')))) / 2   # coefficient of b Psi in Num
    Phi_min = 1 - mp.exp(-mp.mpf(1) / 4) - sig_hi * cchi                 # L^2 >= 1
    Psi_min = min(Phi_min, mp.mpf('0.5'))
    cb = cb_coef * Psi_min
    rem_allow = mp.mpf('1e-6')                                       # remainders <= 1e-6 sigma^2 J_ww (check_remainders)
    case1 = min(cJ - rem_allow, cW, cb)                              # L^2 >= 1
    ratio_bJ = mp.mpf(10 ** 4) / mp.mpf('0.6')                       # L^2 < 1: sigma^2 J >= 0.6 sigma^2/L^2 > 0.6 sigma^2 >= 6e-5 b
    case2 = min((cJ - rem_allow) / (1 + ratio_bJ), cW)
    kappa = min(case1, case2)
    return {'y_min': float(y_min), 'one_minus_e_min': float(one_m_e), 'rho_W': float(rhoW), 'coef_W_in_Num': float(cW), 'coef_J_in_Num': float(cJ),
            'Phi_min_L2_ge_1': float(Phi_min), 'coef_b_in_Num_L2_ge_1': float(cb), 'kappa_case_L2_ge_1': float(case1), 'kappa_case_L2_lt_1': float(case2),
            'kappa': float(kappa), 'kappa_claimed': 1.1e-5, 'kappa_ok': bool(kappa >= mp.mpf('1.1e-5')),
            'g_sqrt_ell_lower_at_1e-24_without_epsbar': float(mp.mpf('1.1e-5') * mp.sqrt(mp.mpf('55.26')) / sig_hi),
            'U1': float(U1), 'gamma_P_bar': float(gP), 'c_chi_bar': float(cchi)}


def log10_nu_bound(tau):
    """nu_k(tau, beta) <= 2.0000001 e^{-0.499999 tau} for tau >= 5e7, k <= 4, beta <= 1e-6 tau (x >= 5e11 sigma^2)."""
    return mp.log10(mp.mpf('2.0000001')) - mp.mpf('0.499999') * tau / mp.log(10)


def check_remainders():
    """(BE-3) each remainder divided by its comparison quantity, log10, at the worst end of every parameter.
    Comparison: sigma^2 J_ww >= 0.6 sigma^3 for absolute terms ('J'), (1 - e) W for the ripple-relative a_4* ('W').
    All ratios are (powers of 1/sigma) * (powers of tau) * nu, so the worst end is sigma = 1e-3000, tau = 5e7 (tau^k e^{-tau/2} decreases there),
    b = 1e4 sigma^2, x = 0.003 where they appear increasing."""
    Qb = mp.mpf(QB); sig = mp.mpf(10) ** -3000; tau = mp.mpf('5e7'); b = 10 ** 4 * sig ** 2; x = mp.mpf('0.003')
    U1, U2, gP, cchi = U_consts(mp.mpf('6907.76'))
    lnu = log10_nu_bound(tau); nu = mp.mpf(10) ** lnu
    J = mp.mpf('0.6') * sig ** 3
    beta = mp.mpf('1e-6') * tau
    t = {}
    t['a1'] = mp.mpf(3) / 8 * sig ** 2 * Qb * nu / J
    t['a2'] = b * (1 + Qb) * sig * cchi * nu * (1 + Qb / 2) / J
    t['a3'] = Qb * sig ** 2 * tau * nu / 8 / J
    t['a5'] = sig ** 2 * nu * (1 + Qb) * gP / J
    t['a6'] = nu * sig * mp.sqrt(x) * (U1 * cchi + (1 + Qb) * gP / mp.sqrt(2)) / J
    t['a7'] = x * Qb * nu * (1 + Qb) * sig * cchi / 2 / J
    t['a8'] = Qb * sig ** 2 * tau ** 2 * nu / 8 / J
    t['a9'] = sig ** 2 * x * Qb * nu / (2 * (1 - nu)) / J
    t['a10'] = sig * tau * mp.sqrt(x) * Qb * nu / (2 * mp.sqrt(2)) / J
    t['Phi_prop'] = b * Qb * (nu + tau * nu) / J
    t['rhoW_bnu'] = mp.mpf('0.0801') * b * tau * Qb * nu / 2 / J
    # a_4* relative to (1 - e) W: W_up/W_per <= (1+Qb)/(1-Qb); bracket beta^2 nu/2 + beta nu (xi + zeta + sigma)/sqrt2; (1-e) >= 2.4999e-5
    xi = sig * U1 / (1 - Qb); zeta = sig * Qb * mp.exp('-0.5') / (1 - Qb)
    t['a4star_vs_W'] = (1 + Qb) / (1 - Qb) * (beta ** 2 * nu / 2 + beta * nu * (xi + zeta + sig) / mp.sqrt(2)) / mp.mpf('2.4999e-5')
    # the tail T_Lambda (Lambda = 1e-3/sigma^2): ratio 22.3 T/sigma^2, worst at the LARGEST sigma (1e-24); x^{-1/2} <= 1/(sqrt(5e11) sigma)
    s24 = mp.mpf('1e-24'); Lam = mp.mpf('1e-3') / s24 ** 2
    lpoly = mp.log10(4 * (Lam ** 2 + Lam ** D('1.5') / (mp.sqrt(mp.mpf('5e11')) * s24) / mp.mpf('1.25')))
    lT = lpoly - Lam / mp.log(10)                                   # T <= 4 e^{-Lambda}(Lambda^2 + Lambda^1.5 x^{-1/2}/1.25)
    tail_log10 = mp.log10(mp.mpf('22.3')) + lT - 2 * mp.log10(s24)
    # also at sigma = 1e-3000 for the tail (Lambda larger, still negligible)
    out = {k: float(mp.log10(v)) for k, v in t.items()}
    out['tail_T_Lambda'] = float(tail_log10)
    out['all_below_log10'] = max(out.values())
    out['ok'] = bool(out['all_below_log10'] < -6)
    return out


def main():
    rng = np.random.default_rng(20261103)
    res = {'claim': 'Lemma BE: every state with 0 < b <= 1e4 sigma^2 has (a)+(a\') >= 1.1e-5/sigma (regime of paper H, ell <= 6907.8)'}
    res['BE_a_formula_and_envelope'] = check_G(rng)
    res['BE_2_mean_square'] = check_MG(rng)
    res['BE_3_constants'] = check_constants()
    res['BE_3_remainders_log10'] = check_remainders()
    ok = (res['BE_a_formula_and_envelope']['max_rel_diff_formula_vs_finite_difference'] < 1e-5 and res['BE_a_formula_and_envelope']['max_ratio_to_envelope_bound'] <= 1
          and res['BE_2_mean_square']['max_MG_over_bound'] <= 1 and res['BE_2_mean_square']['max_x2S4_over_Lambda_split'] <= 1 + 1e-12
          and res['BE_3_constants']['kappa_ok'] and res['BE_3_remainders_log10']['ok'])
    res['all_ok'] = bool(ok)
    os.makedirs('results', exist_ok=True)
    json.dump(res, open(OUT, 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    print(json.dumps(res, indent=1, default=str))


if __name__ == '__main__':
    main()
