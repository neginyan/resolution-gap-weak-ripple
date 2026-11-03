"""Stage CA section 2.1: the reduction as a lemma with explicit error terms (PROOF), and Lemma Z1 with those errors (PROOF), checked against direct
two-dimensional quadrature (NUMERICAL check of the proof) and evaluated at the hardest band points of BZ (NUMERICAL).

Lemma R (lower side of the branches, every periodic window state, every tau; built from v1.5 Lemmas BL0-BL3 (proved) and three new steps).
Notation of App. D: e_n = e^{-n^2 tau/2}, alpha = (1 - b)/sqrt2, V = sigma^2/b, y = alpha^2 V, e = e^{-y/2}, k' = sqrt(x)/sigma, beta = alpha k' V,
nu_k = 2 sum n^k e_n e^{n beta}, omega_2 = 2 sum e_n |1 - n^2 tau|, z = Qbar nu_0/2, eps_1 = nu_0 (1 + 2 Qbar)/(1 - z).  Window: h - 1 = e^{-w^2/2s^2} R(t),
eta = sup|R| (< 1), S0 = sum q0_n^2, S2 = sum n^2 q0_n^2, P = x S2, N_1, N_2 the v1.5 sums (only in remainders).  Gaussian integrals under Gamma = N(0, L^2):
se2 = 1/(1/L^2 + 2/s^2), g_s = sqrt(se2/L^2), int Gamma e^{-w^2/s^2} {1, w^2, w^4, |w|} = g_s {1, se2, 3 se2^2, sqrt(2 se2/pi)},
int Gamma e^{-w^2/s^2}(1 - cos(w/sqrt2)) = g_s (1 - e^{-se2/4}),  D0 = S0 g_s se2/(2 s^4).
  (R1) sandwich [proof]: 1 - eta <= h <= 1 + eta, so J_per := 1/L^2 + int Gamma <h_w^2/h> lies in [1/L^2 + D0/(1 + eta), 1/L^2 + D0/(1 - eta)] and
       W_per := x int Gamma <h_t^2/h> in [P g_s/(2(1 + eta)), P g_s/(2(1 - eta))]  (Parseval: <h_t^2> = S2 e^{-w^2/s^2}/2, <h_w^2> = (w^2/s^4) S0 e^{-w^2/s^2}/2).
  (R2) den <= b + [(1 + nu_0) W_per,up + b tau Qbar nu_2/2]/(1 - z) + sigma^2 (1 + eps_1) J_per,up            [BL1, BL2].
  (R3) A_w >= sigma^2 J_ww - (sigma^2/4) E[s_w^2 phi] with sigma^2 J_ww >= (1 - eps_1) sigma^2 J_per,lo [BL1] and [new, proof]
       E[s_w^2 phi] <= (1/2) E[s_w^2 v'^2] + E[s_w^2 w^2] + 2 E[s_w^2 chi^2]  (phi <= m_z^2/2 <= (v'^2 + R^2)/2, R^2 <= 2 w^2 + 4 chi^2, v' = (1 - b) v),
       E[s_w^2 v'^2] <= (1 - b)^2 V (1 + omega_2) J_per,up/(1 - z)   (BL0(iii) with F = h s_w^2 >= 0 and <h s_w^2> = w^2/L^4 + <h_w^2/h>),
       E[s_w^2 w^2] <= (1 + nu_0)/(1 - z) [3 + 3 se2 D0/(1 - eta)],
       E[s_w^2 chi^2] <= (1 + nu_0)/(1 - z) (3/2)[(xi^2 + zeta^2) J_per,up + (sigma^4/L^4)(3 + 3 se2 D0/(1 - eta))],
       chi = sigma^2(s~ + s_w)/sqrt2, |chi| <= chibar(w) = (xi + zeta + sigma^2 |w|/L^2)/sqrt2, xi = sigma sqrt(x) N_1/(1 - Qbar), zeta = sigma^2 eta e^{-1/2}/(s(1 - eta)).
  (R4) 2 N_phi >= [BL3 right side] + e (1 - nu_0)/(1 + z) * Fterm   [new, proof]: BL3 drops e sigma^2 <theta_tau F (1 - cos psi)> >= 0 from the exact identity
       (eq. BLid); keeping it, with 1 - cos(A + d) >= 1 - cos A - |A||d| - d^2/2 (A = w/sqrt2, d = chi), F = k'^2 h_t^2/h and the sandwich,
       Fterm >= P g_s (1 - e^{-se2/4})/(2(1 + eta)) - P g_s [(c0 sqrt(2 se2/pi) + c1 se2)/sqrt2 + c0^2 + c1^2 se2]/(2(1 - eta)),
       c0 = (xi + zeta)/sqrt2, c1 = sigma^2/(sqrt2 L^2).
       In BL3: W >= [(1 - nu_0) W_per,lo - b tau Qbar nu_2/2]/(1 + z), Psi >= 1 - e(1 - y)(e^{-L^2/4} + r_C), Phi <= 1 - e^{-L^2/4} + r_C,
       r_C = (xi + zeta + sigma^2 sqrt(2/pi)/L)/sqrt2 (|cos psi - cos(w/sqrt2)| <= |chi|, <h> = 1, int Gamma cos(w/sqrt2) = e^{-L^2/4}); gamma_P, E_R (a_1..a_10) as in BL3.
  Then (a) + (a') = (2 A_w + N_phi)/(sigma den) >= Num_lo/(sigma den_up), and g >= (1 - |eps|)((a) + (a')) (Tool 1).
Lemma Z1R (proof): split Num_lo = n_A + n_B and den_up = d_A + d_B (ripple parts: everything proportional to D0 or P; background parts: the rest); the mediant
inequality gives (a) + (a') >= min(T_A, T_B), T_X = n_X/(sigma d_X), and by AM-GM in s^2 (as in BZ's Z1)
  T_B >= (1 - eta)/(1 + eta) * (se2/s^2) sqrt(k1 k2)/( sqrt(mu) (c_0 + c_1 A/(2 k1)) ),  A = (se2/s^2) sigma sqrt(k1 k2/mu),
with k1 = 1 - eps_1 - (ripple part of (sigma^2/2) E[s_w^2 phi])/(2 sigma^2 D0/(1 + eta)), k2 = e (1 - nu_0)/(1 + z) (1 - se2/8 - rho_F),
c_0 = (1 + nu_0)/(1 - z) (1 + eta)... [all explicit in the code], i.e. Z1 of BZ times a factor 1 - delta_B with delta_B explicit.
-> results/ca_reduction.json.  Python 3.11."""
import json, math, sys
import numpy as np
from ca_dec import D, exp, log10, sqrt, expm1, log1p, PI, E, INF

SQ2 = sqrt(2)


def series(tau, beta, kmax=6, nmax=60):
    """nu_k (with e^{n beta}), omega_2; terms beyond nmax are below e^{-tau nmax^2/2 + ...} (tau >= 0.5 here)."""
    nu = [D(0)] * (kmax + 1); om2 = D(0)
    for n in range(1, nmax + 1):
        if n * n * tau / 2 - n * beta > D('1e12'): break                  # the term is below e^{-1e12}
        en = exp(-n * n * tau / 2)
        for k in range(kmax + 1): nu[k] += 2 * n ** k * en * exp(n * beta)
        om2 += 2 * en * abs(1 - n * n * tau)
    return nu, om2


def lemma_R(sig, x, b, L2, s2, S0, S2, eta, Qb, N1, N2, detail=False):
    """Lemma R and Z1R; all inputs mp numbers (sig may be 1e-500).  Returns Num_lo, den_up, T_A, T_B, the Z1R corollary and the pieces."""
    sig, x, b, L2, s2, S0, S2, eta, Qb, N1, N2 = (D(t) for t in (sig, x, b, L2, s2, S0, S2, eta, Qb, N1, N2))
    tau = x / b; V = sig ** 2 / b; al = (1 - b) / SQ2; y = al ** 2 * V; e = exp(-y / 2); kp = sqrt(x) / sig; beta = al * kp * V
    nu, om2 = series(tau, beta)
    z = Qb * nu[0] / 2; eps1 = nu[0] * (1 + 2 * Qb) / (1 - z)
    se2 = 1 / (1 / L2 + 2 / s2); gs = sqrt(se2 / L2); s = sqrt(s2); L = sqrt(L2)
    D0 = S0 * gs * se2 / (2 * s2 ** 2); P = x * S2
    Jlo_G, Jlo_R = 1 / L2, D0 / (1 + eta); Jup_G, Jup_R = 1 / L2, D0 / (1 - eta)
    Wper_lo = P * gs / (2 * (1 + eta)); Wper_up = P * gs / (2 * (1 - eta))
    # (R2) den, split background / ripple
    dA = b + b * tau * Qb * nu[2] / 2 / (1 - z) + sig ** 2 * (1 + eps1) * Jup_G
    dB = (1 + nu[0]) * Wper_up / (1 - z) + sig ** 2 * (1 + eps1) * Jup_R
    # pointwise sizes
    xi = sig * sqrt(x) * N1 / (1 - Qb); zeta = sig ** 2 * eta * exp(D('-0.5')) / (s * (1 - eta))
    # (R3) (sigma^2/2) E[s_w^2 phi] <= (sigma^2/2)[ (1/2) E v'^2 + E w^2 + 2 E chi^2 ], split G / R
    cz = 1 / (1 - z)
    Ev_G = (1 - b) ** 2 * V * (1 + om2) * Jup_G * cz; Ev_R = (1 - b) ** 2 * V * (1 + om2) * Jup_R * cz
    Ew_G = (1 + nu[0]) * cz * 3; Ew_R = (1 + nu[0]) * cz * 3 * se2 * D0 / (1 - eta)
    Ec_G = (1 + nu[0]) * cz * D('1.5') * ((xi ** 2 + zeta ** 2) * Jup_G + sig ** 4 / L2 ** 2 * 3)
    Ec_R = (1 + nu[0]) * cz * D('1.5') * ((xi ** 2 + zeta ** 2) * Jup_R + sig ** 4 / L2 ** 2 * 3 * se2 * D0 / (1 - eta))
    corrA_G = sig ** 2 / 2 * (Ev_G / 2 + Ew_G + 2 * Ec_G); corrA_R = sig ** 2 / 2 * (Ev_R / 2 + Ew_R + 2 * Ec_R)
    twoAw_G = 2 * (1 - eps1) * sig ** 2 * Jlo_G - corrA_G; twoAw_R = 2 * (1 - eps1) * sig ** 2 * Jlo_R - corrA_R
    # (R4) BL3 pieces
    gP = (x * (N2 * (1 + Qb) + N1 ** 2) + sqrt(x) * N1 * (1 + 2 * Qb) / sqrt(E)) / (1 - Qb) ** 2
    cchi = (sqrt(x) * N1 / (1 - Qb) + sqrt(2 / PI) + Qb * exp(D('-0.5')) / (1 - Qb)) / SQ2
    a = [D(3) / 8 * sig ** 2 * Qb * nu[0], b * (1 + Qb) * sig * cchi * nu[0] * (1 + Qb / 2), b * Qb * sig ** 2 * tau ** 2 * nu[2] / (8 * x),
         sig * tau * sqrt(x) * N1 ** 2 * nu[1] / (SQ2 * (1 - Qb)), sig ** 2 * nu[0] * (1 + Qb) * gP,
         nu[1] * (x * N1 * sig * cchi + sig * sqrt(x) * (1 + Qb) * gP / SQ2), x * Qb * nu[2] * (1 + Qb) * sig * cchi / 2,
         Qb * sig ** 2 * tau ** 2 * nu[4] / 8, sig ** 2 * x * Qb * nu[2] / (2 * (1 - nu[0])), sig * tau * sqrt(x) * Qb * nu[3] / (2 * SQ2)]
    sbER = sum(a)
    rC = (xi + zeta + sig ** 2 * sqrt(2 / PI) / L) / SQ2
    eL = exp(-L2 / 4)
    Psi_lo = -expm1(log1p(-y) - y / 2 - L2 / 4) - e * (1 - y) * rC; Phi_up = -expm1(-L2 / 4) + rC    # no cancellation (sigma down to 1e-500)
    Wlo = ((1 - nu[0]) * Wper_lo - b * tau * Qb * nu[2] / 2) / (1 + z)
    coefW = -expm1(-y / 2) - sig ** 2 / (1 - nu[0])
    c0 = (xi + zeta) / SQ2; c1 = sig ** 2 / (SQ2 * L2)
    F_main = P * gs * (-expm1(-se2 / 4)) / (2 * (1 + eta))
    F_rem = P * gs * ((c0 * sqrt(2 * se2 / PI) + c1 * se2) / SQ2 + c0 ** 2 + c1 ** 2 * se2) / (2 * (1 - eta))
    kF = e * (1 - nu[0]) / (1 + z)
    twoN_G = b * Psi_lo / (1 + z) - sig ** 2 * gP * Phi_up / (1 - z) - (b * Qb * (nu[0] + tau * nu[2]) * Phi_up / 2 + sbER) / (1 - z)
    twoN_R = kF * (F_main - F_rem) + (coefW * Wlo if coefW > 0 else coefW * (1 + nu[0]) * Wper_up / (1 - z))
    nA = twoAw_G + twoN_G / 2; nB = twoAw_R + twoN_R / 2
    Num = nA + nB; den = dA + dB
    TA = nA / (sig * dA); TB = nB / (sig * dB); lower = Num / (sig * den)
    # Z1R corollary (AM-GM form): T_B >= (1-eta)/(1+eta) (se2/s^2) sqrt(k1 k2)/(sqrt(mu)(c_0 + c_1 A/(2 k1)))
    mu = x * S2 / S0
    k1 = twoAw_R / (2 * sig ** 2 * D0 / (1 + eta))
    rhoF = F_rem / F_main if F_main > 0 else INF
    k2 = kF * (1 - se2 / 8 - rhoF)
    # n_B >= (1/(1+eta)) [2 k1 sigma^2 D0 + k2 P gs se2/16]  (F_main/2 >= P gs se2 (1 - se2/8)/(16(1+eta)); the W-term of n_B is >= 0 and dropped)
    c0_ = (1 + nu[0]) / (1 - z); c1_ = 1 + eps1                 # d_B = (P gs/(2(1-eta))) [c0_ + c1_ * 2 sigma^2 D0/(P gs)]
    if k1 > 0 and k2 > 0 and coefW >= 0 and Wlo >= 0:
        Aamg = (se2 / s2) * sig * sqrt(k1 * k2 / mu)
        # ratio n_B/(sigma d_B) >= ((1-eta)/(1+eta)) * max(Aamg, 2 k1 delta)/(sigma (c0_ + c1_ delta)), delta = 2 sigma^2 D0/(P gs); min over delta at delta* = Aamg/(2 k1)
        TB_cor = (1 - eta) / (1 + eta) * Aamg / (sig * (c0_ + c1_ * Aamg / (2 * k1)))
    else:
        TB_cor = D(0)
    Z1 = (1 - eta) / (1 + eta) * (se2 / s2) / (sqrt(mu) * (1 + sig / (4 * sqrt(mu))))
    out = {'lower_R': lower, 'T_A': TA, 'T_B': TB, 'T_B_corollary': TB_cor, 'Z1_BZ': Z1, 'delta_B': 1 - TB_cor / Z1 if Z1 > 0 else None,
           'k1': k1, 'k2': k2, 'rho_F': rhoF, 'eps1': eps1, 'nu0': nu[0], 'tau': tau}
    if detail:
        out.update({'Num_lo': Num, 'den_up': den, 'nA': nA, 'nB': nB, 'dA': dA, 'dB': dB, 'sbER': sbER, 'corrA_G': corrA_G, 'corrA_R': corrA_R,
                    'F_main': F_main, 'F_rem': F_rem, 'Psi_lo': Psi_lo, 'Phi_up': Phi_up})
    return out


# ------------------------------------------------------------------ direct two-dimensional quadrature (numerical check of Lemma R)
def direct(sig, x, b, L2, s2, q0, ph, nv=6000, nw=900):
    """exact (a), (a') of the periodic window state p = N_V(v) Gamma(w) h(k'v, w)/Z by quadrature on a fine grid."""
    V = sig ** 2 / b; kp = math.sqrt(x) / sig; L = math.sqrt(L2); s = math.sqrt(s2)
    v = np.linspace(-10 * math.sqrt(V), 10 * math.sqrt(V), nv); w = np.linspace(-10 * L, 10 * L, nw)
    dv = v[1] - v[0]; dw = w[1] - w[0]
    VV, WW = np.meshgrid(v, w, indexing='ij'); T = kp * VV; env = np.exp(-WW ** 2 / (2 * s2))
    R = np.zeros_like(T); Rt = np.zeros_like(T)
    for n, (q, p_) in enumerate(zip(q0, ph), start=1):
        R += q * np.cos(n * T - p_); Rt += -n * q * np.sin(n * T - p_)
    h = 1 + env * R; ht = env * Rt; hw = -(WW / s2) * env * R
    NV = np.exp(-VV ** 2 / (2 * V)) / math.sqrt(2 * math.pi * V); G = np.exp(-WW ** 2 / (2 * L2)) / math.sqrt(2 * math.pi * L2)
    p = NV * G * h; Z = p.sum() * dv * dw; p = p / Z
    sv = -VV / V + kp * ht / h; sw = -WW / L2 + hw / h
    mz = ((VV + sig ** 2 * sv) + (WW + sig ** 2 * sw)) / math.sqrt(2); phi = 1 - np.cos(mz)
    E = lambda f: float((p * f).sum() * dv * dw)
    Jvv = E(sv ** 2); Jww = E(sw ** 2)
    den = sig ** 2 * (Jvv + Jww); Aw = sig ** 2 * E(sw ** 2 * (1 - phi / 4)); Nphi = sig ** 2 / 2 * E(sv ** 2 * phi)
    return {'a': 2 * Aw / (sig * den), 'ap': Nphi / (sig * den), 'sum': (2 * Aw + Nphi) / (sig * den), 'Z': Z}


def window_numbers(q0, ph, x):
    t = np.linspace(0, 2 * math.pi, 200001)
    R = sum(q * np.cos((n + 1) * t - p) for n, (q, p) in enumerate(zip(q0, ph)))
    eta = float(np.max(np.abs(R)))
    n = np.arange(1, len(q0) + 1); q0 = np.array(q0)
    return {'S0': float(np.sum(q0 ** 2)), 'S2': float(np.sum(n ** 2 * q0 ** 2)), 'eta': eta * (1 + 1e-6), 'Qb': float(q0.sum()),
            'N1': float(np.sum(n * q0)), 'N2': float(np.sum(n ** 2 * q0))}


def run_direct_checks(nstates=24, seed=11):
    rng = np.random.default_rng(seed); out = []
    for i in range(nstates):
        sig = 10 ** rng.uniform(-1.6, -0.9); tau = 10 ** rng.uniform(0.6, 2.0); x = 10 ** rng.uniform(-1.0, 0.5)
        b = x / tau; L2 = 10 ** rng.uniform(-2.5, -0.7); s2 = L2 / 4 * 10 ** rng.uniform(-1.3, 0)
        if s2 < sig ** 2: s2 = sig ** 2 * 1.5
        if 4 * s2 > L2: L2 = 4 * s2
        nh = rng.integers(1, 5); q0 = list(rng.uniform(0.02, 0.6 / nh, nh)); ph = list(rng.uniform(0, 2 * math.pi, nh))
        wn = window_numbers(q0, ph, x)
        if wn['eta'] >= 0.9 or b >= 1: continue
        kp = math.sqrt(x) / sig; V = sig ** 2 / b
        nv = int(min(40000, max(4000, 30 * len(q0) * kp * 20 * math.sqrt(V) / (2 * math.pi))))
        d = direct(sig, x, b, L2, s2, q0, ph, nv=nv, nw=1200)
        r = lemma_R(sig, x, b, L2, s2, wn['S0'], wn['S2'], wn['eta'], wn['Qb'], wn['N1'], wn['N2'])
        rec = {'sigma': sig, 'tau': tau, 'x': x, 'b': b, 'L2': L2, 's2': s2, 'q0': q0, 'phases': ph, 'eta': wn['eta'], 'Qbar': wn['Qb'],
               'direct_a_plus_ap': d['sum'], 'lemma_R_lower': float(r['lower_R']), 'min_TA_TB': float(min(r['T_A'], r['T_B'])),
               'ok_R': float(r['lower_R']) <= d['sum'], 'ok_Z1R': float(min(r['T_A'], r['T_B'])) <= d['sum'] * (1 + 1e-9), 'ratio': float(r['lower_R']) / d['sum']}
        out.append(rec)
        print('%2d sigma %.3f tau %6.1f x %.3f eta %.3f: direct %.5f  lemma R %.5f (ratio %.4f)  min(T_A,T_B) %.5f  %s' % (
            i, sig, tau, x, wn['eta'], d['sum'], rec['lemma_R_lower'], rec['ratio'], rec['min_TA_TB'], 'ok' if rec['ok_R'] and rec['ok_Z1R'] else 'FAIL'), flush=True)
    return out


def band_points():
    """the hardest band points of BZ (results/bz_asym.jsonl): the state from by_eval.minimise('lo'), the band numbers from profile_rs (logs via mpmath)."""
    sys.path.insert(0, '.')
    import by_eval as Y
    out = []
    for l in open('results/bz_asym.jsonl'):
        d = json.loads(l); e_ = d['sig_exp']; lsig = -e_ * math.log(10); x = 10.0 ** d['log10_x']
        m = max(2, int(round(d['u_m'] / math.sqrt(x)))); p = Y.profile_rs(x, m)
        r = Y.minimise(lsig, x, p['S0'], p['S2'], p['eta'], 'lo')
        sig = D(10) ** (-e_); L2 = exp(D(r['lL2'])); s2 = exp(D(r['ls2'])); b = exp(D(r['lb']))
        uM = (m + p['W']) * math.sqrt(x)                                     # N_1 <= (m + W) Q, N_2 <= (m + W)^2 Q (crude; remainders only)
        N1 = D(uM) / sqrt(x) * D(p['Q']); N2 = (D(uM) / sqrt(x)) ** 2 * D(p['Q'])
        R = lemma_R(sig, x, b, L2, s2, p['S0'], p['S2'], p['eta'], 0.25, N1, N2)
        ell = e_ * math.log(10)
        rec = {'sig_exp': e_, 'gsl_lo_reduced_formula_BY': r['g'] * math.sqrt(ell), 'gsl_lemma_R': float(R['lower_R']) * math.sqrt(ell),
               'gsl_T_A': float(R['T_A']) * math.sqrt(ell), 'gsl_T_B': float(R['T_B']) * math.sqrt(ell), 'gsl_T_B_corollary': float(R['T_B_corollary']) * math.sqrt(ell),
               'gsl_Z1_BZ': float(R['Z1_BZ']) * math.sqrt(ell), 'delta_B': float(R['delta_B']), 'rho_F': float(R['rho_F']),
               'one_minus_k1': float(1 - R['k1']), 'one_minus_k2': float(1 - R['k2']), 'log10_tau': float(log10(R['tau'])),
               'TA_over_TB': float(R['T_A'] / R['T_B'])}
        out.append(rec)
        print('1e-%d: reduced(lo) %.6f  lemma R %.6f  T_B %.6f  Z1R-cor %.6f  Z1(BZ) %.6f  delta_B %.2e  T_A/T_B %.3g' % (
            e_, rec['gsl_lo_reduced_formula_BY'], rec['gsl_lemma_R'], rec['gsl_T_B'], rec['gsl_T_B_corollary'], rec['gsl_Z1_BZ'], rec['delta_B'], rec['TA_over_TB']), flush=True)
    return sorted(out, key=lambda t: t['sig_exp'])


if __name__ == '__main__':
    res = {}
    if 'direct' in sys.argv[1:] or len(sys.argv) == 1:
        res['direct_checks'] = run_direct_checks()
    if 'band' in sys.argv[1:] or len(sys.argv) == 1:
        res['band_points'] = band_points()
    old = {}
    try: old = json.load(open('results/ca_reduction.json'))
    except FileNotFoundError: pass
    old.update(res)
    json.dump(old, open('results/ca_reduction.json', 'w'), indent=1, default=float)
