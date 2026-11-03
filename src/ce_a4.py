"""Stage CE 2.1: LEMMA A4* (the a_4 remainder of v1.5 BL3 replaced by a RELATIVE bound) and Lemma R with it (Lemma R*).  Python 3.11.

Where a_4 comes from (v1.5 App. D, proof of BL3): the F-term of the exact identity (BLid) is sigma^2 (G[F] - Re A[F e^{i psi}]), F = h'^2/h >= 0,
G[f] = <theta_tau f> (coefficients e_n = e^{-n^2 tau/2}), A[f] = e <Theta f> with Theta the shifted kernel (coefficients e_n e^{-n beta} at frequency n,
n in Z), e = e^{-y/2}.  v1.5 writes  G[F] - Re A[F e^{i psi}] = (1 - e)<theta F> + e<theta F (1 - cos psi)> + e Re<(theta - Theta) F e^{i psi}>  and bounds the
last term by sup F * ||theta - Theta||_A <= k'^2 N_1^2/(1 - Qbar) * beta nu_1, i.e. a_4 = sigma tau sqrt(x) N_1^2 nu_1/(sqrt2 (1 - Qbar)).  This is an ABSOLUTE
term ~ e^{-tau/2}/k' and it makes the background part T_A negative for 20 <= tau <~ 0.48 ell (stage CD).

LEMMA A4* [proof].  For every w, with g = F e^{i psi}, psi = w/sqrt2 + chi(t, w), ghat_n its Fourier coefficients (|ghat_n| <= <|g|> = <F>, since F >= 0):
    Re<(theta - Theta) g> = sum_{n>=1} e_n Re[ (ghat_n + ghat_{-n})(1 - cosh n beta) + (ghat_n - ghat_{-n}) sinh n beta ].
  (i) |(ghat_n + ghat_{-n})(1 - cosh n beta)| <= 2<F> (n beta)^2 e^{n beta}/2      (cosh u - 1 <= (u^2/2) e^u)
      -> sum <= (beta^2 nu_2/2) <F>.
  (ii) split e^{i chi} = 1 + (e^{i chi} - 1).  The '1' part gives e^{i w/sqrt2}(Fhat_n - Fhat_{-n}) sinh n beta = 2i e^{i w/sqrt2} Im(Fhat_n) sinh n beta, whose real part
      is -2 sin(w/sqrt2) Im(Fhat_n) sinh n beta.  F(t, w) is even in w (h - 1 = e^{-w^2/2s^2} R(t)), so its integral against the centred Gamma VANISHES EXACTLY
      (the same parity argument v1.5 uses for a_10).  The rest has Fourier coefficients bounded by <F |e^{i chi} - 1|> <= sup_t|chi| <F>, and sinh u <= u e^u:
      -> sum <= beta nu_1 sup_t|chi(., w)| <F>_w.
  Integrating against Gamma (sigma^2 int Gamma <F> = W_per, per-w <F>_w <= x S2 e^{-w^2/s^2}/(2 sigma^2 (1 - eta)), sup_t|chi| <= chibar(w) of BL3):
    e sigma^2 |int Gamma Re<(theta - Theta) F e^{i psi}>| <= a_4* := e W_per,up [ beta^2 nu_2/2 + beta nu_1 ((xi + zeta)/sqrt2 + sigma^2 sqrt(se2/pi)/L^2) ],
  with beta^2 = tau y (beta = alpha k' V, y = alpha^2 V, tau = k'^2 V).  So a_4 is replaced by a_4*, which is proportional to the ripple (P) and goes into the
  RIPPLE part n_B.  Relative to the W-coefficient (1 - e) >= (y/2)/(1 + y/2) of BL3 its first part is ~ tau nu_2 and its second ~ nu_1 sqrt(tau) chibar/sqrt(y)
  = nu_1 chibar sqrt(2 b tau)/sigma ~ 2 nu_1 U sqrt(x): both are tau-factors (CB), uniformly small for tau >= 20.  Every other a_i is unchanged.
Parts: (1) numerical check of Lemma A4* on random periodic windows (exact Fourier sums + Gauss-Hermite in w); (2) Lemma R* = Lemma R with a_4 -> a_4*,
T_A and T_B checked SEPARATELY on the band windows forced to tau = x/b in [20, 1e6] and on random states.  -> results/ce_a4.json.
"""
import json, math, sys
import numpy as np
import ca_reduction as CR
from ca_dec import D, exp, sqrt, PI, E, log10

SQ2 = sqrt(2)


# ---------------- (1) numerical check of Lemma A4* ----------------
def check_A4(ntrial=300, seed=3):
    rng = np.random.default_rng(seed); worst = 0.0; viol = 0; rows = 0
    N = 512; t = np.arange(N) * 2 * np.pi / N
    gh_x, gh_w = np.polynomial.hermite_e.hermegauss(80)                     # int N(0,1) f = sum w_i f(x_i)/sqrt(2 pi)
    for it in range(ntrial):
        nh = int(rng.integers(1, 30)); n0 = int(rng.integers(1, 40)); ns = np.arange(n0, n0 + nh)
        a = rng.uniform(-1, 1, nh) * (rng.random(nh) < 0.8); ph = rng.uniform(0, 2 * np.pi, nh) if it % 2 else np.zeros(nh)
        R = (a[:, None] * np.cos(ns[:, None] * t + ph[:, None])).sum(0); dR = (-a[:, None] * ns[:, None] * np.sin(ns[:, None] * t + ph[:, None])).sum(0)
        sc = rng.uniform(0.05, 0.9) / max(np.max(np.abs(R)), 1e-12); R *= sc; dR *= sc
        L = 10 ** rng.uniform(-0.5, 0.7); s = L / 2 * 10 ** rng.uniform(-1, 0)
        tau = 10 ** rng.uniform(-0.3, 1.5); beta = 10 ** rng.uniform(-3, 0.3)
        chiamp = 10 ** rng.uniform(-4, 0); chi_t = np.cos(int(rng.integers(1, 9)) * t + rng.uniform(0, 6.3)) * rng.uniform(0, 1)
        nmax = 60; nn = np.arange(1, nmax + 1); en = np.exp(-nn * nn * tau / 2)
        nu1 = 2 * np.sum(nn * en * np.exp(nn * beta)); nu2 = 2 * np.sum(nn * nn * en * np.exp(nn * beta))
        lhs = 0.0; rhs = 0.0
        for xi, wi in zip(gh_x, gh_w):
            w = L * xi; env = math.exp(-w * w / (2 * s * s))
            F = env ** 2 * dR ** 2 / (1 + env * R)
            chi = chiamp * (1 + abs(w)) * chi_t; supchi = np.max(np.abs(chi))
            g = F * np.exp(1j * (w / math.sqrt(2) + chi))
            gh = np.fft.fft(g) / N                                          # ghat_k for e^{-ikt}: g = sum ghat_k e^{ikt} with k = 0..N-1 (negative wrap)
            gp = gh[nn]; gm = gh[-nn]
            val = np.sum(en * ((gp + gm) * (1 - np.cosh(nn * beta)) + (gp - gm) * np.sinh(nn * beta))).real
            mF = float(np.mean(F))
            lhs += wi * val; rhs += wi * (beta ** 2 * nu2 / 2 * mF + beta * nu1 * supchi * mF)
        lhs /= math.sqrt(2 * math.pi); rhs /= math.sqrt(2 * math.pi)
        rows += 1; viol += abs(lhs) > rhs * (1 + 1e-9) + 1e-300; worst = max(worst, abs(lhs) / rhs if rhs > 0 else 0)
    return {'trials': rows, 'violations': int(viol), 'max_ratio_lhs_over_bound': worst}


# ---------------- (2) Lemma R* ----------------
def a4_terms(sig, x, b, L2, s2, S0, S2, eta, Qb, N1):
    sig, x, b, L2, s2, S0, S2, eta, Qb, N1 = (D(v) for v in (sig, x, b, L2, s2, S0, S2, eta, Qb, N1))
    tau = x / b; V = sig ** 2 / b; al = (1 - b) / SQ2; y = al ** 2 * V; e = exp(-y / 2); kp = sqrt(x) / sig; beta = al * kp * V
    nu, om2 = CR.series(tau, beta); z = Qb * nu[0] / 2
    se2 = 1 / (1 / L2 + 2 / s2); gs = sqrt(se2 / L2); s = sqrt(s2); L = sqrt(L2); P = x * S2
    Wper_up = P * gs / (2 * (1 - eta))
    xi = sig * sqrt(x) * N1 / (1 - Qb); zeta = sig ** 2 * eta * exp(D('-0.5')) / (s * (1 - eta))
    a4_old = sig * tau * sqrt(x) * N1 ** 2 * nu[1] / (SQ2 * (1 - Qb))
    a4_new = e * Wper_up * (beta ** 2 * nu[2] / 2 + beta * nu[1] * ((xi + zeta) / SQ2 + sig ** 2 * sqrt(se2 / PI) / L2))
    return a4_old, a4_new, z, {'beta2_over_tau_y': beta ** 2 / (tau * y), 'rel_to_W_coef': a4_new / Wper_up / max(-((-y / 2).exp() - 1), D('1e-9999'))}


def _lemma_R_no_a4():
    """ca_reduction.lemma_R with the a_4 entry of the remainder list set to 0 (same source otherwise), so that nothing has to be added back."""
    import inspect
    src = inspect.getsource(CR.lemma_R)
    old = "sig * tau * sqrt(x) * N1 ** 2 * nu[1] / (SQ2 * (1 - Qb)),"
    assert src.count(old) == 1
    ns = dict(CR.__dict__); exec(src.replace(old, "D(0),").replace("def lemma_R(", "def lemma_R_no_a4("), ns)
    return ns['lemma_R_no_a4']


LEMMA_R_NO_A4 = _lemma_R_no_a4()


def lemma_R_star(sig, x, b, L2, s2, S0, S2, eta, Qb, N1, N2):
    """Lemma R* = Lemma R with a_4 removed from the background part and a_4* subtracted from the ripple part (n_B)."""
    R = LEMMA_R_NO_A4(sig, x, b, L2, s2, S0, S2, eta, Qb, N1, N2, detail=True)
    Rold = CR.lemma_R(sig, x, b, L2, s2, S0, S2, eta, Qb, N1, N2, detail=True)
    a4o, a4n, z, info = a4_terms(sig, x, b, L2, s2, S0, S2, eta, Qb, N1)
    nA = R['nA']; nB = R['nB'] - a4n / (2 * (1 - z))
    sigD = D(sig)
    return {'T_A': nA / (sigD * R['dA']), 'T_B': nB / (sigD * R['dB']), 'lower': (nA + nB) / (sigD * (R['dA'] + R['dB'])),
            'T_A_old': Rold['T_A'], 'T_B_old': Rold['T_B'], 'a4_old': a4o, 'a4_star': a4n, 'a4_star_over_nB': a4n / (2 * (1 - z)) / nB if nB != 0 else None,
            'num': nA + nB, 'den': R['dA'] + R['dB'], 'den_A': R['dA'], 'den_B': R['dB'], **info}


def band_state(e):
    import by_eval as Y
    ref = {json.loads(l)['sig_exp']: json.loads(l) for l in open('results/bz_asym.jsonl')}[e]
    x = 10 ** ref['log10_x']; m = max(2, int(round(ref['u_m'] / math.sqrt(x)))); p = Y.profile_rs(x, m)
    lsig = -e * math.log(10); r = Y.minimise(lsig, x, p['S0'], p['S2'], p['eta'], 'lo')
    uM = (m + p['W']) * math.sqrt(x); N1 = D(uM) / sqrt(D(x)) * D(p['Q']); N2 = (D(uM) / sqrt(D(x))) ** 2 * D(p['Q'])
    return x, p, r, N1, N2


def forced_tau_table():
    out = []
    for e in (24, 48, 96, 150, 300):
        x, p, r, N1, N2 = band_state(e); ell = e * math.log(10); sig = D(10) ** (-e)
        for tau in (20, 40, 80, 330, 1e3, 1e6, None):
            b = exp(D(r['lb'])) if tau is None else D(x) / D(tau)
            for lL2, ls2, lab in ((r['lL2'], r['ls2'], 'optimum L2,s2'), (r['lL2'] + 3, r['ls2'] + 3, 'L2,s2 x e^3'), (r['lL2'] - 3, r['ls2'] - 3, 'L2,s2 x e^-3')):
                if lL2 > ell: continue
                Rs = lemma_R_star(sig, x, b, exp(D(lL2)), exp(D(ls2)), p['S0'], p['S2'], p['eta'], 0.25, N1, N2)
                f = lambda v: float(v) * math.sqrt(ell)
                out.append({'sig_exp': e, 'tau': tau if tau else float(D(x) / b), 'state': lab, 'TA_sqrt_ell_old': f(Rs['T_A_old']), 'TA_sqrt_ell_new': f(Rs['T_A']),
                            'TB_sqrt_ell_old': f(Rs['T_B_old']), 'TB_sqrt_ell_new': f(Rs['T_B']), 'a4star_over_nB': float(Rs['a4_star_over_nB'])})
                print(out[-1], flush=True)
    return out


def random_states(n=3000, seed=5):
    # b > sigma^2 only: Lemma R as coded needs y < 1 (b <~ sigma^2 is the saturation edge of 2.2)
    """random states of every size (edges included): b in [1.1 sigma^2, x/20], L^2 in [4 sigma^2, 1/sigma], s^2 in [sigma^2, L^2/4]; T_A, T_B separately."""
    rng = np.random.default_rng(seed); out = {}
    for e in (48, 300):
        x, p, r, N1, N2 = band_state(e); ell = e * math.log(10); sig = D(10) ** (-e); lsig = -ell
        mins = {'TA_new': 1e300, 'TA_old': 1e300, 'TB_new': 1e300, 'TB_old': 1e300, 'neg_TA_old': 0, 'neg_TA_new': 0, 'max_a4star_over_nB': 0.0}
        for _ in range(n):
            lb = rng.uniform(2 * lsig + 0.1, math.log(x / 20)); lL2 = rng.uniform(math.log(4) + 2 * lsig, -lsig); ls2 = rng.uniform(2 * lsig, lL2 - math.log(4))
            Rs = lemma_R_star(sig, x, exp(D(lb)), exp(D(lL2)), exp(D(ls2)), p['S0'], p['S2'], p['eta'], 0.25, N1, N2)
            ta, tb, tao, tbo = (float(v) * math.sqrt(ell) for v in (Rs['T_A'], Rs['T_B'], Rs['T_A_old'], Rs['T_B_old']))
            mins['TA_new'] = min(mins['TA_new'], ta); mins['TA_old'] = min(mins['TA_old'], tao); mins['TB_new'] = min(mins['TB_new'], tb); mins['TB_old'] = min(mins['TB_old'], tbo)
            mins['neg_TA_old'] += tao < 0; mins['neg_TA_new'] += ta < 0
            if ta < 0: mins['max_log10_b_over_sigma2_where_TA_new_negative'] = max(mins.get('max_log10_b_over_sigma2_where_TA_new_negative', -1e9), (lb - 2 * lsig) / math.log(10))
            if tb < 0:
                mins['neg_TB_new'] = mins.get('neg_TB_new', 0) + 1
                mins['max_log10_dB_over_dA_where_TB_negative'] = max(mins.get('max_log10_dB_over_dA_where_TB_negative', -1e9), float((Rs['den_B'] / Rs['den_A']).log10()) if Rs['den_B'] > 0 else -1e9)
            if Rs['a4_star_over_nB'] is not None: mins['max_a4star_over_nB'] = max(mins['max_a4star_over_nB'], abs(float(Rs['a4_star_over_nB'])))
        mins['n'] = n; out['1e-%d' % e] = mins
        print(e, mins, flush=True)
    return out


if __name__ == '__main__':
    res = {'check_A4': check_A4()}
    print(res['check_A4'], flush=True)
    res['forced_tau'] = forced_tau_table()
    res['random_states'] = random_states()
    json.dump(res, open('results/ce_a4.json', 'w'), indent=1)
