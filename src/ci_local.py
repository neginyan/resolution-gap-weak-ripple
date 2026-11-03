"""Stage CI 2.1: the local-average picture, measured with LOCAL quantities (beats included).  Python 3.11.

What a narrow focus sees.  The Gaussian v-average of a periodic F is <theta_tau(. - t0) F> = (F * N(0, tau))(t0): F smoothed over a width sqrt(tau)
around the focus position t0 (BL0 (i)).  Lemma R_theta (report CI 2.1) keeps these local averages instead of period averages, so its ripple part has
the SAME algebraic form as Lemma R with
    S0 -> S0_theta = 2 <theta R^2>,   S2 -> S2_theta = 2 <theta R'^2>,   mu_theta = x S2_theta/S0_theta,   P_theta = x S2_theta,
and no tau factor F_tau (the theta -> 1 replacement is never made).  What remains open is a LOCAL Z3': P_theta <= Pi*(mu_theta).

This script measures, on windows built from NON-NEGATIVE raw profiles f (so |a_n| <= 2), blurred with K = x (the sharpest allowed),
scaled to sum q0_n <= Qbar:
  (1) rho = P_theta/Pi*(mu_theta) over all focus positions t0 and tau in a list (the local Z3' test, (a) of the instruction);
  (2) the local reduced lower bound g sqrt(ell) (cc_core's x-free main function with (mu_theta, G = P_theta g_s), b = x/tau, minimised over A exactly and
      over r on a grid; the largest gamma with min_r phi >= 0 by bisection), at the worst t0 -- (b) of the instruction, with local quantities.
Families (all with (P2): a fundamental at least as large as every raw coefficient): beats, wave packets, combs of narrow spikes, random bands.  [numerical; not a proof]
-> results/ci_local.json
"""
import json, math, os
import numpy as np
import cc_core as C

QB = 0.25
OUT = 'results/ci_local.json'
NT = 2 ** 17
T = np.linspace(0, 2 * math.pi, NT, endpoint=False)
K = np.fft.rfftfreq(NT, 1.0 / NT)            # wavenumbers 0..NT/2


def ripple_from_raw(f, x):
    """raw f >= 0 (mean normalised to 1) -> blurred ripple R = sum q0_n cos(n t - phi_n), q0_n = |a_n| e^{-n^2 x/2}, scaled to sum q0 <= Qbar."""
    f = f / f.mean()
    c = np.fft.rfft(f) / NT; c[0] = 0.0
    c = c * np.exp(-K ** 2 * x / 2)
    q = 2 * np.abs(c); qs = q.sum()
    eps = min(1.0, QB / qs) if qs > 0 else 1.0       # f -> 1 + eps (f - 1) keeps f >= 0
    c = c * eps
    R = np.fft.irfft(c * NT, NT); R1 = np.fft.irfft(1j * K * c * NT, NT)
    q = 2 * np.abs(c)
    S0 = float(np.sum(q ** 2)); S2 = float(np.sum(K ** 2 * q ** 2))
    return R, R1, S0, S2, float(q.sum()), eps


def smooth(F, tau):
    return np.fft.irfft(np.fft.rfft(F) * np.exp(-K ** 2 * tau / 2), NT)


def p2_window(p, alpha=0.5):
    """f = 1 + alpha cos t + beta p, p zero-mean with max|p| = 1, beta <= 1 - alpha (f >= 0) and beta max_n a_n(p) <= alpha ((P2): |a_n| <= |a_1|)."""
    p = p - p.mean(); p = p / np.max(np.abs(p))
    c = np.fft.rfft(p) / NT; an = 2 * np.abs(c); an[1] = 0.0
    beta = min(1 - alpha, alpha / an.max())
    f = 1 + alpha * np.cos(T) + beta * p
    cf = np.fft.rfft(f) / NT; af = 2 * np.abs(cf)
    assert f.min() >= -1e-12 and af[2:].max() <= af[1] * (1 + 1e-9)
    return f


def families(rng):
    """windows satisfying (P2) (|a_n| <= |a_1|) and non-negativity: a fundamental plus a high structure that beats or is coherent."""
    out = []
    for m in (60, 200):
        for ph in (0.0, 1.3):
            out.append(('beats m=%d phi=%.1f' % (m, ph), p2_window(np.cos(m * T) + np.cos((m + 1) * T - ph)), m))
    for n, w in ((80, 0.3), (200, 0.1), (200, 0.4)):
        B = np.exp(-(np.angle(np.exp(1j * T))) ** 2 / (2 * w * w))
        out.append(('packet n=%d w=%.1f' % (n, w), p2_window(np.cos(n * T) * B), n))
    for n, d in ((40, 0.02), (120, 0.005)):
        out.append(('comb n=%d width=%.3f' % (n, d), p2_window(np.exp(-(1 - np.cos(n * T)) / (n * n * d * d))), n))
    for n, w in ((4000, 0.01), (12000, 0.003), (12000, 0.03)):
        B = np.exp(-(np.angle(np.exp(1j * (T - math.pi / 2)))) ** 2 / (2 * w * w))      # coherent packet at a ZERO of the fundamental
        out.append(('packet at zero n=%d w=%.3f' % (n, w), p2_window(np.cos(n * T) * B), n))
    for m, kk in ((100, 8), (300, 20)):
        for s in range(2):
            ck = rng.normal(size=kk) + 1j * rng.normal(size=kk)
            g = sum(ck[j] * np.exp(1j * (m + j) * T) for j in range(kk))
            out.append(('random band m=%d K=%d #%d' % (m, kk, s), p2_window(np.real(g)), m))
    return out


def gamma_max(mu, G_over_gs, b, lsig, ell, rgrid):
    """largest gamma with min over r of phi(mu, r, b, G = P g_s(r)) >= 0 (cc_core main function, exact min over A); vectorised over points."""
    sig = math.exp(lsig)
    r = rgrid[None, :]; gs = np.sqrt(r / (2 + r))
    mu = mu[:, None]; P = G_over_gs[:, None]; b = b[:, None] if np.ndim(b) else b
    G = P * gs
    lo = np.zeros(mu.shape[0]); hi = np.full(mu.shape[0], 50.0 / math.sqrt(ell))
    for _ in range(60):
        g = 0.5 * (lo + hi); gg = g[:, None]
        alpha = r * (2 - gg * sig) + G * (1 - gg * sig / 2) / (mu * (2 + r))
        phi = 2 * np.sqrt(alpha * (b / (8 * r) + G / (16 * (2 + r)))) + 3 * sig / 8 + G * sig / (16 * b) - gg * b - gg * G / 2
        ok = np.min(phi, axis=1) >= 0
        lo = np.where(ok, g, lo); hi = np.where(ok, hi, g)
    return lo * math.sqrt(ell)


def main():
    rng = np.random.default_rng(20261103)
    res = {'families': [], 'note': 'numbers with LOCAL quantities (beats included); eta = 0 main function; not a proof'}
    sig_list = [(24, 55.26204223), (48, 110.52408446), (96, 221.04816893)]
    rgrid = np.exp(np.linspace(math.log(1e-9), math.log(0.25), 400))
    worst = {'rho_max': 0.0, 'g_min': {str(e): 1e9 for e, _ in sig_list}}
    for name, f, n in families(rng):
        for mu_target in (2.0, 6.0, 20.0):
            x = mu_target / n ** 2
            if x > 0.003: continue
            R, R1, S0, S2, qsum, eps = ripple_from_raw(f.copy(), x)
            if S0 <= 0: continue
            mu_per = x * S2 / S0; P_per = x * S2
            row = {'family': name, 'x': x, 'mu_period': mu_per, 'P_period': P_per, 'P_over_Pistar_period': P_per / math.exp(C.log_Pistar(mu_per)), 'taus': []}
            for tau in (20.0, 1.0, 0.1, 20.0 / n ** 2, 2.0 / n ** 2):
                if tau < x: continue
                S0t = 2 * smooth(R ** 2, tau); S2t = 2 * smooth(R1 ** 2, tau)
                sel = np.arange(0, NT, NT // 512)
                S0t = np.maximum(S0t[sel], 1e-300); S2t = np.maximum(S2t[sel], 0.0)
                mut = np.maximum(x * S2t / S0t, 1e-12); Pt = np.maximum(x * S2t, 1e-300)
                lPi = np.array([C.log_Pistar(max(m, 1e-12)) for m in mut])
                rho = Pt / np.exp(lPi)
                b = x / tau
                tr = {'tau': tau, 'n2tau': n * n * tau, 'b': b, 'rho_max': float(rho.max()), 'mu_theta_range': [float(mut.min()), float(mut.max())],
                      'P_theta_over_P_period_range': [float(Pt.min() / P_per), float(Pt.max() / P_per)]}
                worst['rho_max'] = max(worst['rho_max'], tr['rho_max'])
                gm = {}
                for e, ell in sig_list:
                    if b > 1: continue
                    gv = gamma_max(mut, Pt, b, -e * math.log(10), ell, rgrid)
                    gm[str(e)] = float(gv.min())
                    worst['g_min'][str(e)] = min(worst['g_min'][str(e)], gm[str(e)])
                tr['g_sqrt_ell_min_over_t0'] = gm
                row['taus'].append(tr)
            res['families'].append(row)
            print(name, 'mu %.2f' % mu_per, ['tau %.3g rho %.3g g %s' % (t['tau'], t['rho_max'], {k: round(v, 3) for k, v in t['g_sqrt_ell_min_over_t0'].items()}) for t in row['taus']], flush=True)
    res['worst'] = worst
    os.makedirs('results', exist_ok=True)
    json.dump(res, open(OUT, 'w'), indent=1)
    print(json.dumps(worst, indent=1))


if __name__ == '__main__':
    main()
