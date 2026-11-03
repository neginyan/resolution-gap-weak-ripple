"""Stage BY: are band profiles really below the comb?  The two branches with an UPPER bound of g (side='hi'), the sandwich narrowed by the
actual size of the ripple, the Rudin-Shapiro band (nonnegativity proved), sigma down to 1e-300 in log arithmetic.  NUMERICAL (the sandwich and
the reduction are derivations; see bx_map.py).

Reduction (bx_map.py, separated regime, small angles: every length involved is below 1e-5, so phi = (v' + w)^2/4 to relative 1e-10):
  sigma^2 J_ww = sigma^2 (1/L^2 + S0 g_s se2 c/(2 s^4)),  sigma^2 J_vv = b + x S2 g_s c/2,  A_w = sigma^2 J_ww,
  N_phi = (1/2)[ b (3V + L^2)/4 + x S2 g_s c (V + se2)/8 ],  V = sigma^2/b,  se2 = 1/(1/L^2 + 2/s^2),  g_s = sqrt(se2/L^2),
  (a) = 2 A_w/(sigma den),  (a') = N_phi/(sigma den),  den = sigma^2 (J_vv + J_ww),
  c in [1/(1 + eta), 1/(1 - eta)] with eta = sup_{t,w} |h - 1| (Lemma: <F^2/h> lies between <F^2>/(1 + eta) and <F^2>/(1 - eta)); 'hi' takes
  1/(1 - eta) in the numerators and 1/(1 + eta) in den, 'lo' the reverse, 'mid' c = 1.  bx_map used eta = Qbar (the crude sandwich).
Everything is computed with logarithms (sigma^2 underflows for sigma < 1e-154).
Rudin-Shapiro band (PROOF of nonnegativity): f = 1 + c cos t + c sum_{k<W} e_k cos((m + k) t), e_k the Rudin-Shapiro signs, W = 2^j.  The pair
(P_j, Q_j) satisfies |P|^2 + |Q|^2 = 2W on |z| = 1, so |sum e_k e^{i(m+k)t}| <= sqrt(2W) and f >= 1 - c(1 + sqrt(2W)) >= 0 for c <= 1/(1 + sqrt(2W)).
(P2): |a_n| = c = |a_1|.  The blurred ripple: eta <= c e^{-K/2} + c (2 + sqrt2)[d_0 + int_0^inf e^{-(m + t^2)^2 K/2} dt] (Abel summation with the
Rudin-Shapiro partial sums |sum_{k<n} e_k z^k| <= (2 + sqrt2) sqrt n (Shapiro), d_k = e^{-(m+k)^2 K/2} decreasing, sum_{k>=1} d_k (sqrt(k+1) - sqrt k)
<= d_0 sqrt(pi/(m K))/2, since (m + k)^2 >= m^2 + 2mk).
-> results/by_map.jsonl (one line per (sigma, x)), results/by_compare.json, results/by_direct.json.  Python 3.11."""
import json, math, os, sys
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import quad

LN2 = math.log(2.0)
C_RS = 2 + math.sqrt(2)


def lae(*a):
    """log(sum exp(a_i)) for finite logs."""
    m = max(a); return m + math.log(sum(math.exp(t - m) for t in a))


def gsum(p, m, K, M):
    """sum_{n=m}^{M} n^p e^{-n^2 K}: explicit below 2e6 terms, else the integral with the half end-term corrections (relative error < 1e-6 here)."""
    if M - m < 2_000_000:
        n = np.arange(m, M + 1, dtype=float); return float(np.sum(n ** p * np.exp(-n * n * K)))
    from scipy.special import gammaincc, gamma
    a = (p + 1) / 2
    diff = gammaincc(a, m * m * K) - gammaincc(a, M * M * K)
    val = 0.5 * math.exp(-a * math.log(K) + math.lgamma(a) + math.log(diff)) if diff > 0 else 0.0     # in logarithms (K^{-a} overflows for K < 1e-200)
    return val + 0.5 * (m ** p * math.exp(-m * m * K) + M ** p * math.exp(-M * M * K))


def rs_signs(W):
    n = np.arange(W, dtype=np.int64); cnt = np.zeros(W, dtype=np.int64); v = n.copy()
    while np.any(v):
        cnt += ((v & 3) == 3); v >>= 1
    return np.where(cnt % 2 == 0, 1.0, -1.0)


def profile_rs(x, m, K=None, Qbar=0.25):
    """the Rudin-Shapiro band starting at m (the band has 20 decay lengths, M = m + 20/(m K)); amplitude c = min(nonnegativity, Qbar).
    Returns S0, S2, Q, eta and the two amplitude limits."""
    K = x if K is None else K
    M = m + max(1, int(math.ceil(20 / (m * K))))           # 20 decay lengths: (m + k)^2 >= m^2 + 2mk, so the cut-off is below e^{-20} relative
    j = max(1, math.ceil(math.log2(M - m + 1))); W = 2 ** j
    c_nn = 1 / (1 + math.sqrt(2 * W))                             # f >= 1 - c(1 + sqrt(2W)) >= 0 (proved in the docstring)
    e1 = math.exp(-K / 2)
    band_q = gsum(0, m, K / 2, m + W - 1)
    c_Q = Qbar / (e1 + band_q)
    c = min(c_nn, c_Q)
    S0 = c * c * (e1 * e1 + gsum(0, m, K, m + W - 1)); S2 = c * c * (e1 * e1 + gsum(2, m, K, m + W - 1)); Q = c * (e1 + band_q)
    # sum_{k>=1} d_k/(2 sqrt k) <= int_0^inf d(k)/(2 sqrt k) dk <= d_0 int_0^inf e^{-m k K}/(2 sqrt k) dk = d_0 sqrt(pi/(m K))/2  ((m+k)^2 >= m^2 + 2mk)
    integ = math.exp(-m * m * K / 2) * 0.5 * math.sqrt(math.pi / (m * K))
    eta = c * e1 + c * C_RS * (math.exp(-m * m * K / 2) + integ)
    return {'m': m, 'u_m': m * math.sqrt(x), 'W': W, 'c': c, 'c_nonneg': c_nn, 'c_Qbar': c_Q, 'binding': 'Qbar' if c_Q < c_nn else 'nonneg',
            'S0': S0, 'S2': S2, 'Q': Q, 'eta': eta}


def branches_log(lsig, x, S0, S2, eta, lL2, ls2, lb, side='hi'):
    """(a), (a') from logarithms of sigma, L^2, s^2, b (see the module docstring)."""
    cn, cd = {'mid': (1.0, 1.0), 'lo': (1 / (1 + eta), 1 / (1 - eta)), 'hi': (1 / (1 - eta), 1 / (1 + eta))}[side]
    lse2 = -lae(-lL2, LN2 - ls2); lgs = (lse2 - lL2) / 2
    lV = 2 * lsig - lb
    Jww_n = lae(2 * lsig - lL2, 2 * lsig + math.log(S0 * cn / 2) + lgs + lse2 - 2 * ls2)
    Jww_d = lae(2 * lsig - lL2, 2 * lsig + math.log(S0 * cd / 2) + lgs + lse2 - 2 * ls2)
    Jvv_d = lae(lb, math.log(x * S2 * cd / 2) + lgs)
    lden = lae(Jvv_d, Jww_d)
    lN = math.log(0.5) + lae(lb + math.log(0.25) + lae(math.log(3) + lV, lL2), math.log(x * S2 * cn / 8) + lgs + lae(lV, lse2))
    a = 2 * math.exp(Jww_n - lsig - lden); ap = math.exp(lN - lsig - lden)
    return a, ap


def minimise(lsig, x, S0, S2, eta, side='hi'):
    """min over L^2 <= 1/sigma ((P0)), sigma^2 <= s^2 <= L^2/4 ((P4)), b <= x/80 (separated regime), b >= sigma^2 1e-6."""
    lL2max = -lsig; lbmax = math.log(x / 80); lbmin = 2 * lsig - 6 * math.log(10)
    if lbmax <= lbmin: return None
    def unpack(z):
        lL2 = min(z[0], lL2max); ls2 = min(max(z[1], 2 * lsig), lL2 - math.log(4)); lb = min(max(z[2], lbmin), lbmax)
        return lL2, ls2, lb
    def f(z):
        a, ap = branches_log(lsig, x, S0, S2, eta, *unpack(z), side); return a + ap
    best = None
    for dL in (0.5, 0.0, -2.0):
        for ds in (0.5, -1.0):
            for db in (3.0, 6.0):
                z0 = [0.5 * lsig + dL, lsig + ds, min(lsig + db, lbmax)]
                r = minimize(f, z0, method='Nelder-Mead', options={'xatol': 1e-9, 'fatol': 1e-14, 'maxiter': 2500})
                if best is None or r.fun < best.fun: best = r
    lL2, ls2, lb = unpack(best.x)
    a, ap = branches_log(lsig, x, S0, S2, eta, lL2, ls2, lb, side)
    return {'g': a + ap, 'a': a, 'ap': ap, 'lL2': lL2, 'ls2': ls2, 'lb': lb, 'L2_over_sqrt_sigma': math.exp(lL2 - lsig / 2), 's2_over_sigma': math.exp(ls2 - lsig),
            'b_over_sigma': math.exp(lb - lsig), 'log10_tau': (math.log(x) - lb) / math.log(10)}


def scan(sig_exp, xs_exp, us=None, fn='results/by_map.jsonl', umax=7.0):
    """for sigma = 10^-sig_exp and x = 10^e (e in xs_exp): the best Rudin-Shapiro band over the start u_m, judged by 'hi'; lo and mid at the same profile."""
    have = set()
    if os.path.exists(fn):
        for l in open(fn):
            if l.strip(): d = json.loads(l); have.add((d['sig_exp'], d['x_exp']))
    lsig = -sig_exp * math.log(10); ell = sig_exp * math.log(10)
    us = us if us is not None else np.geomspace(1.0, umax, int(round(28 * math.log(umax) / math.log(7.0))))
    for e in xs_exp:
        if (sig_exp, e) in have: continue
        x = 10.0 ** e; best = None
        for u in us:
            m = max(2, int(round(u / math.sqrt(x)))); p = profile_rs(x, m)
            r = minimise(lsig, x, p['S0'], p['S2'], p['eta'], 'hi')
            if r is None: continue
            if best is None or r['g'] < best['hi']['g']:
                lo = minimise(lsig, x, p['S0'], p['S2'], p['eta'], 'lo'); mid = minimise(lsig, x, p['S0'], p['S2'], p['eta'], 'mid')
                best = {'profile': p, 'hi': r, 'lo': lo, 'mid': mid}
        if best is None: continue
        d = {'sig_exp': sig_exp, 'x_exp': e, 'ell': ell, 'judged_by': 'hi', 'g_hi': best['hi']['g'], 'g_mid': best['mid']['g'], 'g_lo': best['lo']['g'],
             'gsl_hi': best['hi']['g'] * math.sqrt(ell), 'gsl_mid': best['mid']['g'] * math.sqrt(ell), 'gsl_lo': best['lo']['g'] * math.sqrt(ell),
             'theorem9': 1 - 1.7 * math.log(ell) / ell, **best}
        open(fn, 'a').write(json.dumps(d, default=float) + '\n')
        print('sigma 1e-%d x 1e%.1f: g sqrt(ell) hi %.4f mid %.4f lo %.4f  u_m %.2f eta %.2e c-limit %s' % (
            sig_exp, e, d['gsl_hi'], d['gsl_mid'], d['gsl_lo'], best['profile']['u_m'], best['profile']['eta'], best['profile']['binding']), flush=True)


if __name__ == '__main__':
    sig_exp = int(sys.argv[1]); lo_e = float(sys.argv[2]); hi_e = float(sys.argv[3]); step = float(sys.argv[4])
    umax = float(sys.argv[5]) if len(sys.argv) > 5 else 7.0
    scan(sig_exp, list(np.arange(hi_e, lo_e - 1e-9, -step)), umax=umax)
