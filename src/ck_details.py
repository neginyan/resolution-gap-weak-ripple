"""Independent checks of the line-by-line proofs of App. D (Lemmas S, S', A4*, Z3', BE).  Python 3.11, numpy only.
sympy and mpmath are not available in the lab (the package index is not reachable); exact identities are checked with fractions.Fraction on
random rationals, constants with decimal at 60 digits, and Fourier statements with exact trigonometric sums on fine grids (errors ~1e-15).
-> results/ck_details.json.  Every entry has 'ok'."""
import json, math, os, random
from fractions import Fraction as Fr
from decimal import Decimal as D, getcontext
import numpy as np

getcontext().prec = 60
OUT = 'results/ck_details.json'
rnd = random.Random(20261103)


def rq():
    return Fr(rnd.randint(-999, 999), rnd.randint(1, 999))


def trig(nmax, rng, lo=1, hi=None, cos_only=False):
    hi = hi or nmax
    n = np.arange(lo, hi + 1); a = rng.normal(size=len(n)); ph = rng.uniform(0, 2 * np.pi, len(n))
    return n, a, ph


def ev(n, a, ph, t, d=0):
    """sum a cos(n t - ph) and its derivatives."""
    arg = np.outer(t, n) - ph
    if d == 0: return np.cos(arg) @ a
    return (-np.sin(arg) * n) @ a


def check_S(rng, ntrial=300):
    res = {}
    # D.2 Step 2: 1/(1+r) = 1 - r + r^2/(1+r), exact on random rationals with |r| < 1
    ok = True
    for _ in range(2000):
        r = Fr(rnd.randint(-998, 998), 999); ok &= (Fr(1) / (1 + r) == 1 - r + r * r / (1 + r))
    res['step2_identity_exact'] = ok
    # Step 3: the six-product expansion of Y^2 r, exact polynomial identity on random rationals
    ok = True
    for _ in range(2000):
        Yl, Yh, rl, rh = rq(), rq(), rq(), rq()
        lhs = (Yl + Yh) ** 2 * (rl + rh)
        rhs = Yl ** 2 * rl + Yl ** 2 * rh + 2 * Yl * Yh * rl + 2 * Yl * Yh * rh + Yh ** 2 * rl + Yh ** 2 * rh
        ok &= lhs == rhs
    res['step3_expansion_exact'] = ok
    # Step 3 (b), (c): the two products vanish with the gap 2M < m (exact trigonometric sums: period average on a grid finer than the frequencies)
    worst_gap = 0.0; worst_nogap = 0.0
    for _ in range(100):
        M = int(rng.integers(1, 6)); m = 2 * M + int(rng.integers(1, 5)); top = m + int(rng.integers(0, 8))
        t = np.linspace(0, 2 * np.pi, 4 * (top + 2) * 8, endpoint=False)
        nl, al, pl = trig(M, rng); nh, ah, phh = trig(top, rng, lo=m)
        for d in (0, 1):
            yl = ev(nl, al, pl, t, d); yh = ev(nh, ah, phh, t, d); Rl = ev(nl, al, pl, t); Rh = ev(nh, ah, phh, t)
            worst_gap = max(worst_gap, abs(np.mean(yl ** 2 * Rh)), abs(np.mean(yl * yh * Rl)))
        m2 = M + 1                                              # no gap: the products need not vanish (control)
        nh2, ah2, ph2 = trig(m2 + 3, rng, lo=m2); yh2 = ev(nh2, ah2, ph2, t); Rh2 = ev(nh2, ah2, ph2, t)
        worst_nogap = max(worst_nogap, abs(np.mean(ev(nl, al, pl, t) ** 2 * Rh2)))
    res['step3_gap_products_max'] = worst_gap; res['step3_control_without_gap_max'] = worst_nogap
    # Step 3 (e) and Step 4: the bounds with kappa, Minkowski, and Lemma S / S' themselves, with a w-measure (random weights) and 0 <= psi <= 1
    viol = 0; tests = 0
    for _ in range(ntrial):
        gap = rng.uniform() < 0.5
        M = int(rng.integers(1, 5)); m = (2 * M + 1 if gap else M + 1) + int(rng.integers(0, 4)); top = m + int(rng.integers(0, 10))
        t = np.linspace(0, 2 * np.pi, 16 * (top + 2), endpoint=False)
        nl, al, pl = trig(M, rng); nh, ah, phh = trig(top, rng, lo=m)
        al *= rng.uniform(0.01, 1); ah *= rng.uniform(0.01, 1)
        Rl = ev(nl, al, pl, t); Rh = ev(nh, ah, phh, t); R = Rl + Rh
        sc = rng.uniform(0.05, 0.95) / np.max(np.abs(R)); Rl *= sc; Rh *= sc; R *= sc; al = al * sc; ah = ah * sc
        nw = 7; wts = rng.uniform(0.1, 1, nw); psi = rng.uniform(0, 1, nw); phi = rng.normal(size=nw)
        for d in (0, 1):
            yl = ev(nl, al, pl, t, d); yh = ev(nh, ah, phh, t, d); y = yl + yh
            Y2 = np.array([wts[k] * phi[k] ** 2 * np.mean(y ** 2) for k in range(nw)]).sum()
            Ylo2 = sum(wts[k] * phi[k] ** 2 for k in range(nw)) * np.mean(yl ** 2); Yhi2 = Y2 - Ylo2
            mid = sum(wts[k] * phi[k] ** 2 * np.mean(y ** 2 / (1 + psi[k] * R)) for k in range(nw))
            flo = Ylo2 / Y2; fhi = 1 - flo
            el, eh, eta = np.max(np.abs(Rl)), np.max(np.abs(Rh)), np.max(np.abs(R))
            kap = min(el, np.max(np.abs(yh)) * math.sqrt(np.mean(Rl ** 2)) / math.sqrt(np.mean(yh ** 2)))
            E = el * flo + 2 * eh * math.sqrt(flo * fhi) + (eh + kap) * fhi
            if not gap: E += eh * flo + 2 * el * math.sqrt(flo * fhi)
            B = (el * math.sqrt(flo) + kap * math.sqrt(fhi) + eh) ** 2
            lo = Y2 * max(1 / (1 + eta), 1 - E); up = Y2 * min(1 / (1 - eta), 1 + E + B / (1 - eta))
            tests += 1; viol += not (lo <= mid * (1 + 1e-12) and mid <= up * (1 + 1e-12))
            # Step 3(e) second bound and Step 4 middle term on one w slice
            lhs = abs(np.mean(yh ** 2 * Rl)); rhs = np.max(np.abs(yh)) * math.sqrt(np.mean(yh ** 2)) * math.sqrt(np.mean(Rl ** 2))
            viol += lhs > rhs * (1 + 1e-12)
            viol += math.sqrt(np.mean((yh * Rl) ** 2)) > kap * math.sqrt(np.mean(yh ** 2)) * (1 + 1e-12)
    res['lemma_S_Sprime_random'] = {'tests': tests, 'violations': int(viol)}
    res['ok'] = res['step2_identity_exact'] and res['step3_expansion_exact'] and worst_gap < 1e-12 and worst_nogap > 1e-3 and viol == 0
    return res


def check_A4(rng, ntrial=200):
    res = {}
    # D.3: <(theta_tau - Theta) g>_t = sum_n e_n [(g_n + g_-n)(1 - cosh n beta) + (g_n - g_-n) sinh n beta] for complex g
    worst = 0.0
    for _ in range(ntrial):
        K = int(rng.integers(1, 12)); tau = rng.uniform(0.5, 30); beta = rng.uniform(0, 0.3)
        ns = np.arange(-K, K + 1); gh = rng.normal(size=len(ns)) + 1j * rng.normal(size=len(ns))
        t = np.linspace(0, 2 * np.pi, 64 * (K + 1), endpoint=False)
        g = (np.exp(1j * np.outer(t, ns)) @ gh)
        N = 60; nn = np.arange(-N, N + 1); e = np.exp(-nn ** 2 * tau / 2)
        theta = (np.exp(1j * np.outer(t, nn)) @ e).real; Theta = np.exp(1j * np.outer(t, nn)) @ (e * np.exp(nn * beta))
        direct = np.mean((theta - Theta) * g)
        gn = dict(zip(ns, gh)); s = 0
        for n in range(1, K + 1):
            en = math.exp(-n * n * tau / 2)
            s += en * ((gn[n] + gn[-n]) * (1 - math.cosh(n * beta)) + (gn[n] - gn[-n]) * math.sinh(n * beta))
        worst = max(worst, abs(direct - s) / (1 + abs(s)))
    res['pairing_identity_complex_g_max_rel'] = worst
    # parity: int Gamma sin(w/sqrt2) Im F_n(w) dw = 0 when the envelope is even (every state), not when it is shifted (control)
    def parity(shift):
        L2, s2, x = 2.0, 0.3, 0.01; w = np.linspace(-12, 12, 4001); dw = w[1] - w[0]; t = np.linspace(0, 2 * np.pi, 256, endpoint=False)
        n, a, ph = trig(5, rng); a *= 0.2 / np.sum(np.abs(a))
        E = np.exp(-(w - shift) ** 2 / (2 * s2)); R = ev(n, a, ph, t); Rp = ev(n, a, ph, t, 1)
        h = 1 + np.outer(E, R); F = np.outer(E, Rp) ** 2 / h
        Fn = np.mean(F * np.exp(-1j * 3 * t), axis=1)                # n = 3
        Gam = np.exp(-w ** 2 / (2 * L2)) / math.sqrt(2 * math.pi * L2)
        return abs(np.sum(Gam * np.sin(w / math.sqrt(2)) * Fn.imag) * dw) / (np.sum(Gam * np.abs(Fn)) * dw)
    res['parity_even_envelope_rel'] = max(parity(0.0) for _ in range(20)); res['parity_shifted_control_rel'] = parity(0.7)
    # |F_n| <= <F> for F >= 0 (and so |g_n| <= <F> for g = F e^{i psi})
    worst = 0.0
    for _ in range(200):
        t = np.linspace(0, 2 * np.pi, 512, endpoint=False); n, a, ph = trig(8, rng); F = ev(n, a, ph, t) ** 2
        psi = ev(*trig(4, rng), t); g = F * np.exp(1j * psi)
        worst = max(worst, max(abs(np.mean(g * np.exp(-1j * k * t))) for k in range(1, 20)) / np.mean(F))
    res['abs_coeff_over_mean_max'] = worst
    # the two elementary inequalities: cosh u - 1 <= (u^2/2) e^u and sinh u <= u e^u (u >= 0), exact-enough on a grid in decimal
    ok = True
    for k in range(0, 4001):
        u = D(k) / 400
        ok &= (u.exp() + (-u).exp()) / 2 - 1 <= u * u / 2 * u.exp() and (u.exp() - (-u).exp()) / 2 <= u * u.exp()
    res['cosh_sinh_inequalities_decimal'] = ok
    res['ok'] = worst <= 1 + 1e-12 and res['pairing_identity_complex_g_max_rel'] < 1e-10 and res['parity_even_envelope_rel'] < 1e-12 \
        and res['parity_shifted_control_rel'] > 1e-4 and ok
    return res


def check_Z3(rng, ntrial=400):
    res = {}
    # (iv) derivative of (1 + a t)/(1 + t)^2 equals (a - 2 - a t)/(1 + t)^3: exact difference quotient check on rationals (polynomial identity)
    ok = True
    for _ in range(500):
        a, t = abs(rq()), abs(rq())
        f = lambda tt: (1 + a * tt) / (1 + tt) ** 2
        # (f(t+h) - f(t))/h - f'(t) = O(h); use the exact algebraic form: f(t+h) - f(t) = h [a - 2 - a t - h a ... ]; check via two h's
        h1, h2 = Fr(1, 10 ** 6), Fr(1, 10 ** 7)
        d1 = (f(t + h1) - f(t)) / h1; d2 = (f(t + h2) - f(t)) / h2; der = (a - 2 - a * t) / (1 + t) ** 3
        ok &= abs(d2 - der) <= abs(d1 - der) / 9 + Fr(1, 10 ** 30)
    res['iv_derivative_formula'] = ok
    cs = (D('0.0015').exp()) * (1 + D(1).exp().sqrt() + 2 / D(1).exp())
    res['c_s_0.003'] = float(cs); res['mu_threshold_a_le_2'] = float(2 * (cs / 2).ln())
    # the range where t_0 < (a - 2)/a (the monotonicity step does not apply): upper end
    last = 0.0
    for k in range(1, 4000):
        mu = 0.003 + k * 0.0003; t0 = math.e / 2 * (mu - 0.003) * math.exp((mu - 0.003) / 2); a = float(cs) * math.exp(-mu / 2)
        if a > 2 and t0 < (a - 2) / a: last = mu
    res['invalid_range_upper_end'] = last
    # (i)-(iii) on random (P2)/(P3) windows: mean identity, s <= a t, t >= t_0(x), P <= Pi*(mu) for mu >= 2
    import cc_core as C
    bad = 0; tested = 0
    for _ in range(ntrial):
        x = 10 ** rng.uniform(-6, math.log10(0.003)); nmax = int(min(4000, 6 / math.sqrt(x) + 5)); n = np.arange(1, nmax + 1)
        d = np.exp(-(n ** 2 - 1) * x / 2); v = d * rng.uniform(0, 1, nmax) ** rng.uniform(0.2, 5); v[0] = 1.0
        if _ % 2:                                                   # a band at the cap around n^2 x = mu* (to reach mu >= 2 under (P2))
            x = 10 ** rng.uniform(-7, -5); mus = rng.uniform(2.5, 6); n0 = int(math.sqrt(mus / x)); K = int(rng.integers(200, 3000))
            n = np.arange(1, n0 + K + 1); d = np.exp(-(n ** 2 - 1) * x / 2); v = np.zeros(len(n)); v[0] = 1.0
            v[n0:] = d[n0:] * rng.uniform(0.8, 1.0, len(n) - n0)
        y = 0.25 / v.sum() * rng.uniform(0.3, 1); q = y * v
        S0 = np.sum(q ** 2); S2 = np.sum(n ** 2 * q ** 2); mu = x * S2 / S0; P = x * S2
        s = np.sum(v[1:] ** 2); tt = np.sum(v[1:]); a = math.exp(x / 2) * (1 + math.sqrt(math.e) + 2 / math.e) * math.exp(-mu / 2)
        bad += abs(np.sum(v[1:] ** 2 * (n[1:] ** 2 * x - mu)) - (mu - x)) > 1e-9 * (mu + 1)
        bad += s > a * tt * (1 + 1e-12) + 1e-300
        if mu > x: bad += tt < math.e / 2 * (mu - x) * math.exp((mu - x) / 2) * (1 - 1e-12)
        if mu >= 2:
            tested += 1; bad += P > math.exp(C.log_Pistar(mu)) * (1 + 1e-9)
    res['random_windows'] = {'trials': ntrial, 'tested_mu_ge_2': tested, 'violations': int(bad)}
    # Pi*(mu) decreasing for mu >= 2 (used for the left end of a mu cell)
    mus = np.exp(np.linspace(math.log(2), math.log(30000), 20000)); vals = np.array([C.log_Pistar(m) for m in mus])
    res['Pistar_decreasing_on_[2,3e4]'] = bool(np.all(np.diff(vals) < 0))
    res['ok'] = ok and bad == 0 and res['Pistar_decreasing_on_[2,3e4]'] and abs(res['mu_threshold_a_le_2'] - 1.055) < 1e-3 and 0.2 < last < 0.23
    return res


def check_BE():
    res = {}
    # the algebra of the choice of c, exact: e s^2 M/(4c) with c = (b + s^2 J)/(16 s^2) equals 4 e s^4 M/(b + s^2 J);
    # e s^2 2c Phi = e (b + s^2 J) Phi/8
    ok = True
    for _ in range(2000):
        e, s, M, b, J, Phi = (abs(rq()) + Fr(1, 10 ** 6) for _ in range(6))
        c = (b + s * s * J) / (16 * s * s)
        ok &= (e * s * s * M / (4 * c) == 4 * e * s ** 4 * M / (b + s * s * J)) and (e * s * s * 2 * c * Phi == e * (b + s * s * J) * Phi / 8)
    res['choice_of_c_exact'] = ok
    # sigma^2/(b y) = 2/(1-b)^2 with y = (1-b)^2 sigma^2/(2b), exact
    ok2 = True
    for _ in range(1000):
        s, b = abs(rq()) + Fr(1, 10 ** 6), Fr(rnd.randint(1, 998), 1000)
        y = (1 - b) ** 2 * s * s / (2 * b); ok2 &= s * s / (b * y) == 2 / (1 - b) ** 2
    res['sigma2_over_by_exact'] = ok2
    # max_y e^{-y/2}(2 + y) = 2 at y = 0 (derivative -(y/2) e^{-y/2}); decimal grid check
    mx = max(((-D(k) / 200 / 2).exp() * (2 + D(k) / 200)) for k in range(0, 20001))
    res['max_e_y2_times_2py'] = float(mx)
    # rho_W at 60 digits: 48 (1 + Qb) K_1 sigma^2/(1 - b)^2, Lambda = 1e-3/sigma^2, U_1 at ell = 6907.76, sigma = 1e-24, b = 1e4 sigma^2
    Qb = D('0.25'); sig = D('1e-24'); b = 10 ** 4 * sig ** 2; U1 = Qb * (2 * D('6907.76') + D('3.1')).sqrt() + D('1e-6')
    K1s2 = D('1e-3') / (1 - Qb) + sig ** 2 * (U1 ** 2 + D('0.5')) / (1 - Qb) ** 3
    rho = 48 * (1 + Qb) * K1s2 / (1 - b) ** 2
    res['rho_W_60_digits'] = str(rho); res['rho_W_minus_0.08'] = float(rho - D('0.08'))
    # the tail ratio constant 8 g_s T/J^2 <= (8/0.36) T/sigma^2 = 22.2 T/sigma^2 (J >= 0.6 sigma, g_s <= 1, 6/(1 - Qb) = 8)
    res['tail_constant'] = float(D(6) / (1 - Qb) / (D('0.6') ** 2))
    res['ok'] = ok and ok2 and abs(float(mx) - 2) < 1e-30 and D('0.08') < rho <= D('0.0801') and res['tail_constant'] <= 22.3
    return res


def main():
    rng = np.random.default_rng(20261103)
    out = {'S': check_S(rng), 'A4': check_A4(rng), 'Z3': check_Z3(rng), 'BE': check_BE()}
    out['all_ok'] = all(v['ok'] for v in out.values())
    os.makedirs('results', exist_ok=True)
    json.dump(out, open(OUT, 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    print(json.dumps(out, indent=1, default=str))


if __name__ == '__main__':
    main()
