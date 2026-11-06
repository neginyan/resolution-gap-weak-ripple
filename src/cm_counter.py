"""Stage CM 1: is the local Z3' (P_theta <= C Pi*(mu_theta), one universal C) true?  Python 3.11.

Window: f(t) = 1 + a1 cos t + a sum_{n=N1}^{N2} cos(n (t - t0)), t0 = pi/2 (a NODE of the fundamental), a <= a1 ((P2)), f >= 0 checked on a
grid with a Lipschitz margin, K = x (c_b = 1), q_1 = a1 e^{-x/2}, q_n = a e^{-n^2 x/2}, sum q <= Qbar = 1/4 ((P3)).
Focus at t0 with variance tau = x/b.  Exact local quantities (characteristic function of the wrapped normal at integer k is e^{-k^2 tau/2}):
  S0_theta = sum_{n,m} q_n q_m [e^{-(n-m)^2 tau/2} cos(al_n - al_m) + e^{-(n+m)^2 tau/2} cos(al_n + al_m)],  al_n = n t0 - phi_n,
  S2_theta = sum_{n,m} n m q_n q_m [e^{-(n-m)^2 tau/2} cos(al_n - al_m) - e^{-(n+m)^2 tau/2} cos(al_n + al_m)],
mu_theta = x S2/S0, P_theta = x S2.  rho = P_theta/Pi*(mu_theta) with Pi* of cc_core (Lemma Z3', x <= 0.003).
Also the same ratio against the weak bound 2 mu Qbar^2.  The best cases are re-evaluated with 60-digit decimals (Python decimal; cosine by its series).
-> results/cm_counter.json
"""
import json, math, os
import numpy as np
from decimal import Decimal, getcontext
import cc_core as C

QB = 0.25
OUT = 'results/cm_counter.json'


def window(a1, a, N1, N2, t0, x):
    n = np.arange(N1, N2 + 1)
    ns = np.concatenate(([1], n)); phi = np.concatenate(([0.0], n * t0)); araw = np.concatenate(([a1], np.full(len(n), a)))
    q = araw * np.exp(-ns ** 2 * x / 2)
    return ns, phi, araw, q


def f_min(ns, phi, araw, M=1 << 16):
    """min of the raw f = 1 + a1 cos t + a D(t - t0) on a grid of M points, minus the midpoint margin (h/2)^2 max|f''|/2, max|f''| <= sum n^2 a_n.
    D(u) = sum_{n=N1}^{N2} cos(n u) = cos((N1+N2) u/2) sin(W u/2)/sin(u/2), W = N2 - N1 + 1 (closed form; exact)."""
    a1 = araw[0]; a = araw[1]; N1 = int(ns[1]); N2 = int(ns[-1]); W = N2 - N1 + 1; t0 = phi[1] / N1
    t = np.linspace(0, 2 * math.pi, M, endpoint=False); u = t - t0
    su = np.sin(u / 2); small = np.abs(su) < 1e-12
    Dk = np.where(small, float(W), np.cos((N1 + N2) * u / 2) * np.sin(W * u / 2) / np.where(small, 1.0, su))
    f = 1.0 + a1 * np.cos(t) + a * Dk
    l2 = float(np.sum(ns.astype(float) ** 2 * araw)); h = 2 * math.pi / M
    return float(f.min()) - (h / 2) ** 2 * l2 / 2


def local(ns, phi, q, t0, x, tau):
    al = ns * t0 - phi
    D = ns[:, None] - ns[None, :]; S = ns[:, None] + ns[None, :]
    cd = np.cos(al[:, None] - al[None, :]); cs = np.cos(al[:, None] + al[None, :])
    qq = q[:, None] * q[None, :]; nm = ns[:, None] * ns[None, :]
    ed = np.exp(-D ** 2 * tau / 2); es = np.exp(-S ** 2 * tau / 2)
    S0 = float(np.sum(qq * (ed * cd + es * cs))); S2 = float(np.sum(nm * qq * (ed * cd - es * cs)))
    return S0, S2


def ratio(S0, S2, x):
    mu = x * S2 / S0; P = x * S2
    return mu, P, P / math.exp(C.log_Pistar(mu)), P / (2 * mu * QB * QB)


def dcos(z):
    """cos by its series after reduction mod 2 pi (60-digit decimal)."""
    PI = Decimal('3.14159265358979323846264338327950288419716939937510582097494459')
    z = z % (2 * PI)
    if z > PI: z -= 2 * PI
    term = Decimal(1); s = Decimal(1); z2 = z * z; k = 0
    while True:
        k += 2; term = -term * z2 / (k * (k - 1)); s += term
        if abs(term) < Decimal('1e-70'): return s


def mp_check(a1, a, N1, N2, t0, x, tau):
    getcontext().prec = 60
    D = lambda v: Decimal(repr(float(v)))
    PI = Decimal('3.14159265358979323846264338327950288419716939937510582097494459')
    x = D(x); tau = D(tau); t0 = PI / 2
    ns = [1] + list(range(N1, N2 + 1)); phi = [Decimal(0)] + [n * t0 for n in range(N1, N2 + 1)]
    araw = [D(a1)] + [D(a)] * (N2 - N1 + 1)
    q = [ar * (-Decimal(n) ** 2 * x / 2).exp() for n, ar in zip(ns, araw)]
    al = [n * t0 - p for n, p in zip(ns, phi)]
    S0 = Decimal(0); S2 = Decimal(0)
    for i, n in enumerate(ns):
        for j, m in enumerate(ns):
            ed = (-Decimal(n - m) ** 2 * tau / 2).exp(); es = (-Decimal(n + m) ** 2 * tau / 2).exp()
            cd = dcos(al[i] - al[j]); cs = dcos(al[i] + al[j])
            S0 += q[i] * q[j] * (ed * cd + es * cs); S2 += n * m * q[i] * q[j] * (ed * cd - es * cs)
    mu = x * S2 / S0; P = x * S2
    E_ = Decimal(1).exp(); Q = Decimal(QB)
    cs_ = Decimal('0.0015').exp() * (1 + E_.sqrt() + 2 / E_)
    t0_ = E_ / 2 * (mu - Decimal('0.003')) * ((mu - Decimal('0.003')) / 2).exp(); aa = cs_ * (-mu / 2).exp()
    Pi = min(Q ** 2 * mu * (1 + aa * t0_) / (1 + t0_) ** 2, mu * Q ** 2)
    return {'mu_theta': float(mu), 'P_theta': float(P), 'rho': float(P / Pi), 'sum_q': float(sum(q))}


def main():
    res = {'note': 'local Z3 prime test: coherent packet at a node of the fundamental; exact formulas; numerical', 'cases': []}
    t0 = math.pi / 2
    best = None
    # scan: x, band, a1 (a = a1, the largest the packet may be under (P2)); a chosen so that f >= 0 and sum q <= Qbar
    for x in (0.003, 1e-3, 3e-4, 1e-4, 1e-5, 1e-6):
        for mu0 in (3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0):
            N = int(round(math.sqrt(mu0 / x)))
            for width in (1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40):
                N1 = N - width // 2; N2 = N1 + width - 1
                if N1 < 2: continue
                for a1 in (0.1, 0.15, 0.2, 0.25, 0.3, 0.4):
                    # largest a <= a1 with f >= 0: the packet's trough near t0 is about -a*width; bisection on the grid minimum
                    lo, hi = 0.0, a1
                    for _ in range(22):
                        a = 0.5 * (lo + hi)
                        ns, phi, araw, q = window(a1, a, N1, N2, t0, x)
                        if f_min(ns, phi, araw) >= 0 and q.sum() <= QB: lo = a
                        else: hi = a
                    a = lo * (1 - 1e-9)
                    if a <= 0: continue
                    ns, phi, araw, q = window(a1, a, N1, N2, t0, x)
                    for b in (1.0, 0.5, 0.2, 0.1):
                        tau = x / b
                        S0, S2 = local(ns, phi, q, t0, x, tau)
                        if S0 <= 0: continue
                        mu, P, rho, rho_w = ratio(S0, S2, x)
                        case = {'x': x, 'mu0': mu0, 'N1': int(N1), 'N2': int(N2), 'a1': a1, 'a': a, 'b': b, 'tau': tau, 'sum_q': float(q.sum()),
                                'f_min': f_min(ns, phi, araw), 'mu_theta': mu, 'P_theta': P, 'rho': rho, 'rho_weak': rho_w}
                        if mu >= 2 and (best is None or rho > best['rho']): best = case
                        if rho > 10: res['cases'].append(case)
    res['best'] = best
    res['best_mp60'] = mp_check(best['a1'], best['a'], best['N1'], best['N2'], t0, best['x'], best['tau'])
    # the two headline cases: x = 0.003 (the largest x of the class) and a small x
    heads = []
    for x in (0.003, 1e-5):
        cand = [c for c in res['cases'] if c['x'] == x]
        if cand:
            c = max(cand, key=lambda c: c['rho']); c = dict(c); c['mp60'] = mp_check(c['a1'], c['a'], c['N1'], c['N2'], t0, c['x'], c['tau']); heads.append(c)
    res['headline'] = heads
    # rho as a function of mu_theta (the largest rho found in bins of mu_theta), to show the growth
    bins = {}
    for c in res['cases']:
        k = int(c['mu_theta'])
        if k not in bins or c['rho'] > bins[k]['rho']: bins[k] = c
    res['rho_by_mu'] = {str(k): {'rho': v['rho'], 'x': v['x'], 'b': v['b'], 'mu_theta': v['mu_theta']} for k, v in sorted(bins.items())}
    os.makedirs('results', exist_ok=True)
    json.dump(res, open(OUT, 'w'), indent=1)
    print('cases with rho > 10:', len(res['cases']))
    print('best', json.dumps(best, indent=1)); print('mp60', res['best_mp60'])
    for h in heads: print('headline x=%g: mu_theta %.3f rho %.2f (mp60 %.2f) sum_q %.4f f_min %.2e' % (h['x'], h['mu_theta'], h['rho'], h['mp60']['rho'], h['sum_q'], h['f_min']))
    for k, v in res['rho_by_mu'].items(): print('mu_theta ~', k, 'rho_max %.1f' % v['rho'], 'x', v['x'], 'b', v['b'])


if __name__ == '__main__':
    main()
