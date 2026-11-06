"""Stage CO 2.4: adversarial search against the local remainder bounds (I1)-(I10) of CO_探索計画.md, with the deliberately wrong versions (W1)-(W7).
Python 3.11.  Local averages E_rho[g] = int N(0,tau)(s) g(t0 + s) ds by trapezoid on the line (grid resolving the highest harmonic).
Windows: f = 1 + sum a_n cos(n t - phi_n) >= 0 (packets with sum |a_n| <= 1; random bands checked on a grid), blurred q_n = a_n e^{-n^2 x/2}.
-> results/co_search.json
"""
import json, math, sys
import numpy as np

QB_ = 0.25


class Win:
    def __init__(self, ns, a, phi, x=1e-4):
        self.ns = np.asarray(ns, float); self.a = np.asarray(a, float); self.phi = np.asarray(phi, float); self.x = x
        self.q = self.a * np.exp(-self.ns ** 2 * x / 2)
        t = np.linspace(0, 2 * math.pi, 20001)
        self.eta = float(np.max(np.abs(self.R(t)))); self.Qb = float(self.q.sum())

    def R(self, t):
        t = np.asarray(t, float); return (self.q[:, None] * np.cos(self.ns[:, None] * t[None, :] - self.phi[:, None])).sum(0)

    def Rp(self, t):
        t = np.asarray(t, float); return (-self.ns[:, None] * self.q[:, None] * np.sin(self.ns[:, None] * t[None, :] - self.phi[:, None])).sum(0)

    def Rpp(self, t):
        t = np.asarray(t, float); return (-self.ns[:, None] ** 2 * self.q[:, None] * np.cos(self.ns[:, None] * t[None, :] - self.phi[:, None])).sum(0)


def grid(tau, nmax, span=10.0):
    h = min(2 * math.pi / (16 * nmax), math.sqrt(tau) / 50)
    m = int(2 * span * math.sqrt(tau) / h) + 1
    s = np.linspace(-span * math.sqrt(tau), span * math.sqrt(tau), m)
    rho = np.exp(-s ** 2 / (2 * tau)) / math.sqrt(2 * math.pi * tau)
    return s, rho * (s[1] - s[0])


def loc(win, t0, tau):
    s, w = grid(tau, win.ns.max())
    R = win.R(t0 + s); Rp = win.Rp(t0 + s); Rpp = win.Rpp(t0 + s)
    E = lambda g: float((w * g).sum())
    return s, w, R, Rp, Rpp, E


def one_state(win, t0, tau, beta, rng):
    """ratios lhs/rhs of (I1)-(I9) and of the wrong versions at one (window, focus, tau, beta)."""
    s, w, R, Rp, Rpp, E = loc(win, t0, tau)
    eta, Qb, x = win.eta, win.Qb, win.x
    y = beta ** 2 / tau; e = math.exp(-y / 2)
    E0 = E(R * R); E1 = E(Rp * Rp); Rth = E(R)
    out = {}
    # (I1)
    out['I1'] = abs(Rth) / min(eta, Qb * math.exp(-tau / 2))
    # (I2), (W1)
    v = abs(E((s * s / tau) * R)); out['I2'] = v / eta; out['W1'] = v / (eta / 2)
    # (I3), (W2)
    v = abs(E((1 - np.cos(beta * s / tau)) * R)); out['I3'] = v / ((y / 2) * eta); out['W2'] = v / ((y / 2) * eta / 4)
    # (I4), (W6): |E[(1 - cos) R'']| <= (y/2) (Qb/x) b M2 with b = x/tau, i.e. <= (y/2) Qb M2 / tau
    u = np.linspace(0, 200, 200001); M2 = float(np.max(u * np.exp(-u / 2) * (1 + u * np.cosh(beta * np.sqrt(u / tau)))))
    v = abs(E((1 - np.cos(beta * s / tau)) * Rpp)); out['I4'] = v / ((y / 2) * Qb * M2 / tau); out['W6'] = v / ((y / 2) * Qb * (2 / math.e) / tau)
    # (I6), (W4): |E[e^{i beta s/tau}(s/beta - i) R']| <= (1/beta) sqrt(E1 (tau + beta^2))
    v = abs(complex((w * (np.exp(1j * beta * s / tau) * (s / beta - 1j) * Rp)).sum())); out['I6'] = v / (math.sqrt(E1 * (tau + beta ** 2)) / beta); out['W4'] = v / (math.sqrt(E1 * tau / 4) / beta)
    # (I7): E[cos(beta s/tau) R] = Rth - delta_c, |delta_c| <= (y/2) eta  (same as I3); E[sin(beta s/tau) R] finite: record
    out['I7_sin'] = abs(E(np.sin(beta * s / tau) * R))
    # (I9), (W5): E1^{(2tau)} vs the bound (Qb-free form: E1^{(2tau)} <= (S2/2) theta_{2tau}(0) with S2 = sum n^2 q^2 ... here the x-free version)
    s2_, w2_ = grid(2 * tau, win.ns.max()); E1_2 = float((w2_ * win.Rp(t0 + s2_) ** 2).sum())
    S2 = float(np.sum(win.ns ** 2 * win.q ** 2)); th0 = 1 + math.sqrt(math.pi / tau)
    out['I9'] = E1_2 / ((S2 / 2) * th0); out['W5'] = E1_2 / E1 if E1 > 0 else 0.0
    # (I5), (W3): identity E[h d_t(h_w/h)] = E_w E[R'/h] with a random envelope value; and the bound |E_w| sqrt(E1)/(1-eta)
    Ew_over_E = -rng.uniform(0.1, 3.0); Eval = rng.uniform(0.2, 1.0); Ew = Ew_over_E * Eval
    h = 1 + Eval * R; hw = Ew * R; hwt = Ew * Rp; ht = Eval * Rp
    lhs = E(hwt - hw * ht / h); rhs = Ew * E(Rp / h)
    out['I5_id'] = abs(lhs - rhs) / (abs(Ew) * E(np.abs(Rp / h)) + 1e-300)          # relative to the size of the terms (rhs itself may cancel)
    out['I5'] = abs(lhs) / (abs(Ew) * math.sqrt(E1) / (1 - eta)) if E1 > 0 else 0.0
    out['W3'] = abs(lhs + rhs) / (abs(rhs) + 1e-300) if abs(rhs) > 1e-14 else 0.0        # wrong sign: should be ~2, i.e. far from 0
    out.update({'E0': E0, 'E1': E1, 'kappa': E((s * s / tau) * Rp * Rp) / E1 if E1 > 0 else None, 'Rth': Rth, 'eta': eta, 'tau': tau, 'beta': beta})
    return out


def I10_check(win, t0, tau, beta, b, s2, L2, rng):
    """(I10): the b-term lower bound against the exact value, with chi = 0 (the chi part is the pointwise bound |chi| <= chibar of paper F).
    exact: int Gamma b (E[h] - (1-y) Re E[e^{i beta s/tau} h e^{i w/sqrt2}]);  bound: b[Z Psi0 - eta k_E((1-y)(e^{-se_E/4} - e e^{-L2/4}) + y/2)]."""
    s, w, R, Rp, Rpp, E = loc(win, t0, tau)
    y = beta ** 2 / tau; e = math.exp(-y / 2); r = s2 / L2; kE = math.sqrt(r / (1 + r)); seE = s2 / (1 + r)
    wg = np.linspace(-8 * math.sqrt(L2), 8 * math.sqrt(L2), 4001); G = np.exp(-wg ** 2 / (2 * L2)) / math.sqrt(2 * math.pi * L2) * (wg[1] - wg[0])
    Ev = np.exp(-wg ** 2 / (2 * s2))
    cosb = E(np.cos(beta * s / tau) * R); sinb = E(np.sin(beta * s / tau) * R); Rth = E(R)
    # for each w: E[h] = 1 + E Rth;  Re E[e^{i beta s/tau} h e^{iw/sqrt2}] = cos(w/sqrt2)(e + E cosb) - sin(w/sqrt2)(E sinb)
    val = (G * ((1 + Ev * Rth) - (1 - y) * (np.cos(wg / math.sqrt(2)) * (e + Ev * cosb) - np.sin(wg / math.sqrt(2)) * Ev * sinb))).sum()
    Z = 1 + kE * Rth; Psi0 = 1 - (1 - y) * e * math.exp(-L2 / 4)
    bound = Z * Psi0 - win.eta * kE * ((1 - y) * (math.exp(-seE / 4) - e * math.exp(-L2 / 4)) + y / 2)
    wrong = Z * Psi0                                                                             # (W7): the eta k_E term dropped
    return {'I10': bound / val if val > 0 else None, 'W7': wrong / val if val > 0 else None, 'val': val, 'bound': bound, 'Rth': Rth}


def packets(rng):
    fams = []
    for n0, W in ((10, 4), (40, 10), (100, 30), (200, 60), (400, 100)):
        for shape in ('flat', 'gauss'):
            ns = np.arange(n0, n0 + W + 1)
            a = np.ones(len(ns)) if shape == 'flat' else np.exp(-(ns - n0 - W / 2) ** 2 / (2 * (W / 6.7) ** 2))
            for where in ('node', 'crest', 'trough'):
                for amp in (0.45, 0.9):
                    aa = amp * a / a.sum(); phi = np.zeros(len(ns))
                    if where == 'crest': ns2 = np.concatenate(([1], ns)); a2 = np.concatenate(([1 - amp], aa)); phi2 = np.zeros(len(ns2))
                    elif where == 'trough': ns2 = np.concatenate(([1], ns)); a2 = np.concatenate(([1 - amp], aa)); phi2 = np.concatenate(([math.pi], phi))
                    else: ns2 = np.concatenate(([1], ns)); a2 = np.concatenate(([1 - amp], aa)); phi2 = np.concatenate(([math.pi / 2], phi))
                    fams.append(('packet %s n%d W%d %s amp%.2f' % (shape, n0, W, where, amp), Win(ns2, a2, phi2)))
    # two packets, combs
    ns = np.concatenate(([1], np.arange(30, 41), np.arange(80, 91))); a = np.concatenate(([0.2], 0.4 * np.ones(11) / 11, 0.4 * np.ones(11) / 11))
    fams.append(('two packets aligned', Win(ns, a, np.zeros(len(ns)))))
    phi = np.concatenate(([0.0], np.zeros(11), math.pi * np.ones(11))); fams.append(('two packets opposed', Win(ns, a, phi)))
    ns = np.arange(1, 201); a = np.where(ns % 7 == 1, 1.0, 0.0); a = 0.9 * a / a.sum(); fams.append(('comb 7', Win(ns, a, np.zeros(200))))
    return fams


def random_band(rng, nmax=120):
    for _ in range(100):
        k = rng.integers(2, 40); ns = np.sort(rng.choice(np.arange(1, nmax + 1), k, replace=False)); ns[0] = 1
        a = rng.uniform(0, 1, k); a = a / a.sum() * rng.uniform(0.3, 1.0); phi = rng.uniform(0, 2 * math.pi, k)
        w = Win(ns, a, phi)
        t = np.linspace(0, 2 * math.pi, 40001)
        if np.min(1 + (a[:, None] * np.cos(ns[:, None] * t[None, :] - phi[:, None])).sum(0)) >= 0: return w
    return w


def main(seed=20261104, nrand=150):
    rng = np.random.default_rng(seed)
    fams = packets(rng) + [('random band %d' % i, random_band(rng)) for i in range(nrand)]
    keys = ['I1', 'I2', 'I3', 'I4', 'I5', 'I5_id', 'I6', 'I9', 'W1', 'W2', 'W3', 'W4', 'W5', 'W6']
    worst = {k: (-1.0, None) for k in keys}; n = 0
    for name, win in fams:
        nmax = win.ns.max()
        taus = [1e-4, 1e-3, 1e-2, 0.1, 1.0, 5.0, 20.0] + list(4.0 / nmax ** 2 * np.array([0.25, 1, 4]))
        for tau in taus:
            if tau > 20: continue
            for off in (0.0, 1.0, 2.0, 5.0):
                t0 = off * math.sqrt(tau)
                for beta in (1e-6 * tau, 0.3 * math.sqrt(tau)):
                    o = one_state(win, t0, tau, beta, rng); n += 1
                    for k in keys:
                        v = o.get(k)
                        if v is not None and v > worst[k][0]: worst[k] = (v, {'window': name, 'tau': tau, 'off': off, 'beta': beta, 'kappa': o['kappa'], 'eta': o['eta']})
    # (I10) on random (b, s2, L2) and the windows with the largest |R_theta|
    w10 = (-1.0, None); w7 = (-1.0, None)
    for name, win in fams[:40] + fams[-30:]:
        for tau in (1e-3, 0.1, 1.0, 10.0):
            for off in (0.0, 1.0, 3.0):
                for _ in range(2):
                    b = 10 ** rng.uniform(-6, -1); s2 = 10 ** rng.uniform(-4, 0); L2 = s2 / 10 ** rng.uniform(-3, math.log10(0.25)); L2 = max(L2, 4 * s2)
                    sig2 = b * 10 ** rng.uniform(-3, 0); y = (1 - b) ** 2 * sig2 / (2 * b); beta = math.sqrt(tau * y)
                    o = I10_check(win, off * math.sqrt(tau), tau, beta, b, s2, L2, rng)
                    if o['I10'] is not None and o['I10'] > w10[0]: w10 = (o['I10'], {'window': name, 'tau': tau, 'off': off, 'b': b, 's2': s2, 'L2': L2, 'y': y, 'Rth': o['Rth']})
                    if o['W7'] is not None and o['W7'] > w7[0]: w7 = (o['W7'], {'window': name, 'tau': tau, 'off': off, 'b': b, 's2': s2, 'L2': L2, 'y': y, 'Rth': o['Rth']})
    worst['I10'] = w10; worst['W7'] = w7
    res = {'n_states': n, 'n_windows': len(fams), 'worst': {k: {'max_ratio': v[0], 'at': v[1]} for k, v in worst.items()}}
    ok = all(worst[k][0] <= 1 + 1e-9 for k in ('I1', 'I2', 'I3', 'I4', 'I5', 'I6', 'I9', 'I10')) and worst['I5_id'][0] < 1e-6
    broken = {k: worst[k][0] > 1 for k in ('W1', 'W2', 'W4', 'W5', 'W6', 'W7')}; broken['W3'] = worst['W3'][0] > 1.5
    res['all_true_ok'] = ok; res['wrong_versions_broken'] = broken
    for k in keys + ['I10', 'W7']: print('%-6s max ratio %.4g  at %s' % (k, worst[k][0], worst[k][1]), flush=True)
    print('true inequalities ok:', ok, ' wrong versions broken:', broken)
    json.dump(res, open('results/co_search.json', 'w'), indent=1, default=float)
    return res


if __name__ == '__main__':
    main(nrand=int(sys.argv[1]) if len(sys.argv) > 1 else 150)
