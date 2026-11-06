"""Stage CO 2.1: Lemma R_theta (the local route with the F term of Lemma CN-2 and the remainders N1-N16 of co_remainders.py) as a lower bound
of (a) + (a') = (2 A_w + N_phi)/(sigma den), checked against the direct two-dimensional quadrature of ca_reduction.direct on random states with a
NARROW focus (tau from 0.3 to 20) -- the numerical check of the derivation (signs, factors, the local identity).  Python 3.11.
Here the window quantities are the actual ones of the window (eta = sup|R|, N1, N2, E0, E1, kappa, R_theta, E1^{(2tau)} exact), so this checks the
inequalities, not the window-generic constants (those are checked in co_remainders.py / cl_link).
-> results/co_rtheta.json
"""
import json, math, sys
import numpy as np
from ca_dec import D, exp, log, sqrt, expm1, log1p, PI, E as EE
import ca_reduction as CA
import cm_local_bound as LB
SQ2 = sqrt(2)


def local_numbers(q0, ph, x, tau):
    """E0, E1, kappa, R_theta, E1^{(2 tau)} at the focus t0 = 0 (the state of ca_reduction.direct: R = sum q cos(n t - phi))."""
    ns = np.arange(1, len(q0) + 1, dtype=float); q = np.asarray(q0, float); al = -np.asarray(ph, float)
    E0, E1 = LB.local_quantities(ns, q, al, x, tau); _, E1_2 = LB.local_quantities(ns, q, al, x, 2 * tau)
    s = np.linspace(-12 * math.sqrt(tau), 12 * math.sqrt(tau), 400001); rho = np.exp(-s ** 2 / (2 * tau)) / math.sqrt(2 * math.pi * tau) * (s[1] - s[0])
    R = sum(qq * np.cos(n * s - p) for n, qq, p in zip(ns, q, ph)); Rp = sum(-n * qq * np.sin(n * s - p) for n, qq, p in zip(ns, q, ph))
    kap = float((rho * (s * s / tau) * Rp * Rp).sum()) / E1; Rth = float((rho * R).sum())
    return E0, E1, kap, Rth, E1_2


def M2_const(beta, tau):
    u = np.linspace(0, 400, 400001); c = float(beta / sqrt(tau))
    return float(np.max(u * np.exp(-u / 2) * (1 + u * np.cosh(c * np.sqrt(u))))) * (1 + 1e-6)


def lemma_R_theta(sig, x, b, L2, s2, eta, Qb, N1, N2, E0, E1, kap, Rth, E1_2, detail=False):
    """Num_lo/(sigma den_up) of the local route; all the remainders of co_remainders.py with the window's own quantities."""
    sig, x, b, L2, s2, eta, Qb, N1, N2, E0, E1, kap, Rth, E1_2 = (D(t) for t in (sig, x, b, L2, s2, eta, Qb, N1, N2, E0, E1, kap, Rth, E1_2))
    tau = x / b; y = (1 - b) ** 2 * sig ** 2 / (2 * b); e = exp(-y / 2); beta = sqrt(tau * y)
    s = sqrt(s2); L = sqrt(L2); r = s2 / L2; se2 = s2 / (2 + r); gs = sqrt(se2 / L2); kE = sqrt(r / (1 + r)); seE = s2 / (1 + r)
    Wt = x * gs * E1; D0 = E0 * gs * se2 / s2 ** 2
    xi = sig * sqrt(x) * N1 / (1 - Qb); zeta = sig ** 2 * eta * exp(D('-0.5')) / (s * (1 - eta)); xie = xi / SQ2
    c0 = (xi + zeta) / SQ2; c1 = sig ** 2 / (SQ2 * L2)
    rC = (xi + zeta + sig ** 2 * sqrt(2 / PI) / L) / SQ2; rCE = kE * (xi + zeta + sig ** 2 * sqrt(2 * seE / PI) / L2) / SQ2
    gP = (x * (N2 * (1 + Qb) + N1 ** 2) + sqrt(x) * N1 * (1 + 2 * Qb) / sqrt(EE)) / (1 - Qb) ** 2
    eL = exp(-L2 / 4); eSE = exp(-seE / 4); ese2 = exp(-se2 / 4)
    Psi0 = -expm1(log1p(-y) - y / 2 - L2 / 4)
    aRth = abs(Rth)                                           # the proof uses |R_theta| <= eta; here the actual value (<= eta) for the Z-pieces
    # ---- main (local sandwiches)
    Jw_lo = 1 / L2 + D0 / (1 + eta); Jw_up = 1 / L2 + D0 / (1 - eta)
    Lam1 = aRth * kE * (2 + r) / ((1 + r) * L2)
    twoAw_main = 2 * sig ** 2 * (Jw_lo - Lam1)
    # A_w corrections (phi <= (v'^2 + R_z^2)/2, R_z^2 <= 2 w^2 + 4 chi^2)
    J2 = (1 / L2) * (1 + eta * kE * r / (1 + r)) + 2 * eta * kE * seE / (L2 * s2) + eta ** 2 * gs * se2 / (s2 ** 2 * (1 - eta))
    Jup = (1 / L2) * (1 + eta * kE * r / (1 + r)) + 2 * eta * kE * seE / (L2 * s2) + D0 / (1 - eta)
    M4 = 3 * (1 + eta * kE * r ** 2 / (1 + r) ** 2) + 6 * eta * kE * r / (1 + r) ** 2 + 3 * se2 * D0 / (1 - eta)
    corrA = sig ** 4 / (4 * b) * (1 - b) ** 2 * J2 + sig ** 2 / 2 * M4 + sig ** 2 * D('1.5') * ((xi ** 2 + zeta ** 2) * Jup + sig ** 4 / L2 ** 2 * M4)
    twoAw = twoAw_main - corrA
    # ---- (B) b-term
    Bt = b * (Psi0 - aRth * kE * Psi0 - eta * kE * ((1 - y) * (eSE - e * eL) + y / 2) - (1 + eta) * (1 - y) * rC)
    # ---- (F) Lemma CN-2 (chibar_e: the F_rem form of Lemma R with the local W)
    Frem = Wt / (1 - eta) * ((c0 * sqrt(2 * se2 / PI) + c1 * se2) / SQ2 + c0 ** 2 + c1 ** 2 * se2)
    Ft = Wt / (1 + eta) * ((-expm1(-se2 / 4)) + ese2 * ((y / 2) * kap - xie * sqrt(y * kap) - D('0.51') * y ** 2 * E1_2 / E1)) - Frem
    # ---- (X) cross term
    Ggt = x * (N2 * kE + N1 ** 2 * gs / (1 - eta)) / (1 - eta) + sig * sqrt(x) * N1 * kE * sqrt(2 * seE / PI) * (1 + eta / (1 - eta)) / (s2 * (1 - eta))
    GgtPhi = (x * N2 * kE * (1 - e * eSE) / (1 - eta) + x * N1 ** 2 * gs * (1 - e * ese2) / (1 - eta) ** 2
              + sig * sqrt(x) * N1 * (1 + eta / (1 - eta)) / (s2 * (1 - eta)) * ((1 - e) * kE * sqrt(2 * seE / PI) + e * kE * sqrt(2 / PI) * seE * sqrt(seE) / 2))
    GgtChi = Ggt * (xi + zeta) / SQ2 + sig ** 2 * gP * sqrt(2 / PI) / (SQ2 * L)
    Xt = -sig ** 2 * (1 - b) * (Wt / (1 - eta) + D('0.184') * b * kE) - sig ** 3 * (1 - b) * sqrt(x * E1) * kE * sqrt(2 * seE / PI) / (s2 * (1 - eta)) \
         - sig ** 2 * (1 - b) * (1 + eta) * (GgtPhi + GgtChi)
    # ---- (L) last line
    M2 = D(M2_const(beta, tau))
    Lt = -(sqrt(1 + y) * (sqrt(2 * x * E1 * b) * rCE + sig * sqrt(b) * (1 + eta) * Ggt) + D('0.184') * b * kE * (1 - eSE) + (y / 2) * b * Qb * kE * M2)
    twoN = Bt + Ft + Xt + Lt
    Num = twoAw + twoN / 2
    den = b * (1 + aRth * kE) + Wt / (1 - eta) + D('0.184') * b * kE + sig ** 2 * (Jw_up + Lam1)
    out = {'lower': Num / (sig * den), 'Num': Num, 'den': den}
    if detail: out.update({'twoAw_main': twoAw_main, 'corrA': corrA, 'B': Bt, 'F': Ft, 'X': Xt, 'L': Lt, 'Frem': Frem})
    return out


def packet_window(rng, n0, W, where, amp=0.12, tot=0.24):
    """fundamental + a flat packet of harmonics n0..n0+W at the node / crest / trough of the fundamental (focus at t = 0): q0 list by harmonic index."""
    q0 = [0.0] * (n0 + W); ph = [0.0] * (n0 + W)
    q0[0] = tot - amp; ph[0] = {'node': math.pi / 2, 'crest': 0.0, 'trough': math.pi}[where]
    for n in range(n0, n0 + W): q0[n - 1] = amp / W
    return q0, ph


def run(nstates=30, seed=7, npack=6):
    rng = np.random.default_rng(seed); out = []
    for i in range(nstates * 3 + npack):
        if len(out) >= nstates + npack: break
        sig = 10 ** rng.uniform(-1.6, -0.9); tau = 10 ** rng.uniform(-0.5, 1.3); x = 10 ** rng.uniform(-1.5, 0.3)
        b = x / tau; L2 = 10 ** rng.uniform(-2.5, -0.7); s2 = L2 / 4 * 10 ** rng.uniform(-1.3, 0)
        if s2 < sig ** 2: s2 = sig ** 2 * 1.5
        if 4 * s2 > L2: L2 = 4 * s2
        if i < npack:
            n0, W = ((6, 4), (10, 6), (15, 5))[i % 3]; where = ('node', 'crest', 'trough')[i % 3]
            q0, ph = packet_window(rng, n0, W, where, amp=0.12 if i < 3 else 0.2); tau = (0.3, 1.0, 3.0)[i % 3] * (1 if i < 3 else 2); x = 10 ** rng.uniform(-1.5, -0.5); b = x / tau
        else:
            nh = rng.integers(1, 5); q0 = list(rng.uniform(0.02, 0.6 / nh, nh)); ph = list(rng.uniform(0, 2 * math.pi, nh))
        wn = CA.window_numbers(q0, ph, x)
        if wn['eta'] >= 0.9 or b >= 0.5: continue
        kp = math.sqrt(x) / sig; V = sig ** 2 / b
        nv = int(min(80000, max(6000, 40 * len(q0) * kp * 20 * math.sqrt(V) / (2 * math.pi))))
        d = CA.direct(sig, x, b, L2, s2, q0, ph, nv=nv, nw=1500)
        E0, E1, kap, Rth, E1_2 = local_numbers(q0, ph, x, tau)
        lr = lemma_R_theta(sig, x, b, L2, s2, wn['eta'], wn['Qb'], wn['N1'], wn['N2'], E0, E1, kap, Rth, E1_2, detail=True)
        R = CA.lemma_R(sig, x, b, L2, s2, wn['S0'], wn['S2'], wn['eta'], wn['Qb'], wn['N1'], wn['N2'])
        rec = {'sigma': sig, 'tau': tau, 'x': x, 'b': b, 'L2': L2, 's2': s2, 'q0': q0, 'phases': ph, 'eta': wn['eta'], 'E0': E0, 'E1': E1, 'kappa': kap, 'R_theta': Rth,
               'direct': d['sum'], 'lemma_R_theta': float(lr['lower']), 'ratio_theta': float(lr['lower']) / d['sum'], 'lemma_R_period': float(R['lower_R']),
               'ratio_period': float(R['lower_R']) / d['sum'], 'ok': float(lr['lower']) <= d['sum'] * (1 + 1e-9),
               'pieces': {k: float(v) for k, v in lr.items() if k not in ('lower',)}}
        out.append(rec)
        print('%2d sig %.3f tau %6.2f x %.3f b %.3f eta %.3f kappa %.2f Rth %+.3f: direct %.5f  R_theta %.5f (ratio %.4f)  R(period) %.5f (ratio %.3f) %s' % (
            len(out), sig, tau, x, b, wn['eta'], kap, Rth, d['sum'], rec['lemma_R_theta'], rec['ratio_theta'], rec['lemma_R_period'], rec['ratio_period'], 'ok' if rec['ok'] else 'FAIL'), flush=True)
    res = {'states': out, 'all_ok': all(t['ok'] for t in out), 'max_ratio_theta': max(t['ratio_theta'] for t in out), 'min_ratio_theta': min(t['ratio_theta'] for t in out)}
    json.dump(res, open('results/co_rtheta.json', 'w'), indent=1, default=float)
    print('all ok:', res['all_ok'], ' ratio_theta in [%.4f, %.4f]' % (res['min_ratio_theta'], res['max_ratio_theta']))
    return res


if __name__ == '__main__':
    run(int(sys.argv[1]) if len(sys.argv) > 1 else 30)
