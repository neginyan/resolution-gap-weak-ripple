"""Stage CN 2.2: the two-focus lemma (CN-1), its corollary (the proved joint bound Psi_2), and independent checks of every step of
Lemmas CN-1 and CN-2 (exact arithmetic, 60-digit decimals, quadrature on random windows).  Python 3.11.

Lemma CN-1 (two foci).  With E1 = <theta_tau R'^2>, kappa = <(s^2/tau) theta_tau R'^2>/E1 (unfolded to the line), and delta^2 > kappa tau:
    E1 <= sqrt(tau'/tau) e^{delta^2/(2 tau')} E1^{(tau')} / (1 - kappa tau/delta^2)      for every tau' > 0.
Corollary (Psi_2).  With delta^2 = K kappa tau, tau' = L delta^2, b' = b/(L K kappa), A = sqrt(L K kappa) e^{1/(2L)}/(1 - 1/K), and CM-1 at tau',
E0^{(tau')} <= min(E0/sqrt(L K kappa), 1/16):  x E1 <= A [ (1 + b'/2) E0^ ln(2 B0(E0^)/E0^) + 0.0325 b' ],  E0^ := min(E0/sqrt(LK kappa), 1/16).
Psi_2(mu, b, kappa) := the largest P = 2 mu E0 satisfying this, minimised over K in {2..64}, L in {1..16} with L K kappa <= 1 (tau' <= tau).
-> results/cn_lemmas.json ; exit 1 if a check fails.
"""
import json, math, sys
from decimal import Decimal, getcontext
import numpy as np
import cm_local_bound as LB
import cn_search as S

QB = 0.25
getcontext().prec = 60
res = {}; ok_all = True


def rec(name, ok, **kw):
    global ok_all
    ok_all &= bool(ok); res[name] = {'ok': bool(ok), **kw}
    print('[%s] %s' % ('OK ' if ok else 'FAIL', name), {k: (float(v) if isinstance(v, (int, float, np.floating)) else v) for k, v in kw.items()}, flush=True)


def B0(E0, b): return (1 + math.sqrt(E0) + (E0 / 2) ** (1.0 / (2 + b))) ** 2


def rhs2(P, mu, b, kap, K, L):
    E0 = P / (2 * mu); bp = b / (L * K * kap); A = math.sqrt(L * K * kap) * math.exp(1 / (2 * L)) / (1 - 1 / K)
    Eh = min(E0 / math.sqrt(L * K * kap), 1.0 / 16)
    return A * ((1 + bp / 2) * Eh * math.log(2 * B0(Eh, bp) / Eh) + 0.0325 * bp) / E0 * 1.0    # = (x E1 bound)/E0 ; compare with mu


def Psi2(mu, b, kap, Ks=(2, 4, 8, 16, 32, 64), Ls=(1, 2, 4, 8, 16)):
    """largest P with mu <= rhs2(P) for the best (K, L) (RHS decreasing in P).  Returns P (not capped by 2 mu Qbar^2)."""
    best = S.psi_raw(mu, b)                 # with no admissible (K, L) (L K kappa > 1: tau' > tau) the corollary says nothing beyond CM-1
    for K in Ks:
        for L in Ls:
            if L * K * kap > 1: continue    # the corollary needs tau' <= tau (step (2))
            lo, hi = -700.0, 3.0
            if rhs2(math.exp(hi), mu, b, kap, K, L) >= mu: continue
            for _ in range(120):
                m = 0.5 * (lo + hi)
                if rhs2(math.exp(m), mu, b, kap, K, L) >= mu: lo = m
                else: hi = m
            best = min(best, math.exp(lo))
    return best


if __name__ == '__main__':
    D = Decimal
    # ---- CN-1 steps
    # (1) kernel comparison on |s| <= delta: w_tau^2 <= sqrt(tau'/tau) e^{delta^2/2tau'} w_tau'^2, i.e. e^{-s^2/2tau} <= e^{delta^2/2tau' - s^2/2tau'}
    worst = D(-1)
    for tau, tp, delta in ((D('0.01'), D('0.0004'), D('0.02')), (D(1), D('0.05'), D('0.3')), (D('1e-6'), D('1e-8'), D('5e-5'))):
        for i in range(0, 1001):
            s = delta * i / 1000
            lhs = (-s * s / (2 * tau)).exp(); rhs = (delta * delta / (2 * tp) - s * s / (2 * tp)).exp()
            worst = max(worst, lhs / rhs)
    rec('CN-1 (1) e^{-s^2/2tau} <= e^{(delta^2 - s^2)/2tau\'} on |s| <= delta (60 digits, three cases)', worst <= 1, max_ratio=str(worst)[:12])
    # (2) e^{-s^2/2tau'} <= e^{-s^2/2tau} for tau' <= tau (E0^{(tau')} <= sqrt(tau/tau') E0): trivial sign check
    rec('CN-1 (2) w_tau\'^2 <= sqrt(tau/tau\') w_tau^2 for tau\' <= tau', all((-D(i) ** 2 / (2 * D('0.3'))).exp() <= (-D(i) ** 2 / (2 * D(1))).exp() for i in range(0, 50)))
    # (3) Chebyshev and the whole lemma on random windows (line quadrature)
    rng = np.random.default_rng(5); worst_lem = 0.0; worst_cor = 0.0; rows = []
    for trial in range(40):
        x = float(np.exp(rng.uniform(math.log(1e-5), math.log(3e-3)))); N = int(rng.integers(5, int(min(2000, 4 / math.sqrt(x))))); W = int(rng.integers(1, 150))
        ns = np.concatenate(([1], np.arange(N, N + W))); tc = rng.uniform(0, 2 * math.pi); ph = np.concatenate(([0.0], ns[1:] * tc + rng.normal(0, 0.3, W)))
        ar = S.normalise(ns, np.concatenate(([rng.uniform(0.02, 0.25)], rng.uniform(0.2, 1, W))), ph, x)
        b = float(np.exp(rng.uniform(math.log(1e-3), 0))); tau = x / b; t0 = tc + rng.normal(0, 0.3 * math.sqrt(tau))
        r = S.evaluate(ns, ar, ph, x, t0, b)
        if r is None: continue
        kap = r['kappa']; E1 = r['E1']; E0 = r['E0']
        for K, L in ((4, 2), (16, 4)):
            delta2 = K * kap * tau; tp = L * delta2
            if tp > tau: continue               # the lemma is used with tau' <= tau only
            E0p, E1p, _ = S.quantities(ns, ar, ph, x, t0, tp)
            bound = math.sqrt(tp / tau) * math.exp(delta2 / (2 * tp)) * E1p / (1 - kap * tau / delta2)
            worst_lem = max(worst_lem, E1 / bound)
            worst_cor = max(worst_cor, E0p / (math.sqrt(tau / tp) * E0))
        # the corollary (Psi_2) and the search candidate on the same state
        P = r['P']; mu = r['mu']
        rows.append({'mu': mu, 'kappa': kap, 'P': P, 'Psi2': Psi2(mu, b, kap), 'Psi_sqrt2k': S.psi_raw(mu, b) * min(1, math.sqrt(2 * kap)), 'b': b})
    rec('CN-1 lemma numerically on 40 random states x 2 (K,L): max E1/bound = %.3f <= 1' % worst_lem, worst_lem <= 1 + 1e-9, max_ratio=worst_lem)
    rec('CN-1 (2) numerically: max E0(tau\')/(sqrt(tau/tau\') E0) = %.3f <= 1' % worst_cor, worst_cor <= 1 + 1e-9, max_ratio=worst_cor)
    rec('Psi_2 (proved joint bound) holds on these states: max P/Psi_2 = %.3g' % max(r['P'] / r['Psi2'] for r in rows), all(r['P'] <= r['Psi2'] for r in rows))
    # (4) the size of the proved joint bound against the needed form Psi sqrt(2 kappa): C_2 = Psi_2/(Psi sqrt kappa)
    tab = []
    for mu in (5.8, 8.0, 10.0, 12.0):
        for kap in (0.3, 0.1, 0.01, 1e-3, 1e-4, 1e-5):
            b = 1e-20; p2 = Psi2(mu, b, kap); pr = S.psi_raw(mu, b)
            tab.append({'mu': mu, 'kappa': kap, 'Psi2_over_Psi': p2 / pr, 'C2_Psi2_over_Psi_sqrtkappa': p2 / (pr * math.sqrt(kap))})
    res['C2_table'] = tab
    print('proved joint bound (two foci + CM-1): Psi_2/Psi and C_2 = Psi_2/(Psi sqrt kappa) at b -> 0:', [(t['mu'], t['kappa'], round(t['Psi2_over_Psi'], 3), round(t['C2_Psi2_over_Psi_sqrtkappa'], 1)) for t in tab])
    small = [t for t in tab if t['kappa'] <= 1e-4]; mid = [t for t in tab if t['kappa'] >= 0.01]
    rec('proved joint bound: no gain over CM-1 for kappa >= 0.01 (Psi_2 = Psi); for kappa <= 1e-4 a sqrt(kappa) law with C_2 in [20, 60] (the envelope needs C <= 1.9 at ell = 55.26): the gap is documented, not closed',
        all(abs(t['Psi2_over_Psi'] - 1) < 1e-9 for t in mid) and all(20 <= t['C2_Psi2_over_Psi_sqrtkappa'] <= 60 for t in small))
    # ---- CN-2 steps
    # (5) A[g] = int N_V e^{i alpha v} g dv  equals  E_rho[e^{i beta s/tau} g]  (change of variable s = k' v; beta/tau = alpha/k'): numerically
    for (V, alpha, kp) in ((0.01, 0.3, 2.0), (1e-4, 0.05, 30.0)):
        tau = kp * kp * V; beta = alpha * kp * V
        v = np.linspace(-12 * math.sqrt(V), 12 * math.sqrt(V), 200001); dv = v[1] - v[0]; g = np.cos(3 * kp * v) + 0.5 * np.sin(7 * kp * v)
        lhs = np.sum(np.exp(-v * v / (2 * V)) / math.sqrt(2 * math.pi * V) * np.exp(1j * alpha * v) * g) * dv
        s = kp * v; ds = s[1] - s[0]
        rhs = np.sum(np.exp(-s * s / (2 * tau)) / math.sqrt(2 * math.pi * tau) * np.exp(1j * beta * s / tau) * g) * ds
        rec('CN-2 (i) A[g] = E_rho[e^{i beta s/tau} g] (V=%g)' % V, abs(lhs - rhs) < 1e-9, lhs=complex(lhs).real, rhs=complex(rhs).real)
    # (6) Gaussian average: int Gamma E^2 cos(w/sqrt2 + c) dw = (int Gamma E^2) e^{-se2/4} cos c, se2 = (1/L^2 + 2/s^2)^-1
    for (L2, s2, c) in ((1.0, 0.1, 0.4), (3.0, 0.5, -1.1)):
        se2 = 1 / (1 / L2 + 2 / s2); w = np.linspace(-14 * math.sqrt(L2), 14 * math.sqrt(L2), 400001); dw = w[1] - w[0]
        wt = np.exp(-w * w / (2 * L2)) * np.exp(-w * w / s2)
        lhs = np.sum(wt * np.cos(w / math.sqrt(2) + c)) * dw; rhs = np.sum(wt) * dw * math.exp(-se2 / 4) * math.cos(c)
        rec('CN-2 (ii) tilted Gaussian average of cos(w/sqrt2 + c) (L^2=%g, s^2=%g)' % (L2, s2), abs(lhs - rhs) < 1e-9 * abs(rhs), lhs=lhs, rhs=rhs)
    # (7) 1 - cos z >= z^2/2 - z^4/24 and |cos(a + chi) - cos a| <= |chi|
    zs = np.linspace(-6, 6, 12001); ok7 = bool(np.all(1 - np.cos(zs) >= zs ** 2 / 2 - zs ** 4 / 24 - 1e-12))
    a = rng.uniform(-math.pi, math.pi, 10000); ch = rng.uniform(-1, 1, 10000); ok7 &= bool(np.all(np.abs(np.cos(a + ch) - np.cos(a)) <= np.abs(ch) + 1e-12))
    rec('CN-2 (iii) 1 - cos z >= z^2/2 - z^4/24 and |cos(a+chi) - cos a| <= |chi|', ok7)
    # (8) (s^4/tau^2) e^{-s^2/2tau} <= (64/e^2) e^{-s^2/4tau}  (the s^4 term's kernel bound; constant 64/e^2 = 8.66, times sqrt2 for the 2tau normalisation = 12.25)
    zz = [D(i) / 100 for i in range(0, 6001)]; E_ = D(1).exp()
    m8 = max(4 * z * z * (-z).exp() / (-z / 2).exp() for z in zz)          # z = s^2/2tau: (s^4/tau^2) e^{-s^2/2tau} = 4 z^2 e^{-z}; divided by e^{-z/2}
    rec('CN-2 (iv) max 4 z^2 e^{-z/2} = 64/e^2 (60 digits)', m8 <= 64 / E_ / E_ + D('1e-50'), max_on_grid=str(m8)[:12], bound=str(64 / E_ / E_)[:12])
    json.dump(res, open('results/cn_lemmas.json', 'w'), indent=1, default=float)
    print('ALL OK' if ok_all else 'SOME CHECKS FAILED'); sys.exit(0 if ok_all else 1)
