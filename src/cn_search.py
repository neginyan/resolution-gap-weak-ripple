"""Stage CN 2.1: try to break the candidate joint lemma ("a bright local ripple is spread over the focus") before proving anything.
Plan: CN_探索計画.md.  Python 3.11.

Local quantities at focus t0, variance tau = x/b (exact Fourier formulas; wrapped normal has characteristic function e^{-k^2 tau/2}):
  E0 = <theta R^2>, E1 = <theta R'^2>, E1ring = <(theta + 2 tau d theta/d tau) R'^2> = E1 + 2 tau dE1/dtau  (kernel = sum_j (s-2pi j)^2/tau * Gaussian >= 0)
  mu = x E1/E0, P = 2 x E1, kappa = E1ring/E1.
Candidates and deliberately wrong versions: see CANDS below.  Families: coherent packets at the positivity limit (node/crest/trough of the
fundamental; focus at/off the packet; widths narrower and wider than sqrt tau; cos and sin phases), two packets, combs, random bands, hill climb.
-> results/cn_search.json
"""
import json, math, os, sys, time
import numpy as np
import cm_local_bound as LB

QB = 0.25
OUT = 'results/cn_search.json'
NG = 1 << 16


def psi_raw(mu, b):
    """the CM-1 bound Psi(mu, b) WITHOUT the cap 2 mu Qbar^2 (largest P with mu <= rhs(P); rhs decreasing in P)."""
    lo, hi = -700.0, 3.0
    if LB.rhs(math.exp(hi), mu, b) >= mu: return math.exp(hi)
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if LB.rhs(math.exp(m), mu, b) >= mu: lo = m
        else: hi = m
    return math.exp(lo)


def psi_min(mu, b):
    return min(2 * mu * QB * QB, psi_raw(mu, b))


# ---- candidates: ratio = LHS/RHS; > 1 breaks the statement -------------------------------------------------------------------
def c1(mu, P, kap, b):   # C1': the joint lemma; below mu = 2 the weak bound 2 mu Qbar^2 (no kappa factor: it is a bound on E0 alone)
    return P / min(2 * mu * QB * QB, psi_raw(mu, b) * min(1.0, math.sqrt(2 * kap)))
def c2(mu, P, kap, b): return P / (1.2 * mu * math.exp(-mu) * min(1.0, math.sqrt(2 * kap))) if mu >= 2 else 0.0
def w1(mu, P, kap, b): return 0.5 / kap                                     # "kappa >= 1/2 always"
def w2(mu, P, kap, b): return P / min(2 * mu * QB * QB, psi_raw(mu, b) * min(1.0, 2 * kap))      # kappa instead of sqrt kappa
def w3(mu, P, kap, b): return P / (0.3 * mu * math.exp(-mu) * min(1.0, math.sqrt(2 * kap))) if mu >= 2 else 0.0
def w4(mu, P, kap, b): return P / min(2 * mu * QB * QB, psi_raw(mu, b) * min(1.0, math.sqrt(2 * kap)) / 10)
CANDS = {'C1': c1, 'C2': c2, 'W1': w1, 'W2': w2, 'W3': w3, 'W4': w4}


# ---- windows -------------------------------------------------------------------------------------------------------------------
def f_min(ns, araw, phi, M=None):
    """min of f = 1 + sum araw cos(n t - phi) on M grid points, minus the midpoint margin (h/2)^2 sum n^2 a_n / 2 (FFT synthesis)."""
    if M is None:
        M = NG
        while M < 16 * int(np.max(ns)): M *= 2                 # grid fine enough for the highest harmonic (16 points per period)
    c = np.zeros(M // 2 + 1, complex); c[ns] = araw * np.exp(-1j * phi) / 2
    f = 1.0 + np.fft.irfft(c * M, M)                        # irfft(X, M) = (1/M)[X_0 + 2 Re sum_{n>=1} X_n e^{int}] -> sum a_n cos(n t - phi_n)
    h = 2 * math.pi / M
    return float(f.min()) - (h / 2) ** 2 * float(np.sum(ns.astype(float) ** 2 * araw)) / 2


def normalise(ns, araw, phi, x):
    """enforce (P2) (fundamental largest), f >= 0 and (P3) by scaling the ripple (keeps f >= 0 and (P2)); returns araw."""
    araw = np.asarray(araw, float).copy()
    if ns[0] != 1: raise ValueError('first harmonic must be the fundamental')
    araw[0] = max(araw[0], araw.max())
    m = f_min(ns, araw, phi)
    if m < 0: araw *= 0.999 / (1 - m)
    s = float(np.sum(araw * np.exp(-ns ** 2 * x / 2)))
    if s > QB: araw *= 0.999 * QB / s
    return araw


def quantities(ns, araw, phi, x, t0, tau):
    q = araw * np.exp(-ns ** 2 * x / 2); al = ns * t0 - phi
    D = ns[:, None] - ns[None, :]; S = ns[:, None] + ns[None, :]
    cd = np.cos(al[:, None] - al[None, :]); cs = np.cos(al[:, None] + al[None, :])
    qq = q[:, None] * q[None, :]; nm = ns[:, None] * ns[None, :]
    ed = np.exp(-D ** 2 * tau / 2); es = np.exp(-S ** 2 * tau / 2)
    E0 = 0.5 * float(np.sum(qq * (ed * cd + es * cs)))
    E1 = 0.5 * float(np.sum(nm * qq * (ed * cd - es * cs)))
    E1r = 0.5 * float(np.sum(nm * qq * ((1 - D ** 2 * tau) * ed * cd - (1 - S ** 2 * tau) * es * cs)))
    return E0, E1, E1r


def evaluate(ns, araw, phi, x, t0, b):
    tau = x / b
    E0, E1, E1r = quantities(ns, araw, phi, x, t0, tau)
    if E0 <= 0 or E1 <= 0 or E1r <= 0: return None
    mu = x * E1 / E0; P = 2 * x * E1; kap = E1r / E1
    return {'mu': mu, 'P': P, 'kappa': kap, 'E0': E0, 'E1': E1, 'E1ring': E1r, 'x': x, 'b': b, 't0': t0, 'ratios': {k: f(mu, P, kap, b) for k, f in CANDS.items()}}


class Best:
    def __init__(self): self.b = {k: (0.0, None) for k in CANDS}
    def add(self, r, desc):
        if r is None: return
        for k, v in r['ratios'].items():
            if v > self.b[k][0]: self.b[k] = (v, {**{kk: vv for kk, vv in r.items() if kk != 'ratios'}, 'desc': desc, 'ratio': v})


# ---- families ------------------------------------------------------------------------------------------------------------------
def packet(a1, N1, N2, tc, x, phase='cos'):
    ns = np.concatenate(([1], np.arange(N1, N2 + 1))); W = N2 - N1 + 1
    ph = np.concatenate(([0.0], ns[1:] * tc + (0.0 if phase == 'cos' else math.pi / 2)))
    lo, hi = 0.0, a1
    for _ in range(22):
        a = 0.5 * (lo + hi); ar = np.concatenate(([a1], np.full(W, a)))
        if f_min(ns, ar, ph) >= 0 and np.sum(ar * np.exp(-ns ** 2 * x / 2)) <= QB: lo = a
        else: hi = a
    return ns, np.concatenate(([a1], np.full(W, lo * (1 - 1e-9)))), ph


def fam_packets(best, log):
    n = 0
    for x in (3e-3, 1e-3, 1e-4, 1e-5, 1e-6):
        for mu0 in (3.0, 5.0, 8.0, 12.0):
            N = int(round(math.sqrt(mu0 / x)))
            for W in (1, 2, 4, 8, 16, 32, 64, 128, 256):
                if N - W // 2 < 2 or W > N: continue
                for tc_name, tc in (('node', math.pi / 2), ('crest', 0.0), ('trough', math.pi)):
                    for phase in ('cos', 'sin'):
                        for a1 in (0.05, 0.2):
                            ns, ar, ph = packet(a1, N - W // 2, N - W // 2 + W - 1, tc, x, phase)
                            if ar[1] <= 0: continue
                            for b in (1.0, 0.3, 0.1, 0.01):
                                tau = x / b
                                for dt_name, dt in (('centre', 0.0), ('+sqrt', math.sqrt(tau)), ('+2sqrt', 2 * math.sqrt(tau)), ('far node', math.pi / 2 - tc if tc != math.pi / 2 else math.pi)):
                                    r = evaluate(ns, ar, ph, x, tc + dt, b); n += 1
                                    best.add(r, 'packet x=%g mu0=%g W=%d at %s %s a1=%g focus %s b=%g' % (x, mu0, W, tc_name, phase, a1, dt_name, b))
    log('packets: %d evaluations' % n)


def fam_narrow(best, log):
    """packets narrower than the focus at a small blur x (so that e^{-n^2 x/2} is flat across a band of thousands of harmonics): the family that
    decides between sqrt(kappa) and kappa at mu_theta >= 2."""
    n = 0; starts = []
    for x in (1e-8, 1e-7):
        for mu0 in (5.0, 8.0, 12.0):
            N = int(round(math.sqrt(mu0 / x)))
            for W in (250, 500, 1000):
                for tc_name, tc in (('node', math.pi / 2), ('crest', 0.0)):
                    for phase in ('cos', 'sin'):
                        ns, ar, ph = packet(0.2, N - W // 2, N - W // 2 + W - 1, tc, x, phase)
                        if ar[1] <= 0: continue
                        for b in (1e-6, 1e-5, 1e-4, 1e-3):
                            tau = x / b
                            for dt_name, dt in (('centre', 0.0), ('+sqrt', math.sqrt(tau))):
                                r = evaluate(ns, ar, ph, x, tc + dt, b); n += 1
                                best.add(r, 'narrow packet x=%g mu0=%g W=%d at %s %s focus %s b=%g' % (x, mu0, W, tc_name, phase, dt_name, b))
                                if r and r['mu'] >= 2 and dt == 0.0 and b == 1e-5 and W <= 500: starts.append((ns, ar, ph, x, tc, b))
    log('narrow packets (small x): %d evaluations' % n)
    return starts


def gpacket(a1, N, W, tc, x, sig_n):
    """packet with a Gaussian envelope in n (compact envelope in t): the family that separates the sqrt(kappa) law from the kappa law."""
    ns = np.concatenate(([1], np.arange(N - W // 2, N - W // 2 + W))); env = np.exp(-((ns[1:] - N) / sig_n) ** 2 / 2)
    ph = np.concatenate(([0.0], ns[1:] * tc))
    lo, hi = 0.0, 10.0
    for _ in range(26):
        a = 0.5 * (lo + hi); ar = np.concatenate(([a1], a * env))
        if ar[1:].max() <= a1 and f_min(ns, ar, ph) >= 0 and np.sum(ar * np.exp(-ns ** 2 * x / 2)) <= QB: lo = a
        else: hi = a
    return ns, np.concatenate(([a1], lo * env * (1 - 1e-9))), ph


def fam_gauss_narrow(best, log):
    n = 0; starts = []
    for x in (1e-7, 1e-8):
        for mu0 in (5.0, 8.0, 12.0):
            N = int(round(math.sqrt(mu0 / x)))
            for W in (1000, 2000, 4000):
                sig_n = W / 6.7
                for tc_name, tc in (('node', math.pi / 2), ('crest', 0.0)):
                    ns, ar, ph = gpacket(1.0 / sig_n, N, W, tc, x, sig_n)
                    if ar[1:].max() <= 0: continue
                    for b in (1e-6, 1e-7, 3e-8):
                        tau = x / b
                        for dt_name, dt in (('centre', 0.0), ('+sqrt', math.sqrt(tau))):
                            r = evaluate(ns, ar, ph, x, tc + dt, b); n += 1
                            best.add(r, 'gauss packet x=%g mu0=%g W=%d sig_n=%g at %s cos focus %s b=%g' % (x, mu0, W, sig_n, tc_name, dt_name, b))
                            if r and r['mu'] >= 2 and dt == 0.0 and b == 1e-6 and W <= 1000: starts.append((ns, ar, ph, x, tc, b))
    log('gaussian-envelope narrow packets: %d evaluations' % n)
    return starts


def fam_two_packets(best, log):
    n = 0
    for x in (1e-4, 1e-6):
        for mu0 in (5.0, 10.0):
            N = int(round(math.sqrt(mu0 / x)))
            for Wn, Ww in ((64, 4), (128, 2), (256, 8)):
                for b in (1.0, 0.1):
                    tau = x / b
                    if Wn < 1 / math.sqrt(tau): continue
                    ns = np.concatenate(([1], np.arange(N - Wn // 2, N - Wn // 2 + Wn), np.arange(2 * N, 2 * N + Ww)))
                    tc = math.pi / 2
                    for sh in (0.5, 1.0, 2.0):
                        ph = np.concatenate(([0.0], ns[1:1 + Wn] * tc, ns[1 + Wn:] * tc))
                        ar = np.concatenate(([0.1], np.full(Wn, 1.0 / Wn * sh), np.full(Ww, 1.0 / Ww)))
                        ar = normalise(ns, ar, ph, x)
                        r = evaluate(ns, ar, ph, x, tc, b); n += 1
                        best.add(r, 'two packets x=%g mu0=%g Wn=%d Ww=%d sh=%g b=%g' % (x, mu0, Wn, Ww, sh, b))
    log('two packets: %d evaluations' % n)


def fam_combs(best, log):
    n = 0
    for x in (3e-3, 1e-4, 1e-6):
        for w in (0.1 * x, x, 10 * x, 100 * x):
            M = int(min(3000, 6 / math.sqrt(w)))
            ns = np.arange(1, M + 1); ar = 2 * np.exp(-ns ** 2 * w / 2); ph = np.zeros(M)
            ar = normalise(ns, ar, ph, x)
            for b in (1.0, 0.1, 0.01):
                tau = x / b
                for t0_name, t0 in (('tooth', 0.0), ('between', math.pi), ('sqrt tau off', math.sqrt(tau))):
                    r = evaluate(ns, ar, ph, x, t0, b); n += 1
                    best.add(r, 'comb x=%g w/x=%g M=%d focus %s b=%g' % (x, w / x, M, t0_name, b))
    log('combs: %d evaluations' % n)


def fam_random(best, log, rng, n_win=3000):
    n = 0
    for i in range(n_win):
        x = float(np.exp(rng.uniform(math.log(1e-6), math.log(3e-3))))
        N = int(rng.integers(5, int(min(3000, 6 / math.sqrt(x))))); W = int(rng.integers(1, 200))
        ns = np.concatenate(([1], np.arange(N, N + W)))
        kind = rng.integers(3)
        if kind == 0: ph = rng.uniform(0, 2 * math.pi, len(ns)); ar = np.concatenate(([rng.uniform(0.02, 0.25)], rng.uniform(0, 1, W)))
        elif kind == 1: tc = rng.uniform(0, 2 * math.pi); ph = np.concatenate(([0.0], ns[1:] * tc)); ar = np.concatenate(([rng.uniform(0.02, 0.25)], rng.uniform(0.5, 1, W)))
        else: tc = rng.uniform(0, 2 * math.pi); ph = np.concatenate(([0.0], ns[1:] * tc + rng.normal(0, 0.3, W))); ar = np.concatenate(([rng.uniform(0.02, 0.25)], np.exp(-((np.arange(W) - W / 2) / (W / 4 + 1)) ** 2)))
        ar = normalise(ns, ar, ph, x)
        b = float(np.exp(rng.uniform(math.log(1e-3), 0)))
        t0 = rng.uniform(0, 2 * math.pi) if kind == 0 else (tc + rng.normal(0, math.sqrt(x / b)) if rng.random() < 0.7 else tc)
        r = evaluate(ns, ar, ph, x, t0, b); n += 1
        best.add(r, 'random kind=%d x=%g N=%d W=%d b=%g' % (kind, x, N, W, b))
    log('random: %d evaluations' % n)


def hill(best, log, rng, cand, start, steps=200):
    """perturb raw amplitudes (log), phases, focus position and log tau from `start`, keep the class, maximise the ratio of `cand`."""
    ns, ar, ph, x, t0, b = start
    cur = evaluate(ns, ar, ph, x, t0, b); curv = cur['ratios'][cand] if cur else 0.0
    for s in range(steps):
        sc = 0.3 * (1 - s / steps) + 0.02
        ar2 = ar * np.exp(rng.normal(0, sc, len(ar)) * (rng.random(len(ar)) < 0.3)); ph2 = ph + rng.normal(0, sc, len(ph)) * (rng.random(len(ph)) < 0.3)
        t02 = t0 + rng.normal(0, sc * math.sqrt(x / b)); b2 = float(min(1.0, max(1e-3, b * math.exp(rng.normal(0, sc)))))
        ar2 = normalise(ns, ar2, ph2, x)
        r = evaluate(ns, ar2, ph2, x, t02, b2)
        if r and r['ratios'][cand] > curv: ar, ph, t0, b, cur, curv = ar2, ph2, t02, b2, r, r['ratios'][cand]
    best.add(cur, 'hill climb for %s from %s' % (cand, 'start'))
    return curv, cur


def main():
    t_start = time.time(); logs = []
    def log(m): logs.append(m); print(m, flush=True)
    rng = np.random.default_rng(20261105)
    best = Best()
    fam_packets(best, log); narrow_starts = fam_narrow(best, log); gauss_starts = fam_gauss_narrow(best, log); fam_two_packets(best, log); fam_combs(best, log); fam_random(best, log, rng)
    # hill climbs: from the worst case of each candidate (if it is a packet/random we need its window; rebuild some starts)
    starts = []
    for x, mu0, W, tc, phase, b in ((1e-6, 12.0, 256, math.pi / 2, 'sin', 1.0), (1e-5, 8.0, 64, math.pi / 2, 'cos', 1.0), (1e-4, 5.0, 16, math.pi / 2, 'sin', 0.1),
                                    (3e-3, 5.0, 8, math.pi / 2, 'cos', 1.0), (1e-6, 5.0, 128, 0.0, 'sin', 1.0), (1e-5, 12.0, 32, math.pi, 'cos', 0.3)):
        N = int(round(math.sqrt(mu0 / x))); ns, ar, ph = packet(0.1, N - W // 2, N - W // 2 + W - 1, tc, x, phase); starts.append((ns, ar, ph, x, tc, b))
    for i in range(14):
        x = float(np.exp(rng.uniform(math.log(1e-6), math.log(3e-3)))); N = int(rng.integers(5, int(min(3000, 6 / math.sqrt(x))))); W = int(rng.integers(1, 200))
        ns = np.concatenate(([1], np.arange(N, N + W))); tc = rng.uniform(0, 2 * math.pi); ph = np.concatenate(([0.0], ns[1:] * tc + rng.normal(0, 0.2, W)))
        ar = normalise(ns, np.concatenate(([0.1], rng.uniform(0.3, 1, W))), ph, x); starts.append((ns, ar, ph, x, tc, float(np.exp(rng.uniform(math.log(1e-2), 0)))))
    starts += narrow_starts[:4] + gauss_starts[:3]
    hc = {}
    for cand in CANDS:
        vals = [hill(best, log, rng, cand, st)[0] for st in starts]
        hc[cand] = max(vals); log('hill climb %s: max ratio %.4g over %d starts' % (cand, max(vals), len(starts)))
    out = {'plan': 'CN_探索計画.md', 'seconds': time.time() - t_start, 'log': logs,
           'worst': {k: {'ratio': v[0], 'case': v[1]} for k, v in best.b.items()}, 'hill_climb_max': hc,
           'verdict': {k: ('BROKEN' if best.b[k][0] > 1 else 'survived') for k in CANDS}}
    for k in CANDS: log('%s: worst ratio %.4g -> %s  (%s)' % (k, best.b[k][0], out['verdict'][k], (best.b[k][1] or {}).get('desc', '')))
    os.makedirs('results', exist_ok=True); json.dump(out, open(OUT, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
