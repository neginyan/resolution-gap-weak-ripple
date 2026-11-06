"""Stage CQ: the search plan's families for the new inequalities (CQ_探索計画.md), with the wrong versions W1-W6.  Python 3.11.
  (F1) identity: the generated functions in the body form equal the v1.0 / CP functions (cl_link.phi_cells, cp_link.phi_cells_beta, cl_cells.cell_ok);
  (F2) components: Pi* decreasing (K1), the skip rule (K5);
  (F3) box soundness: random tail boxes ((i) and (ii)) and wedge boxes; a box accepted at a random target must have every sampled state positive
       (v1.0 point functions at the state's own mu, G, ell, b), with G in {0, G_max(state)} and random fractions;
  (F4) tau factors (K6): 60-digit decimal differences on a log grid of tau in [20, 1e8].
-> results/cq_checks.json"""
import json, math, sys, time
from decimal import Decimal, getcontext
import numpy as np
import cl_link as K, cp_link as P, cl_cells as LC, cc_cells as CC, cq_core as Q, cq_phi as F

A = lambda v: np.array([float(v)])
LQB2 = 2 * math.log(0.25)


def identity(n=3000, seed=1):
    rng = np.random.default_rng(seed); out = {'phi_b_maxdiff': 0.0, 'phi_bl_maxdiff': 0.0, 'sign_dis': 0, 'ok_i_dis': 0, 'boxes': 0}
    for ea in (55.26, 300.26, 950.8):
        eb = ea + 10
        lm_a = rng.uniform(math.log(1e-4), math.log(4 * eb), n); lm_b = lm_a + rng.choice([1e-3, 0.1, 1.0], n)
        lr_a = rng.uniform(-3 * eb, math.log(0.25) - 1, n); lr_b = np.minimum(lr_a + rng.choice([1e-3, 0.5, 3.0], n), math.log(0.25))
        lb_a = rng.uniform(math.log(1e4) - 2 * ea, math.log(0.003 / 20) - 1, n); lb_b = lb_a + rng.choice([1e-3, 0.5, 2.0], n)
        lB_a = rng.uniform(math.log(1e4), math.log(1e4) + 40, n); lB_b = lB_a + rng.choice([1e-3, 0.5, 2.0], n)
        mode = np.zeros(n, int); gt = Q.g_terms(mode, lm_a, lm_b, lr_a, lr_b); lG = gt[0]
        for Dd in (10.0, 300.0):
            a = K.phi_cells(lm_a, lm_b, lr_a, lr_b, lb_a, lb_b, ea, eb, Dd, lG); b = F.phi_b(lm_a, lm_b, lr_a, lr_b, lb_a, lb_b, ea, eb, Dd, gt)
            f = np.isfinite(a) & np.isfinite(b); out['phi_b_maxdiff'] = max(out['phi_b_maxdiff'], float(np.max(np.abs(a[f] - b[f])))); out['sign_dis'] += int(np.sum((a > 1e-10) != (b > 1e-10)))
            a = P.phi_cells_beta(lm_a, lm_b, lr_a, lr_b, lB_a, lB_b, ea, eb, Dd, lG); b = F.phi_bl(lm_a, lm_b, lr_a, lr_b, lB_a, lB_b, np.full(n, ea), np.full(n, eb), Dd, gt)
            f = np.isfinite(a) & np.isfinite(b); out['phi_bl_maxdiff'] = max(out['phi_bl_maxdiff'], float(np.max(np.abs(a[f] - b[f])))); out['sign_dis'] += int(np.sum((a > 1e-10) != (b > 1e-10)))
            out['boxes'] += 2 * n
    for ea in (55.26, 3000.0, 6890.0):
        eb = ea + 10
        lm_a = rng.uniform(math.log(1e-4), math.log(4 * eb), n); lm_b = lm_a + rng.choice([1e-3, 0.1, 1.0], n)
        lr_a = rng.uniform(-3 * eb, math.log(0.25) - 1, n); lr_b = np.minimum(lr_a + rng.choice([1e-3, 0.5, 3.0], n), math.log(0.25))
        lb_a = rng.uniform(-3 * eb, -1, n); lb_b = lb_a + rng.choice([1e-3, 0.5, 2.0], n)
        for dl in (0.002, 0.2, 2.0):
            o = LC.cell_ok(lm_a, lm_b, lr_a, lr_b, lb_a, lb_b, ea, eb, dl)
            nn = Q.ok_i(np.column_stack([lm_a, lm_b, lr_a, lr_b, lb_a, lb_b, np.zeros(n)]), ea, eb, dl); out['ok_i_dis'] += int(np.sum(o != nn)); out['boxes'] += n
    return out


def K1(n=20000, seed=2):
    """Pi*(mu) <= Pi*(mu_a) for mu >= mu_a >= 2, and the upper-tail G used by the code (g_terms, mode UPPER) bounds G_max(mu) for every mu in the box."""
    rng = np.random.default_rng(seed); la = rng.uniform(math.log(2), math.log(4000), n); lm = la + rng.exponential(30, n); lr = rng.uniform(-300, math.log(0.25), n)
    gt = Q.g_terms(np.full(n, Q.UPPER), la, np.full(n, np.inf), lr, lr)
    lGmax = CC.lPistar_vec(lm) + Q.lgs(lr)
    v_G = int(np.sum(lGmax > gt[1] + 1e-12))
    v_ratio = int(np.sum(lGmax - math.log(2) - lm - np.log(2 + np.exp(lr)) > gt[3] + 1e-12))
    return {'states': n, 'violations_G': v_G, 'violations_G_over_mu': v_ratio}


def point_ii(lmu, lr, lb, ell, lG, Dd):
    return float(K.phi_cells(A(lmu), A(lmu), A(lr), A(lr), A(lb), A(lb), ell, ell, Dd, lG)[0])


def lGmax(lmu, lr):
    lP = float(CC.lPistar_vec(A(lmu))[0]) if lmu >= math.log(2) else lmu + LQB2
    return lP + float(Q.lgs(A(lr))[0])


def F3_ii(nbox=300, npts=24, seed=3, which=('lower', 'upper')):
    rng = np.random.default_rng(seed); viol = 0; acc = 0; pts = 0
    for _ in range(nbox):
        ea = rng.uniform(55.26, 950.0); eb = ea + rng.choice([0.0, 1.0, 10.0]); kind = rng.choice(list(which))
        lr_a = rng.uniform(-3 * eb, math.log(0.25) - 1); lr_b = min(lr_a + rng.choice([1e-3, 0.3]), math.log(0.25))
        lb_a = rng.uniform(math.log(1e4) - 2 * ea, math.log(0.003 / 20) - 1); lb_b = lb_a + rng.choice([1e-3, 0.3])
        if kind == 'lower':
            lm_b = rng.uniform(-30.0, 0.0); lm_a = -np.inf if rng.uniform() < 0.5 else lm_b - rng.uniform(0.5, 20); md = Q.LOWER
        else:
            lm_a = math.log(Q.mu_hi_of(eb)) + rng.uniform(0, 3); lm_b = np.inf; md = Q.UPPER
        Dd = float(10.0 ** rng.uniform(1, 3.5))
        X = np.array([[lm_a, lm_b, lr_a, lr_b, lb_a, lb_b, md]])
        if not Q.ok_ii(X, ea, eb, Dd)[0]: continue
        acc += 1
        for _ in range(npts):
            e = rng.uniform(ea, eb) if eb > ea else ea; r = rng.uniform(lr_a, lr_b); b = rng.uniform(lb_a, lb_b)
            if kind == 'lower': m = (rng.uniform(lm_b - 60, lm_b) if not np.isfinite(lm_a) else rng.uniform(lm_a, lm_b))
            else: m = lm_a + rng.exponential(20)
            for u in (0.0, 1.0, rng.uniform()):
                lG = -np.inf if u == 0 else lGmax(m, r) + math.log(u)
                pts += 1
                if not (point_ii(m, r, b, e, lG, Dd) > 0): viol += 1
    return {'boxes': nbox, 'accepted': acc, 'points': pts, 'violations': viol}


def F3_i(nbox=300, npts=24, seed=4):
    rng = np.random.default_rng(seed); viol = 0; acc = 0; pts = 0
    for _ in range(nbox):
        ea = rng.uniform(55.26, 6800.0); eb = ea + rng.choice([0.0, 1.0, 10.0]); kind = rng.choice(['lower', 'upper', 'bzero'])
        lr_a = rng.uniform(-3 * eb, math.log(0.25) - 1); lr_b = min(lr_a + rng.choice([1e-3, 0.3]), math.log(0.25))
        lb_a = rng.uniform(-3 * eb, -1.0); lb_b = lb_a + rng.choice([1e-3, 0.3]); md = Q.BODY
        if kind == 'lower': lm_b = rng.uniform(-30.0, 0.0); lm_a = -np.inf if rng.uniform() < 0.5 else lm_b - rng.uniform(0.5, 20); md = Q.LOWER
        elif kind == 'upper': lm_a = math.log(Q.mu_hi_of(eb)) + rng.uniform(0, 3); lm_b = np.inf; md = Q.UPPER
        else: lm_a = rng.uniform(math.log(1e-4), math.log(4 * eb)); lm_b = lm_a + 0.1; lb_b = -700.0 * 3 - 3 * eb; lb_a = -np.inf
        dl = float(10.0 ** rng.uniform(-3, 0.5))
        if not Q.ok_i(np.array([[lm_a, lm_b, lr_a, lr_b, lb_a, lb_b, md]]), ea, eb, dl)[0]: continue
        acc += 1
        for _ in range(npts):
            e = rng.uniform(ea, eb) if eb > ea else ea; r = rng.uniform(lr_a, lr_b)
            b = rng.uniform(lb_b - 5000, lb_b) if not np.isfinite(lb_a) else rng.uniform(lb_a, lb_b)
            if kind == 'lower': m = (rng.uniform(lm_b - 60, lm_b) if not np.isfinite(lm_a) else rng.uniform(lm_a, lm_b))
            elif kind == 'upper': m = lm_a + rng.exponential(20)
            else: m = rng.uniform(lm_a, lm_b)
            pts += 1
            if not LC.cell_ok(A(m), A(m), A(r), A(r), A(b), A(b), e, e, dl)[0]: viol += 1     # the v1.0 test at the state (G = G_max(mu))
    return {'boxes': nbox, 'accepted': acc, 'points': pts, 'violations': viol}


def F3_w(nbox=300, npts=24, seed=5, big_beta=False):
    rng = np.random.default_rng(seed); viol = 0; acc = 0; pts = 0
    for _ in range(nbox):
        ea = rng.uniform(55.26, 950.0); wl = rng.choice([1e-3, 0.3, 1.0]); eb = ea + wl
        lB_a = rng.uniform(math.log(1e4), math.log(1e4) + 20) if not big_beta else rng.uniform(60, 2 * ea - 20); lB_b = lB_a + rng.choice([1e-3, 0.3])
        lr_a = rng.uniform(-3 * eb, math.log(0.25) - 1); lr_b = min(lr_a + rng.choice([1e-3, 0.3]), math.log(0.25))
        lm_a = rng.uniform(math.log(1e-4), math.log(4 * eb)); lm_b = lm_a + rng.choice([1e-3, 0.3]); md = Q.BODY
        Dd = float(10.0 ** rng.uniform(1, 3.5))
        X = np.array([[lm_a, lm_b, lr_a, lr_b, lB_a, lB_b, ea, eb, md]])
        if not Q.ok_w(X, Dd)[0]: continue
        acc += 1
        for _ in range(npts):
            e = rng.uniform(ea, eb); r = rng.uniform(lr_a, lr_b); B = rng.uniform(lB_a, lB_b); m = rng.uniform(lm_a, lm_b)
            if B - 2 * e > math.log(0.003 / 20): continue
            for u in (0.0, 1.0):
                lG = -np.inf if u == 0 else lGmax(m, r)
                pts += 1
                if not (point_ii(m, r, B - 2 * e, e, lG, Dd) > 0): viol += 1
    return {'boxes': nbox, 'accepted': acc, 'points': pts, 'violations': viol}


def adversarial_w(nbox=120, seed=6):
    """boxes thin in everything but ell (width 1), placed near the envelope of the (ii) inequality (mu ~ ell/4.6, b/sigma^2 ~ e^{ell - 3}, r ~ e^{-ell/2},
    G = G_max), with D tuned so that the state at the box's (ell_a, beta_a, r_b, mu_a) corner is just negative: a sound box test must reject.
    Detects gamma or window constants taken at the wrong ell end (W3, W4).  Returns the number of false acceptances."""
    rng = np.random.default_rng(seed); false_acc = 0; tried = 0
    for _ in range(nbox):
        ea = rng.uniform(55.26, 900.0); eb = ea + 1.0
        lB_a = ea - 3 + rng.uniform(-3, 3); lB_b = lB_a + 1e-6
        lr_a = -ea / 2 - 2 + rng.uniform(-3, 3); lr_b = lr_a + 1e-6; lm_a = math.log(ea / 4.6) + rng.uniform(-0.5, 0.5); lm_b = lm_a + 1e-6
        lG = lGmax(lm_a, lr_b)
        f = lambda D: point_ii(lm_a, lr_b, lB_a - 2 * ea, ea, lG, D)
        lo, hi = -2.0, 2e4
        if not (f(lo) > 0 and f(hi) < 0): continue
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if f(mid) > 0: lo = mid
            else: hi = mid
        Dd = hi + 1e-9 * max(1.0, abs(hi)); tried += 1
        if Q.ok_w(np.array([[lm_a, lm_b, lr_a, lr_b, lB_a, lB_b, ea, eb, Q.BODY]]), Dd)[0]: false_acc += 1
    return {'tried': tried, 'false_acceptances': false_acc}


def phi_i_exact(mu, rho, r, b, sig, gam):
    """the exact reduced function of certificate (i) (cc_core) at G = mu rho, G/mu = rho (finite mu >= 0)."""
    G = mu * rho
    a = r * (2 - gam * sig) + rho * (1 - gam * sig / 2) / (2 + r); c = b / (8 * r) + G / (16 * (2 + r))
    return 2 * math.sqrt(a * c) + 3 * sig / 8 + G * sig / (16 * b) - gam * b - gam * G / 2


def K2_corner(nrect=3000, npts=30, seed=8, only_mu_b=False):
    """(K2) for the exact function of (i): on a rectangle (mu, rho) in [mu_a, mu_b] x [0, rho_cell] (mu_b < 2) the minimum is at a corner.
    Wrong version W1: only the corners with rho = rho_cell (the G = 0 corner forgotten).  (Dropping the mu_a corners instead is harmless here:
    rho = 0 gives G = 0 at every mu, and adding rho only raises the alpha factor, so the G = 0 corner always lies below phi(mu_a, rho).)
    Returns violations (interior value below the corner minimum)."""
    rng = np.random.default_rng(seed); viol = 0
    for _ in range(nrect):
        ell = rng.uniform(55.26, 6900); sig = math.exp(-ell); gam = (1 + rng.uniform(0, 5)) / math.sqrt(ell)
        if rng.uniform() < 0.5:     # generic
            r = math.exp(rng.uniform(-3 * min(ell, 200), math.log(0.25))); b = math.exp(rng.uniform(-2 * min(ell, 300), 0)); mb = math.exp(rng.uniform(-30, math.log(2)))
        else:                       # where phi increases with mu (r ~ 1/4, b small): the corner mu_a matters
            r = math.exp(rng.uniform(-3, math.log(0.25))); b = math.exp(rng.uniform(-14, -2)); mb = math.exp(rng.uniform(-2, math.log(2)))
        rho = 0.0625 * math.sqrt(r / (2 + r)); ma = 0.0 if rng.uniform() < 0.5 else mb * rng.uniform(0, 1)
        qs = (rho,) if only_mu_b else (0.0, rho)
        cmin = min(phi_i_exact(m, q, r, b, sig, gam) for m in (ma, mb) for q in qs)
        scale = abs(cmin) + 3 * sig / 8 + gam * b + 2 * math.sqrt((r * 2 + rho) * (b / (8 * r) + 1)) + 1e-300
        for _ in range(npts):
            m = rng.uniform(ma, mb); q = rng.uniform(0, rho)
            if phi_i_exact(m, q, r, b, sig, gam) < cmin - 1e-12 * scale: viol += 1
    return {'rectangles': nrect, 'points': nrect * npts, 'violations': viol}


def K4(n=20000, seed=9):
    """(K4) components as the generated code takes them (cq_phi._ells): gamma(ell) <= gamma_box and U1, U2, gamma_P, c_chi(ell) <= their box values
    for every ell in the box."""
    rng = np.random.default_rng(seed); ea = rng.uniform(55.26, 950, n); eb = ea + rng.choice([1e-3, 0.3, 1.0, 10.0], n); e = ea + rng.uniform(0, 1, n) * (eb - ea)
    lg, lu = F._ells(ea, eb); Dd = 10.0
    g = lambda l: (1 + (2.1 + Dd) / l + 2.5e-4) / np.sqrt(l)
    vg = int(np.sum(g(e) > g(lg) * (1 + 1e-15)))
    W = F.window_consts_v(e); Wb = F.window_consts_v(lu)
    vu = int(sum(np.sum(a > b_ * (1 + 1e-15)) for a, b_ in zip(W, Wb)))
    return {'boxes': n, 'violations_gamma': vg, 'violations_window_consts': vu}


def K5(nbox=4000, seed=7):
    rng = np.random.default_rng(seed); viol = 0; skipped = 0
    for _ in range(nbox):
        ca = rng.uniform(55.26, 950.0); cb = ca + 10; ea = rng.uniform(ca, cb - 0.01); eb = min(cb, ea + rng.choice([0.01, 1.0, 5.0]))
        lB_a = rng.uniform(math.log(1e4), math.log(1e4) + 2 * (cb - ca)); lB_b = lB_a + rng.choice([0.01, 1.0])
        X = np.array([[0, 1, 0, 1, lB_a, lB_b, ea, eb, 0]])
        if not Q.skip_w(X, ca)[0]: continue
        skipped += 1
        for _ in range(20):
            e = rng.uniform(ea, eb); B = rng.uniform(lB_a, lB_b)
            if B - 2 * e < math.log(1e4) - 2 * ca - 1e-12: viol += 1
    return {'skipped_boxes': skipped, 'violations': viol}


def K6(npts=4000, wrong=False):
    """tau^j nu_k(tau, 1e-6 tau) (j, k as used) and omega_2 decrease on [20, 1e8] (60 digits); wrong version W6: the same with beta = 0.6 tau (c_1 < 0)."""
    getcontext().prec = 60
    t0 = 20.0; bt = Decimal('0.6') if wrong else Decimal('1e-6')
    tmax = Decimal(1e8) if not wrong else Decimal(2000)
    taus = [Decimal(t0) * (tmax / Decimal(t0)) ** (Decimal(i) / Decimal(npts)) for i in range(npts + 1)]
    def nu(k, t, j):
        s = Decimal(0)
        for nn in range(1, 60):
            ex = -Decimal(nn * nn) * t / 2 + Decimal(nn) * bt * t
            if ex < -400 or ex > 400: break
            s += Decimal(nn) ** k * ex.exp()
        return 2 * s * t ** j
    def om2(t):
        s = Decimal(0)
        for nn in range(1, 60):
            ex = -Decimal(nn * nn) * t / 2
            if ex < -400: break
            s += ex.exp() * abs(1 - Decimal(nn * nn) * t)
        return 2 * s
    fs = {'nu0': lambda t: nu(0, t, 0), 'tau nu2': lambda t: nu(2, t, 1), 'tau^2 nu4': lambda t: nu(4, t, 2), 'sqrt(tau) nu1': lambda t: nu(1, t, Decimal('0.5')),
          'tau^1.5 nu3': lambda t: nu(3, t, Decimal('1.5')), 'omega2': om2}
    out = {}
    for name, f in fs.items():
        v = [f(t) for t in taus]; inc = sum(1 for a, b in zip(v, v[1:]) if b > a)
        out[name] = inc
    return {'grid': npts, 'tau_from': t0, 'beta_over_tau': str(bt), 'increases': out}


def run_all(scale=1.0):
    K.set_tau_min(20.0); r = {}; t0 = time.time()
    r['F1_identity'] = identity(int(2000 * scale)); print('F1', r['F1_identity'], flush=True)
    r['K1'] = K1(); r['K5'] = K5(); r['K2'] = K2_corner(int(3000 * scale)); r['K4'] = K4(); print('K1', r['K1'], 'K5', r['K5'], 'K2', r['K2'], 'K4', r['K4'], flush=True)
    r['F3_ii'] = F3_ii(int(300 * scale)); print('F3_ii', r['F3_ii'], flush=True)
    r['F3_i'] = F3_i(int(300 * scale)); print('F3_i', r['F3_i'], flush=True)
    r['F3_w'] = F3_w(int(300 * scale)); print('F3_w', r['F3_w'], flush=True)
    r['adv_w'] = adversarial_w(int(120 * scale)); print('adv_w', r['adv_w'], flush=True)
    r['K6'] = K6(); print('K6', r['K6'], flush=True)
    wr = {}
    wr['W1'] = K2_corner(int(3000 * scale), seed=13, only_mu_b=True)
    Q.WRONG['W2'] = True; wr['W2'] = {'K1': K1(5000, seed=14), 'F3_ii': F3_ii(int(150 * scale), which=('upper',), seed=15)}; Q.WRONG['W2'] = False
    F.CQV['W3'] = True; wr['W3'] = {'K4': K4(5000, seed=16), 'adv': adversarial_w(int(60 * scale), seed=16)}; F.CQV['W3'] = False
    F.CQV['W4'] = True; wr['W4'] = {'K4': K4(5000, seed=17), 'adv': adversarial_w(int(60 * scale), seed=17)}; F.CQV['W4'] = False
    Q.WRONG['W5'] = True; wr['W5'] = K5(4000, seed=18); Q.WRONG['W5'] = False
    wr['W6'] = K6(400, wrong=True)
    r['wrong'] = wr; print('wrong', wr, flush=True)
    br = {'W1': wr['W1']['violations'] > 0, 'W2': wr['W2']['K1']['violations_G'] + wr['W2']['K1']['violations_G_over_mu'] + wr['W2']['F3_ii']['violations'] > 0,
          'W3': wr['W3']['K4']['violations_gamma'] + wr['W3']['adv']['false_acceptances'] > 0, 'W4': wr['W4']['K4']['violations_window_consts'] + wr['W4']['adv']['false_acceptances'] > 0, 'W5': wr['W5']['violations'] > 0, 'W6': sum(wr['W6']['increases'].values()) > 0}
    r['wrong_broken'] = br; r['seconds'] = round(time.time() - t0)
    return r


if __name__ == '__main__':
    sc = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    r = run_all(sc); print('broken', r['wrong_broken'], flush=True); json.dump(r, open('results/cq_checks.json', 'w'), indent=1)
