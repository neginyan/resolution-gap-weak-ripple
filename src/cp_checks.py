"""Stage CP: the search plan's families (F1) point agreement and (F2) box soundness for cp_link.phi_cells_beta, with the wrong versions V1-V4.
Python 3.11.  -> results/cp_checks.json"""
import json, math, sys, time
import numpy as np
import cl_link as K, cp_link as P, cc_cells as CC

def lGc(lmu_a, lmu_b, lr_b):
    lPi = np.where(lmu_a >= math.log(2.0), CC.lPistar_vec(np.maximum(lmu_a, math.log(2.0))), lmu_b + 2 * math.log(0.25)) + 1e-6
    return lPi + 0.5 * (lr_b - np.log(2 + np.exp(lr_b)))

def old_point(lmu, lr, lB, ell, G, Dd=10.0):
    """cl_link on a zero-width box at sigma = e^-ell, b = beta sigma^2 (the reference)."""
    lb = lB - 2 * ell
    return K.phi_cells(lmu, lmu, lr, lr, lb, lb, ell, ell, Dd, G)

def new_box(lmu_a, lmu_b, lr_a, lr_b, lB_a, lB_b, ea, eb, G, Dd=10.0):
    return P.phi_cells_beta(lmu_a, lmu_b, lr_a, lr_b, lB_a, lB_b, ea, eb, Dd, G)

def states(n, rng):
    ell = rng.uniform(55.26, 960.8, n)
    lbmax = math.log(0.003 / 20)
    lBhi = lbmax + 2 * ell
    thick = rng.uniform(0, 1, n) < 0.5
    lB = np.where(thick, rng.uniform(math.log(1e4), math.log(1e6), n), rng.uniform(math.log(1e4), lBhi))
    lmu = rng.uniform(math.log(1e-4), np.log(np.maximum(4 * ell, 50)))
    lr = rng.uniform(-3 * ell, math.log(0.25))
    return ell, lB, lmu, lr

def F1(n=10000, seed=1):
    rng = np.random.default_rng(seed); ell, lB, lmu, lr = states(n, rng); worst = 0.0; sgn = 0
    for i in range(0, n, 500):
        sl = slice(i, i + 500)
        for G in ('0', 'c'):
            for e in np.unique(np.round(ell[sl], 6)):
                pass
        # evaluate per state (ell differs): loop in small groups of equal ell is not possible -> per state
    res = []
    for j in range(n):
        e = float(ell[j]); a = np.array([lmu[j]]); r = np.array([lr[j]]); B = np.array([lB[j]])
        for G in (-np.inf, lGc(a, a, r)):
            vo = float(old_point(a, r, B, e, G)[0]); vn = float(new_box(a, a, r, r, B, B, e, e, G)[0])
            res.append((vo, vn))
    res = np.array(res); fin = np.isfinite(res).all(axis=1)
    d = np.abs(res[fin, 0] - res[fin, 1]); sign_dis = int(np.sum((res[:, 0] > 1e-10) != (res[:, 1] > 1e-10)))
    return {'n_states': n, 'pairs': len(res), 'max_abs_diff_normalised': float(d.max()) if len(d) else None, 'sign_disagreements': sign_dis,
            'both_minus_inf': int(np.sum(~fin)), 'frac_positive': float(np.mean(res[:, 1] > 1e-10))}

def F2(nbox=1000, npts=64, seed=2):
    """random boxes: lower bound > 1e-10 must imply every sampled point positive; record violations and the worst (box bound - min point)."""
    rng = np.random.default_rng(seed); viol = 0; tested = 0; pos_boxes = 0; negpts = 0
    for k in range(nbox):
        ea = rng.uniform(55.26, 950.0); wl = rng.choice([0.0, 0.1, 1.0, 10.0]); eb = ea + wl
        lbmax = math.log(0.003 / 20)
        small = rng.uniform() < 0.6
        lB_a = rng.uniform(math.log(1e4), math.log(1e7)) if small else rng.uniform(math.log(1e4), lbmax + 2 * eb - 1)
        lB_b = lB_a + rng.choice([1e-3, 0.1, 1.0, 5.0])
        lm_a = rng.uniform(math.log(1e-3), math.log(4 * eb)); lm_b = lm_a + rng.choice([1e-3, 0.1, 0.5])
        lr_a = rng.uniform(-3 * eb, math.log(0.25) - 1); lr_b = min(lr_a + rng.choice([1e-3, 0.5, 3.0]), math.log(0.25))
        A = lambda v: np.array([v])
        Dd = float(rng.choice([10.0, 10.0 ** rng.uniform(1, 4)]))     # larger targets make the points of mixed sign (power of the test)
        Gc = lGc(A(lm_a), A(lm_b), A(lr_b))
        for G in (-np.inf, Gc):
            vb = float(new_box(A(lm_a), A(lm_b), A(lr_a), A(lr_b), A(lB_a), A(lB_b), ea, eb, G, Dd)[0])
            if not (vb > 1e-10) or not np.isfinite(vb): continue
            pos_boxes += 1
            for _ in range(npts):
                e = rng.uniform(ea, eb) if wl > 0 else ea
                B = rng.uniform(lB_a, lB_b)
                if B - 2 * e > lbmax: continue                       # not a state of the theorem (b > b_max)
                m = rng.uniform(lm_a, lm_b); r = rng.uniform(lr_a, lr_b)
                Gp = G if not np.isfinite(G) else G               # G_cell is an upper end; the point uses the same G (affine-in-G argument)
                vo = float(old_point(A(m), A(r), A(B), e, Gp, Dd)[0]); tested += 1

                if not (vo > 0): viol += 1
    return {'boxes': nbox, 'positive_boxes': pos_boxes, 'points_tested': tested, 'violations': viol}

def J_components(nbox=2000, npts=50, seed=7):
    """(J2)-(J5) as sampled inequalities.  Boxes: random ell cells (width 0..10), random beta cells; points: admissible states (b <= b_max).
    J2: y_true <= y_max, y_true >= y_min;  J3/J4: sigma <= sigma_eff, b <= b_hi;  J5: for random coefficient sets c_p (p in -1..3, mixed signs),
    exp(P) - exp(N) <= sum_p c_p sigma^p at every sampled sigma.  Returns violation counts."""
    rng = np.random.default_rng(seed); lbmax = math.log(0.003 / K.TAU_MIN); v = {'J2': 0, 'J3J4': 0, 'J5': 0}; pts = 0
    for _ in range(nbox):
        ea = rng.uniform(55.26, 950.0); eb = ea + rng.choice([0.0, 0.5, 10.0])
        top = rng.uniform() < 0.3
        lB_a = rng.uniform(lbmax + 2 * ea - 3, lbmax + 2 * eb) if top else rng.uniform(math.log(1e4), math.log(1e8))
        lB_b = lB_a + rng.choice([1e-3, 0.5, 3.0])
        A = lambda x: np.array([x])
        ls_lo, ls_hi, empty, lb_hi, ly_max, ly_min = P.box_bounds(A(lB_a), A(lB_b), ea, eb)
        # random coefficient sets for J5 (logs in [-5, 5], powers -1..3, each power with a positive and/or a negative term)
        pos = [(A(rng.uniform(-5, 5)), p) for p in (-1, 0, 1, 2, 3) if rng.uniform() < 0.7]
        neg = [(A(rng.uniform(-5, 5)), p) for p in (-1, 0, 1, 2, 3) if rng.uniform() < 0.7]
        lP_, lN_ = P.net_powers(pos, neg, ls_lo, ls_hi, 1); lP0, lN0 = float(lP_[0]), float(lN_[0])
        for _ in range(npts):
            ls = rng.uniform(-eb, -ea) if eb > ea else -ea; lB = rng.uniform(lB_a, lB_b); lb = lB + 2 * ls
            if lb > lbmax: continue
            pts += 1
            if bool(empty[0]): v['J3J4'] += 1; continue
            b = math.exp(lb); ly = 2 * math.log1p(-b) - math.log(2) - lB
            if not (ly <= float(ly_max[0]) + 1e-12 and ly >= float(ly_min[0]) - 1e-12): v['J2'] += 1
            if not (ls <= float(ls_hi[0]) + 1e-12 and ls >= float(ls_lo[0]) - 1e-12 and lb <= float(lb_hi[0]) + 1e-12): v['J3J4'] += 1
            ex = [float(t[0][0]) + t[1] * ls for t in pos + neg] + [x for x in (lP0, lN0) if np.isfinite(x)]; M = max(ex)
            E = lambda x: math.exp(x - M) if np.isfinite(x) else 0.0
            val = sum(E(float(t[0][0]) + t[1] * ls) for t in pos) - sum(E(float(t[0][0]) + t[1] * ls) for t in neg)
            sc = sum(E(float(t[0][0]) + t[1] * ls) for t in pos + neg)
            if E(lP0) - E(lN0) > val + 1e-12 * sc: v['J5'] += 1
    return {'boxes': nbox, 'points': pts, 'violations': v}


def run_wrong_components():
    out = {}
    for vv in ('V1', 'V2', 'V3', 'V4'):
        for kk in P.VARIANT: P.VARIANT[kk] = (kk == vv)
        jc = J_components(600, 20, seed=21); f1 = F1(300, seed=22) if vv == 'V4' else None
        out[vv] = {'components': jc['violations'], 'F1_max_diff': f1 and f1['max_abs_diff_normalised'],
                   'broken': bool(sum(jc['violations'].values()) > 0 or (f1 and f1['max_abs_diff_normalised'] > 0.05))}
    for kk in P.VARIANT: P.VARIANT[kk] = False
    return out


if __name__ == '__main__':
    K.set_tau_min(20.0); t0 = time.time()
    n1 = int(sys.argv[1]) if len(sys.argv) > 1 else 10000; n2 = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    r = {'F1': F1(n1), 'F2': F2(n2), 'J': J_components(2000, 50)}
    print('F1', r['F1'], '\nF2', r['F2'], '\nJ', r['J'], flush=True)
    r['wrong'] = run_wrong_components(); print('wrong', r['wrong'], flush=True); r['seconds'] = round(time.time() - t0)
    json.dump(r, open('results/cp_checks.json', 'w'), indent=1)
