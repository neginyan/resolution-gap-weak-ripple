"""Stage CP 2.1: certificate (ii) of paper H (cl_link.py) rebuilt with boxes in (log mu, log r, log beta), beta = b/sigma^2.  Python 3.11.

WHY.  cl_link.certify boxes b in [b_min, b_max] with b_min = 1e4 sigma_top^2 (sigma_top = e^{-ell_a}), while Theorem H1 and Lemma BE split the states
at b = 1e4 sigma^2 with the state's own sigma: the wedge 1e4 sigma^2 <= b < 1e4 sigma_top^2 inside each column was covered by neither (stage CO,
sec. 1b).  With beta as the box variable the boundary beta = 1e4 is a face of the boxes and the wedge disappears.

THE SAME INEQUALITY.  phi_cells_beta bounds exactly the expression of cl_link.phi_cells (Num_R* - gamma sigma den_R*, divided by sigma, in the
variable a' = 4 r A/sigma), term by term; only the bookkeeping of b changes:
  (J1) every term c b^k sigma^p becomes c beta^k sigma^(p + 2k)  (identity);
  (J2) y = (1 - b)^2 sigma^2/(2 b) = (1 - b)^2/(2 beta) <= 1/(2 beta_a), y >= (1 - b_hi)^2/(2 beta_b)   (sigma-free);
  (J3) the theorem's states have b <= b_max = 0.003/tau_min, so a box with beta >= beta_a only contains sigma <= sigma_eff = min(sigma_top, sqrt(b_max/beta_a));
       a box with sigma_eff < sigma_bot is empty (accepted, recorded);
  (J4) b <= b_hi = min(beta_b sigma_eff^2, b_max);
  (J5) terms of equal power p are netted at equal sigma; a positive net is then bounded by min(sigma_bot^p, sigma_eff^p), a negative net by
       max(sigma_bot^p, sigma_eff^p) (for p < 0 the ends swap).  Powers now range over -1..3.
Everything else (the term list, the H bracket, the a' range, the B1/B2 branches, keepF, G in {0, G_cell}) is cl_link's.  The H bracket, formed at
equal sigma in cl_link, is here netted power by power (b-part and sigma^2 gamma_P part are both power 2 in beta).

ROUNDING (the rule of cl_link, rederived).  Each term log is a sum of k <= 17 operands (one fewer than in cl_link for the y-terms, since
2 log sigma - log b became -log beta); the operands are bounded by B = 2(max|lmu| + max|lr| + max|lbeta|) + 3 ell_b + 600 as before (|p log sigma| <= 3 ell_b
for |p| <= 3).  log sigma_eff is computed with 3 operations and then moved outwards by 1e-9 (|error| <= 4u * 1000 < 1e-12), so sigma_eff is an upper
bound of the true one.  The perturbation eta_t = S(25u(|l_t| + 3B + y_c B) + 1000u), S = 10, is therefore unchanged.  Acceptance: normalised minimum > 1e-10.

CLI: python cp_link.py ELL_FROM ELL_TO STEP [D] [ln beta_min] [tau_min]   -> $CP_LINK_FN (default results/cp_link_cells.jsonl), resumable.
"""
import json, math, os, sys, time
import numpy as np
import cc_cells as CC
import cl_link as K

FN = os.environ.get('CP_LINK_FN', 'results/cp_link_cells.jsonl')
U = K.U; SAFE = K.SAFE; QB = K.QB; ETAB = K.ETAB; LN2 = K.LN2; CZETA = K.CZETA
HSHIFT_LO = True     # (J6); False reproduces the first trial (H shift sigma_top only)
VARIANT = {'V1': False, 'V2': False, 'V3': False, 'V4': False}    # deliberately wrong versions (search plan sec. 3); all False in the certificate


def lse(*a): return K.lse(*a)


def box_bounds(lB_a, lB_b, ell_a, ell_b):
    """(J2)-(J4): log sigma ends (bottom, effective top), emptiness, log b_hi, log y_max, log y_min for beta boxes."""
    n = len(lB_a); lbmax = math.log(0.003 / K.TAU_MIN)
    ls_lo = np.full(n, -float(ell_b))
    ls_hi = np.minimum(-ell_a, 0.5 * (lbmax - lB_a)) + 1e-9                     # (J3), moved outwards
    empty = ls_hi < ls_lo - 1e-9
    ls_hi = np.maximum(ls_hi, ls_lo)
    lb_hi = np.minimum(lB_b + 2 * (ls_lo if VARIANT['V3'] else ls_hi), lbmax)   # (J4)
    ly_max = -LN2 - (lB_b if VARIANT['V2'] else lB_a)                           # (J2)
    ly_min = 2 * np.log1p(-np.exp(lb_hi)) - LN2 - lB_b
    return ls_lo, ls_hi, empty, lb_hi, ly_max, ly_min


def net_powers(pos, neg, ls_lo, ls_hi, n):
    """(J5): terms (log coefficient, power p) netted at equal sigma per p; positive nets at min(sigma_bot^p, sigma_top^p), negative at the max.
    Returns (log total positive, log total negative)."""
    NEG = -np.inf
    def worst_pos(p): return np.minimum(p * ls_lo, p * ls_hi)
    def worst_neg(p):
        if VARIANT['V1'] and p < 0: return p * ls_hi
        return np.maximum(p * ls_lo, p * ls_hi)
    P, N = [], []
    for p in sorted(set([t[1] for t in pos] + [t[1] for t in neg])):
        lp = lse(*[t[0] for t in pos if t[1] == p]) if any(t[1] == p for t in pos) else np.full(n, NEG)
        ln = lse(*[t[0] for t in neg if t[1] == p]) if any(t[1] == p for t in neg) else np.full(n, NEG)
        d = np.where(lp > ln, lp + np.log1p(-np.exp(np.minimum(ln - lp, 0.0))), NEG)
        e_ = np.where(ln >= lp, ln + np.log1p(-np.exp(np.minimum(lp - ln, 0.0))), NEG)
        e_ = np.where((ln >= lp) & ~np.isfinite(lp), ln, e_)
        P.append(d + worst_pos(p)); N.append(e_ + worst_neg(p))
    if not P: return np.full(n, NEG), np.full(n, NEG)
    return lse(*P), lse(*N)


def phi_cells_beta(lmu_a, lmu_b, lr_a, lr_b, lB_a, lB_b, ell_a, ell_b, Dd, lG):
    n = len(lmu_a); NEG = -np.inf
    NU, Z, CZ, EPS1, FAC = K.NU, K.Z, K.CZ, K.EPS1, K.FAC
    G0 = np.isscalar(lG) and not np.isfinite(lG)
    ls_lo, ls_hi, empty, lb_hi, ly_max, ly_min = box_bounds(lB_a, lB_b, ell_a, ell_b)
    gam = (1 + (2.1 + Dd) / ell_a + 2.5e-4) / math.sqrt(ell_a); lgam = math.log(gam)
    U1, U2, gP, cchi = K.window_consts(ell_b); X1 = U1 / (1 - QB); lX1 = math.log(X1)
    r_lo = np.exp(lr_a); r_hi = np.exp(lr_b)
    b_hi = np.exp(lb_hi)
    y_max = np.exp(np.minimum(ly_max, 700)); y_min = np.exp(ly_min)
    e_min = np.exp(-y_max / 2)
    kF_lo = e_min * (1 - NU['nu0']) / (1 + Z)
    w0_max = 1.5 * y_max
    if G0:
        lg0_lo = lg0_hi = np.full(n, NEG); lGv = NEG
    else:
        lGv = lG; lg0_lo = lG - LN2 - lmu_b - np.log(2 + r_hi); lg0_hi = lG - LN2 - lmu_a - np.log(2 + r_lo)
    z0 = np.zeros(n)
    sh2 = 0 if VARIANT['V4'] else 2                                             # the power shift of (J1) per power of b
    Lp, Ln, Cp, Cn, Sn = [], [], [], [], []
    kv = (1 + NU['om2']) * CZ / 4; kw = 1.5 * (1 + NU['nu0']) * CZ; kq = CZ * QB * (NU['nu0'] + NU['tnu2'])
    # ---- L (cl_link lines: Lp, Ln)
    Lp += [(z0 + math.log((1 - EPS1) / 2), 0), (math.log(1 - EPS1) + lg0_lo - LN2 - lr_b + math.log(FAC['low_S0']), 0)]
    Ln += [(math.log(kv / 4) - lB_a, 2 - sh2), (math.log(kv / 4) - lB_a + lg0_hi - lr_a, 2 - sh2),
           (z0 + math.log(kw * (X1 ** 2 + CZETA) / 4), 2), (math.log(kw * (X1 ** 2 + CZETA) / 4) + lg0_hi - lr_a, 2),
           (math.log(3 * kw / 4) + lr_b, 2), (math.log(3 * kw / 4) + lr_b + lg0_hi - np.log(2 + r_lo), 2),
           (z0 + lgam + math.log((1 + EPS1) / 4), 1), (lgam + math.log(2 * (1 + EPS1) / 4) + lg0_hi - lr_a, 1)]
    LnB1 = []
    # ---- H bracket, power by power: b/(2(1+z)(1+w0)^2) [p 2] + kF G r/(4(2+r)) [p 0] - kq b/2 [p 2] - sigma^2 gamma_P cz [p 2]
    Hp = [(lB_a - 2 * np.log1p(w0_max) - math.log(2 * (1 + Z)), sh2)]
    HpG = [] if G0 else [(np.log(kF_lo) + lGv + lr_a - np.log(2 + r_lo) - math.log(4) + math.log(FAC['low_S2']), 0)]
    Hn = [(lB_b + math.log(kq / 2), sh2), (z0 + math.log(gP * CZ), 2)]
    # ---- C
    Cp += [(math.log(3 / (8 * (1 + Z))) + 2 * np.log1p(-b_hi) - np.log1p(w0_max), 1)]
    Cn += [(z0 + math.log(kw), 1), (math.log(kw) + lg0_hi - np.log(2 + r_lo), 1)]
    has = (1 - y_min) > 0
    Cn += [(np.where(has, lB_b + np.log(np.maximum(1 - y_min, 1e-300)), NEG) + lX1 - math.log(2 * math.sqrt(2) * (1 + Z)), sh2)]
    Cn += [(z0 + math.log(gP * CZ * X1 / (2 * math.sqrt(2))), 2)]
    Cn += [(lB_b + math.log(kq * X1 / (4 * math.sqrt(2))), sh2)]
    sb = 0.5 * lB_b; sh1 = sh2 / 2
    Cn += [(z0 + math.log(CZ / 2 * 3 / 8 * QB * NU['nu0']), 1),
           (lB_b + math.log(CZ / 2 * (1 + QB) * cchi * NU['nu0'] * (1 + QB / 2)), sh2),
           (z0 + math.log(CZ / 2 * QB * NU['tnu2'] / 8), 1),
           (z0 + math.log(CZ / 2 * NU['nu0'] * (1 + QB) * gP), 1),
           (sb + math.log(CZ / 2 * NU['stnu1'] * (U1 * cchi + (1 + QB) * gP / math.sqrt(2))), sh1),
           (lB_b + math.log(CZ / 2 * 0.5 * QB * (1 + QB) * cchi * NU['tnu2']), sh2),
           (z0 + math.log(CZ / 2 * QB * NU['t2nu4'] / 8), 1),
           (lB_b + math.log(CZ / 2 * QB * NU['tnu2'] / (2 * (1 - NU['nu0']))), 1 + sh2),
           (sb + math.log(CZ / 2 * QB * NU['t15nu3'] / (2 * math.sqrt(2))), sh1)]
    Cn += [(z0 + lgam + lB_b + math.log(1 + NU['tnu2'] * QB / (2 * (1 - Z))), sh2)]
    # ---- S
    lR1 = np.log(math.sqrt(CZETA) + np.sqrt(2 * r_hi / math.pi)) - 0.5 * LN2
    Sn += [(np.where(has, lB_b + np.log(np.maximum(1 - y_min, 1e-300)), NEG) + lR1 - math.log(2 * (1 + Z)) - LN2 - 0.5 * lr_a, 0.5 + sh2),
           (math.log(gP * CZ / 2) + lR1 - LN2 - 0.5 * lr_a, 2.5),
           (lB_b + math.log(kq / 4) + lR1 - LN2 - 0.5 * lr_a, 0.5 + sh2)]
    CnB2, SnB2, CnB1, Mn = [], [], [], []
    if not G0:
        kr = lGv - math.log(4 * (1 - ETAB))
        Cp += [(lGv + math.log(FAC['low_S2']) + math.log((1 - NU['nu0']) / (4 * (1 + Z))) + 2 * np.log1p(-b_hi) - math.log(4) - lB_b - np.log1p(y_max / 2), 1 - sh2)]
        Cn += [(z0 + lGv + math.log(1 / (4 * (1 + Z))), 1), (z0 + math.log(NU['tnu2'] * QB / (16 * (1 + Z))), 1)]
        lW = lGv - math.log(2 * (1 - ETAB))
        Cn += [(lW + math.log(CZ / 2 * NU['tnu2'] / 2) - LN2 - lB_a, 1 - sh2),
               (lW + math.log(CZ / 2 * NU['stnu1']) - 0.5 * (LN2 + lB_a) + lX1 - 0.5 * LN2, 1 - sh1)]
        lB21 = np.log(math.sqrt(CZETA / 2) + r_hi / np.sqrt(math.pi * (2 + r_lo)))
        Sn += [(lW + math.log(CZ / 2 * NU['stnu1']) - 0.5 * (LN2 + lB_a) + lB21 - LN2 - 0.5 * lr_a, 1.5 - sh1)]
        CnB2 += [(z0 + kr + lX1 - 0.5 * LN2, 0)]
        SnB2 += [(kr + lB21 - LN2 - 0.5 * lr_a, 0.5)]
        lb10 = lse(math.log(math.sqrt(CZETA)) - 0.5 * np.log(2 * math.pi * (2 + r_lo)), lr_b - np.log(2 * (2 + r_lo)), z0 + 2 * lX1)
        CnB1 += [(kr + lb10, 1)]
        LnB1 += [(kr + np.logaddexp(z0 + math.log(CZETA), 2 * lr_b - np.log(2 * (2 + r_lo))) - math.log(4) - lr_a, 2)]
        Mn += [(kr + lX1 - 0.5 * np.log(2 * math.pi * (2 + r_lo)) + LN2 + 0.5 * lr_b, 0.5)]
        Cn += [(z0 + lgam + math.log((1 + NU['nu0']) * CZ / 2) + lGv + math.log(FAC['up_S2']), 0)]
    # ---- rounding (cl_link's rule)
    Bc = 2 * (np.maximum(np.abs(lmu_a), np.abs(lmu_b)) + np.maximum(np.abs(lr_a), np.abs(lr_b)) + np.maximum(np.abs(lB_a), np.abs(lB_b))) + 3 * ell_b + 600
    yc = np.where(y_max < 1500, y_max, 0.0)
    def eta(l):
        la = np.where(np.isfinite(l), np.abs(l), 0.0)
        return SAFE * (25 * U * (la + 3 * Bc + yc * Bc) + 1000 * U)
    def adj(lst, sgn):
        return [((np.asarray(t[0], dtype=float) + np.zeros(n)) + sgn * eta(np.asarray(t[0], dtype=float) + np.zeros(n)), t[1]) for t in lst]
    Lp = adj(Lp, -1); Cp = adj(Cp, -1); Hp = adj(Hp, -1); HpG = adj(HpG, -1)
    Ln = adj(Ln, 1); Cn = adj(Cn, 1); Hn = adj(Hn, 1); Sn = adj(Sn, 1); LnB1 = adj(LnB1, 1); CnB2 = adj(CnB2, 1); SnB2 = adj(SnB2, 1); CnB1 = adj(CnB1, 1); Mn = adj(Mn, 1)
    def worst_neg(p): return np.maximum(p * ls_lo, p * ls_hi) if not (VARIANT['V1'] and p < 0) else p * ls_hi
    def net(pos, neg): return net_powers(pos, neg, ls_lo, ls_hi, n)
    def wn(lst): return lse(*[t[0] + worst_neg(t[1]) for t in lst]) if lst else np.full(n, NEG)
    lbep_all, lben_all = net(Hp + HpG, Hn); lbep_noF, lben_noF = net(Hp, Hn)
    # (J6) option 'lo': the H bracket of power 2, c sigma^2/(a' + sigma), is INCREASING in sigma (x^p/(a + x) increases for p >= 1), so a positive net is
    # bounded below by c sigma_bot^2/(a' + sigma_bot); the power-0 G part (>= 0) is dropped.  Valid only where every power net of the bracket is positive.
    lbep_lo, lben_lo = net_powers(Hp, Hn, ls_lo, ls_lo, n)
    lbep_lo = np.where(np.isfinite(lben_lo), -np.inf, lbep_lo); lben_lo = np.full(n, NEG)
    expos = lambda l: np.exp(np.minimum(l, 700.0))
    exneg = lambda l: np.exp(np.minimum(l, 709.0)) + np.where(np.isfinite(l), 1e-300, 0.0)
    final = np.full(n, -np.inf)
    for keepF, hopt in [(k, h) for k in ((True, False) if not G0 else (False,)) for h in (('hi', 'lo') if HSHIFT_LO else ('hi',))]:
        res = {}
        lbep, lben = ((lbep_all, lben_all) if keepF else (lbep_noF, lben_noF)) if hopt == 'hi' else (lbep_lo, lben_lo)
        ls_sh = ls_hi if hopt == 'hi' else ls_lo
        for branch in (('B2', 'B1') if keepF else ('none',)):
            Lneg = Ln + (LnB1 if branch == 'B1' else [])
            Cneg = Cn + (CnB1 if branch == 'B1' else CnB2 if branch == 'B2' else [])
            Sneg = Sn + (SnB2 if branch == 'B2' else [])
            lalp, laln = net(Lp, Lneg); lkap, lkan = net(Cp, Cneg)
            lS = wn(Sneg); lM = wn(Mn) if branch == 'B1' else np.full(n, NEG)
            ok = (laln < lalp) & (lben < lbep)
            lal = lalp + np.log1p(-np.where(laln < lalp, np.exp(np.minimum(laln - lalp, 0.0)), 0.0))
            lbe = lbep + np.log1p(-np.where(lben < lbep, np.exp(np.minimum(lben - lbep, 0.0)), 0.0))
            lsc = 0.5 * (lal + lbe); lu = 0.5 * (lbe - lal)
            lSu = lS + 0.5 * lu; lMu = lM - 0.5 * lu
            Lst = np.maximum.reduce([lsc, lkap, np.where(np.isfinite(lkan), lkan, -np.inf), np.where(np.isfinite(lSu), lSu, -np.inf)])
            Lst = np.where(np.isfinite(Lst), Lst, 0.0)
            w = np.exp(lsc - Lst)
            kpos = expos(lkap - Lst); kneg = exneg(lkan - Lst); kap = kpos - kneg
            Sp = exneg(lSu - Lst); Mp = exneg(lMu - Lst)
            sh = np.exp(np.clip(ls_sh - lu, -700.0, 700.0))                     # H shift: sigma_top (option hi) or sigma_bot (option lo, J6)
            lp0 = math.log(16) + lr_a + 2 * ls_lo - lu
            p = np.where(lp0 < -700, 0.0, np.exp(np.minimum(lp0, 700.0)))
            lq = math.log(4) + lr_b - ls_lo - lu
            q = np.where(lq > 700, np.inf, np.exp(np.maximum(lq, -700.0)))
            best = {}
            lp_ = np.where(p > 0, np.log(np.maximum(p, 1e-308)), -np.inf); lq_ = lq; lsh_ = ls_sh - lu
            for tstar in ((1e-9, 1e-3, 0.1, 1.0) if branch != 'none' else (None,)):
                ts = np.clip(np.full(n, tstar), p, q) if tstar is not None else q
                lts = np.clip(np.full(n, math.log(tstar)), lp_, lq_) if tstar is not None else lq_
                vb = np.full(n, -np.inf)
                for frac in (0.0, 1e-2, 0.3):
                    up = ts if branch == 'B2' else q
                    with np.errstate(all='ignore'):
                        Sterm = (np.where(w > 0, Sp * Sp / (4 * frac * np.maximum(w, 1e-300)), np.where(Sp > 1e-299, np.inf, 0.0)) if frac > 0 else Sp * np.sqrt(up))
                        Mterm = Mp / np.sqrt(ts) if branch == 'B1' else 0.0
                    Kk = kap - Sterm - Mterm
                    mag0 = kpos + kneg + Sterm + Mterm
                    if branch == 'B1':
                        v = K.interval_min(w * (1 - frac), w, Kk, ts, q, sh, mag0)
                    else:
                        v = K.interval_min(w * (1 - frac), w, Kk, p, ts if branch == 'B2' else q, sh, mag0)
                    vb = np.maximum(vb, v)
                lqb = lts if branch == 'B2' else lq_
                lpos = np.logaddexp(lsc - np.logaddexp(lqb, lsh_), lkap)
                negs = [lkan, lS + 0.5 * (lu + lqb)]
                if branch == 'B1': negs.append(lM - 0.5 * (lu + lts))
                lneg = lse(*[np.where(np.isfinite(x), x, -np.inf) for x in negs])
                with np.errstate(all='ignore'):
                    vlog = np.where(lpos > lneg, -np.expm1(lneg - lpos), -1.0)
                vb = np.maximum(vb, np.where(np.isfinite(lpos), vlog, -np.inf))
                best[tstar] = np.where(ok, vb, -np.inf)
            res[branch] = best
        if keepF:
            for tstar in res['B1']:
                final = np.maximum(final, np.minimum(res['B1'][tstar], res['B2'][tstar]))
        else:
            final = np.maximum(final, res['none'][None])
    return np.where(empty, np.inf, final)


def cell_ok(lmu_a, lmu_b, lr_a, lr_b, lB_a, lB_b, ell_a, ell_b, Dd):
    lPi = np.where(lmu_a >= math.log(2.0), CC.lPistar_vec(np.maximum(lmu_a, math.log(2.0))), lmu_b + 2 * math.log(QB)) + 1e-6
    lGc = lPi + 0.5 * (lr_b - np.log(2 + np.exp(lr_b)))
    v0 = phi_cells_beta(lmu_a, lmu_b, lr_a, lr_b, lB_a, lB_b, ell_a, ell_b, Dd, -np.inf)
    v1 = phi_cells_beta(lmu_a, lmu_b, lr_a, lr_b, lB_a, lB_b, ell_a, ell_b, Dd, lGc)
    return (v0 > 1e-10) & (v1 > 1e-10)


def initial_boxes(ell_a, ell_b, lBmin=math.log(1e4), nmu=16, nr=24, nB=24):
    """the initial table: mu in [1e-4, max(4 ell_b, 50)] (16 boxes) plus the two mu tails, r in [sigma_bot^3, 1/4] (24), beta in [beta_min, b_max/sigma_bot^2] (24)."""
    lsig3 = -3 * ell_b; lBmax = math.log(0.003 / K.TAU_MIN) + 2 * ell_b
    mu_lo, mu_hi = 1e-4, max(4.0 * ell_b, 50.0)
    lm = np.linspace(math.log(mu_lo), math.log(mu_hi), nmu + 1); lr = np.linspace(lsig3, math.log(0.25), nr + 1); lB = np.linspace(lBmin, lBmax, nB + 1)
    A = np.array([(a, b, c, d, e, f) for a, b in zip(lm[:-1], lm[1:]) for c, d in zip(lr[:-1], lr[1:]) for e, f in zip(lB[:-1], lB[1:])])
    tails = [(-30.0, math.log(mu_lo), lsig3, math.log(0.25), lBmin, lBmax), (math.log(mu_hi), math.log(mu_hi) + 60, lsig3, math.log(0.25), lBmin, lBmax)]
    return np.vstack([A, np.array(tails)])


def certify(ell_a, ell_b, Dd, lBmin=math.log(1e4), max_cells=400_000_000, chunk=200_000):
    stack = [initial_boxes(ell_a, ell_b, lBmin)]; total = 0; t0 = time.time(); nempty = 0
    while stack:
        X = stack.pop()
        if len(X) > chunk: stack.extend(np.array_split(X, int(math.ceil(len(X) / chunk)))); continue
        total += len(X)
        if total > max_cells:
            return {'ell_a': ell_a, 'ell_b': ell_b, 'D': Dd, 'closed': False, 'cells': total, 'reason': 'cell budget', 'seconds': round(time.time() - t0, 1)}
        with np.errstate(all='ignore'):
            ok = cell_ok(X[:, 0], X[:, 1], X[:, 2], X[:, 3], X[:, 4], X[:, 5], ell_a, ell_b, Dd)
        bad = X[~ok]
        if not len(bad): continue
        wdt = np.stack([bad[:, 1] - bad[:, 0], (bad[:, 3] - bad[:, 2]) * 0.7, (bad[:, 5] - bad[:, 4]) * 0.7], axis=1)
        k = np.argmax(wdt, axis=1); ms = float(np.min(np.max(wdt, axis=1)))
        if ms < 1e-6:
            j = int(np.argmin(np.max(wdt, axis=1)))
            return {'ell_a': ell_a, 'ell_b': ell_b, 'D': Dd, 'closed': False, 'cells': total, 'reason': 'cells too small', 'stuck_cell': bad[j].tolist(),
                    'stuck_log_beta': float(bad[j, 4]), 'seconds': round(time.time() - t0, 1)}
        new = []
        for jj in range(3):
            sel = bad[k == jj]
            if not len(sel): continue
            mid = (sel[:, 2 * jj] + sel[:, 2 * jj + 1]) / 2
            L = sel.copy(); R = sel.copy(); L[:, 2 * jj + 1] = mid; R[:, 2 * jj] = mid; new += [L, R]
        stack.append(np.vstack(new))
    return {'ell_a': ell_a, 'ell_b': ell_b, 'D': Dd, 'lBmin': lBmin, 'tau_min': K.TAU_MIN, 'closed': True, 'cells_evaluated': total, 'seconds': round(time.time() - t0, 1),
            'rule': 'CL rounding-covering, S = 10, beta boxes'}


if __name__ == '__main__':
    ef, et, st = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]); Dd = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0
    lB0 = float(sys.argv[5]) if len(sys.argv) > 5 else math.log(1e4); tmin = float(sys.argv[6]) if len(sys.argv) > 6 else 20.0
    K.set_tau_min(tmin)
    have = set()
    if os.path.exists(FN):
        for l in open(FN):
            if l.strip(): d = json.loads(l); have.add((round(d['ell_a'], 6), d['D'], round(d.get('lBmin', 0), 6), d.get('tau_min')))
    ea = ef
    while ea < et - 1e-9:
        eb = min(et, ea + st)
        if (round(ea, 6), Dd, round(lB0, 6), tmin) not in have:
            d = certify(ea, eb, Dd, lB0); d['lBmin'] = lB0; d['tau_min'] = tmin
            open(FN, 'a').write(json.dumps(d) + '\n')
            print('ell [%.2f, %.2f] D %g tau >= %g beta >= e^%.2f: %s cells %s %.0fs %s' % (ea, eb, Dd, tmin, lB0, 'CLOSED' if d['closed'] else 'OPEN (' + d.get('reason', '') + ')',
                  d.get('cells_evaluated', d.get('cells')), d['seconds'], d.get('stuck_log_beta', '')), flush=True)
        ea = eb
