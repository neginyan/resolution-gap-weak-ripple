"""Stage CI 2.2: the link cells of ce_link.py with LEMMA A* inserted (a_5, a_6 without the window-generic gamma_P).  Python 3.11.

CHANGE AGAINST ce_link.py (everything else identical).  Lemma A* (report CH 2.1, CI 2): in BL3 every remainder nu * sup_t|g| may use nu * <|g|>
(the kernels Theta - 1, theta_tau - Theta, d_t Theta have absolutely summable Fourier series with the same sums; |g_n| <= <|g|>).  Hence
  a_5* = sigma^2 nu_0 int Gamma <h |sigma^2 G|> <= sigma^2 nu_0 sqrt(M_G),
  a_6* = nu_1 [ x int Gamma <|h_t|> chi_bar + sigma sqrt(x) int Gamma <h |sigma^2 G|>/sqrt2 ]
       <= nu_1 sigma sqrt(x) [ sqrt(G/2) sqrt(2 c_a^2 + 1/4) + sqrt(M_G/2) ],
with chi_bar <= sigma c_a + sigma^2 |w|/(sqrt2 L^2), c_a = (U1 + Qb e^{-1/2})/((1 - Qb) sqrt2), int Gamma E chi_bar <= sqrt(g_s) sqrt(int Gamma chi_bar^2),
and (Lemma BE-2) M_G <= 1.5 (K_1 G + g_s T/(1 - Qb)), K_1 = Lambda/(1-Qb) + (U1^2 + 1/2)/(1-Qb)^3, Lambda = 4 ell + 60,
T <= 4 e^{-Lambda} (Lambda^2 + Lambda^1.5 e^{ell}/(1.25 sqrt(5e11))), so sqrt(M_G) <= sqrt(1.5 K_1 G) + sqrt(2 T).
The sqrt(G) is made affine in G on each cell by sqrt(G) <= (G + t)/(2 sqrt t), t = G_cell/16 (valid for every G >= 0), so the
G in {0, G_cell} argument of ce_link.py is unchanged.  -> results/ci_link_cells.jsonl.

--- original docstring of ce_link.py ---
Stage CE 2.2: the link to Lemma R as a COMPUTER-ASSISTED PROOF.  Python 3.11.

CLAIM (certified on cells).  For every state with tau = x/b >= tau_min, x <= 0.003, x >= 5e11 sigma^2 (so beta <= 1e-6 tau, as in CB), b in
[B0 sigma^2, 1.5e-4], every L^2 <= 1/sigma, sigma^2 <= s^2 <= L^2/4, every mu > 0 and every ripple G = P g_s in [0, Pi*(mu) g_s(r)] (x-free Lemma Z3'),
    Num_R* - gamma sigma den_R* >= 0,   gamma = T_D(ell)/sqrt(ell),  T_D = 1 + (2.1 + D)/ell + 2.5e-4  (D = 10),
where Num_R*/den_R* is the FULL lower bound of Lemma R* (Lemma R of CA with a_4 -> a_4*, ce_a4.py): every remainder of Lemma R is kept.  Hence
g >= (1 - eps) (lambda/Upsilon) T_D/sqrt(ell) with the sandwich factors restored by the scaling argument (see below).  The states with b < B0 sigma^2
(the b-edge) are NOT covered here.

WINDOW CONSTANTS (remainders only).  Qbar = 1/4.  In remainder terms eta <= etabar = 1/2 (factors 1/(1 - etabar) = 2, zeta = sigma sqrt(c_zeta A),
c_zeta = etabar^2/(e (1-etabar)^2)).  q0_n <= 2 e^{-n^2 x/2} (|a_n| <= 2 by non-negativity, blur K >= x) and sum q0_n <= Qbar give, with
U^2 = ln(1/x) + 30 <= 2 ell + 3.1:  U1 = sqrt(x) N_1 <= Qbar sqrt(2 ell + 3.1) + 1e-6,  U2 = x N_2 <= Qbar (2 ell + 3.1) + 1e-4  (the tails over
u_n = n sqrt(x) > U of u e^{-u^2/2} and u^2 e^{-u^2/2}, decreasing there, are bounded by the integral plus one term: < 2e^{-15}(1 + U sqrt x) etc.).
gamma_P and c_chi of BL3 follow from U1, U2.  Every remainder carries nu_k(tau, beta) tau^j with j <= 2; these decrease on tau >= 20 (checked,
ce_link_checks.py), so they are taken at tau = tau_min; x = b tau inside a_6, a_7, a_9, a_10.

MAIN-FACTOR SCALING.  The sandwich factors of the main terms are set to 1 in the cells; for a window with lower factors >= lambda and upper factors
<= Upsilon (eta sandwich, or Lemma S / S'), Num >= lambda (positive part) - (negative remainders, already at etabar) >= lambda (cell expression) and
den_(S2 + background) <= Upsilon den_mid with Upsilon = up_S2 (the S0 part of den carries its factor <= 1/(1 - etabar) = 2 inside the cells),
so g >= (lambda/Upsilon) gamma, lambda = min(low_S0, low_S2).

THE REDUCED FORM (proof of each line in the comments of phi_cells).  With A = sigma^2/s^2 in [4 sigma^3, 1], r = s^2/L^2, c = sigma^2/(4r) (L^2/4 = c/A),
a' = 4 r A/sigma (so A = sigma a'/(4r) and c/(A + c) = sigma/(a' + sigma)), Phi/sigma is bounded below by
    alpha a' + beta/(a' + sigma) + kappa - S sqrt(a') [- M/sqrt(a') on a' >= a*],
every coefficient being a finite sum of (sigma-free monotone piece) * sigma^p, p >= 0.  Key inequalities: Psi_lo >= w/(1+w) - e(1-y)^+ r_C with
w = 3y/2 + L^2/4 ((1-y)e^{-y/2-L^2/4} <= e^{-w}); w/(1+w) >= w0/(1+w0) + (c/(1+w0)^2)/(A + c); 1 - e^{-L^2/4} <= 2c/(A + c);
1 - e^{-se2/4} >= c1/(A + c1) >= c1/(A + c) (c1 = sigma^2/(4(2+r)) <= c); r_C = rho0 + rho1 sqrt(A); F_rem <= (G/(2(1-etabar))) min(B1, B2)
(B2 from |cos(a+d) - cos a| <= |d|, B1 the v1.5 form); (1 - e) >= (y/2)/(1 + y/2); A <= 1 for A^2 terms; the F-term e sigma^2 <theta F (1 - cos psi)> is
>= 0, so it may be dropped together with its remainder (option keepF = False).  Same powers of sigma are netted at EQUAL sigma before the worst end of
the ell cell is taken; the beta bracket (b-term, F-main, minus sigma^2 gamma_P cz and the nu-part of Phi_up) is formed at equal sigma and then bounded
with shift sigma_hi.  The minimum over a' uses S sqrt(t) <= theta t + S^2/(4 theta) or S sqrt(t) <= S sqrt(t_max); the a'-range is split at a* (B2 below,
B1 above; the same a* on both sides); all in logarithms, normalised by the largest part (underflow only drops non-negative parts; negative parts are
never rounded to 0).  Acceptance: the minimum divided by the total magnitude of its terms exceeds 1e-10.
G: for each A the lower bound is affine in G, so its minimum over A is concave in G and G in {0, G_cell} suffices; at G = 0 the keepF options coincide,
so at G_cell either option may be used.  Cells in (log mu, log r, log b); tails mu in [e^-30, 1e-4] and [4 ell, 4 ell e^60].
-> results/ce_link_cells.jsonl (one line per ell cell; resumable).  CLI: python ce_link.py ELL_FROM ELL_TO STEP [D] [log B0] [tau_min].
"""
import json, math, os, sys, time
import numpy as np
import cc_cells as CC, cc_core as C

FN = 'results/ci_link_cells.jsonl'
QB = 0.25; ETAB = 0.5; LN2 = math.log(2.0)
CZETA = ETAB ** 2 / (math.e * (1 - ETAB) ** 2)


def nus_tau(t=20.0):
    n = np.arange(1, 200); bt = 1e-6 * t
    nu = lambda k: float(2 * np.sum(n ** k * np.exp(-n * n * t / 2 + n * bt)))
    om2 = float(2 * np.sum(np.exp(-n * n * t / 2) * np.abs(1 - n * n * t)))
    return {'nu0': nu(0), 'nu1': nu(1), 'nu2': nu(2), 'nu3': nu(3), 'nu4': nu(4), 'tnu2': t * nu(2), 't2nu4': t * t * nu(4), 'stnu1': math.sqrt(t) * nu(1),
            't15nu3': t ** 1.5 * nu(3), 'om2': om2}


TAU_MIN = 20.0
FAC = {'low_S0': 1.0, 'low_S2': 1.0, 'up_S2': 1.0}     # main sandwich factors per channel (1 = mid; the scaling argument restores them); set per window in cf_windows.py
NU = nus_tau(TAU_MIN)
Z = QB * NU['nu0'] / 2; CZ = 1 / (1 - Z); EPS1 = NU['nu0'] * (1 + 2 * QB) / (1 - Z)


def set_tau_min(t):
    """every tau-weighted nu is evaluated at tau = tau_min (they decrease on tau >= 20, checked in ce_link_checks)."""
    global TAU_MIN, NU, Z, CZ, EPS1
    TAU_MIN = float(t); NU = nus_tau(TAU_MIN); Z = QB * NU['nu0'] / 2; CZ = 1 / (1 - Z); EPS1 = NU['nu0'] * (1 + 2 * QB) / (1 - Z)


def window_consts(ell):
    U1 = QB * math.sqrt(2 * ell + 3.1) + 1e-6; U2 = QB * (2 * ell + 3.1) + 1e-4
    gP = (U2 * (1 + QB) + U1 ** 2 + U1 * (1 + 2 * QB) / math.sqrt(math.e)) / (1 - QB) ** 2
    cchi = (U1 / (1 - QB) + math.sqrt(2 / math.pi) + QB * math.exp(-0.5) / (1 - QB)) / math.sqrt(2)
    return U1, U2, gP, cchi


def lse(*arrs):
    a = np.stack(arrs); m = np.max(a, axis=0); m = np.where(np.isfinite(m), m, 0.0)
    return m + np.log(np.sum(np.exp(a - m), axis=0))


def phi_cells(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, Dd, lG, lGt=None):
    """vectorised lower bound of min over A of Phi/sigma, Phi = Num_R* - gamma sigma den_R*, scaled by its largest positive term (sign is what matters).
    Normalisation: a' = 4 r A/sigma (A = sigma a'/(4r)), so that the leading terms are sigma-free (as in CC); every term is
    (sigma-free coefficient) * sigma^p, and sigma^p is taken at the worst end of the ell cell (p >= 0 throughout).  Types (proof in the module docstring):
      L: k a',  H: k/(a' + sigma),  C: k,  S: -k sqrt(a'),  M: -k/sqrt(a') (B1 branch, a' >= a*).
    The H-bracket is formed at equal sigma (b-term, F-main term, minus the sigma^2 gamma_P and nu-terms of Phi_up) and then bounded with shift sigma_hi
    (if the bracket can be negative: shift sigma_lo)."""
    n = len(lmu_a); NEG = -np.inf
    G0 = np.isscalar(lG) and not np.isfinite(lG)
    s_lo, s_hi = math.exp(-ell_b), math.exp(-ell_a)
    ls_lo, ls_hi = -ell_b, -ell_a
    gam = (1 + (2.1 + Dd) / ell_a + 2.5e-4) / math.sqrt(ell_a); lgam = math.log(gam)
    U1, U2, gP, cchi = window_consts(ell_b); X1 = U1 / (1 - QB); lX1 = math.log(X1)
    r_lo = np.exp(lr_a); r_hi = np.exp(lr_b); b_lo = np.exp(lb_a); b_hi = np.exp(lb_b)
    ly_max = 2 * ls_hi - LN2 - lb_a; ly_min = 2 * np.log1p(-b_hi) + 2 * ls_lo - LN2 - lb_b
    y_max = np.exp(np.minimum(ly_max, 700)); y_min = np.exp(ly_min)
    e_min = np.exp(-y_max / 2)
    kF_lo = e_min * (1 - NU['nu0']) / (1 + Z)
    w0_max = 1.5 * y_max
    if G0:
        lg0_lo = lg0_hi = np.full(n, NEG); lGv = NEG
    else:
        lGv = lG; lg0_lo = lG - LN2 - lmu_b - np.log(2 + r_hi); lg0_hi = lG - LN2 - lmu_a - np.log(2 + r_lo)
    z0 = np.zeros(n)
    Lp, Ln, Cp, Cn, Sn = [], [], [], [], []
    kv = (1 + NU['om2']) * CZ / 4; kw = 1.5 * (1 + NU['nu0']) * CZ; kq = CZ * QB * (NU['nu0'] + NU['tnu2'])
    # ---- L: alpha A/sigma = alpha a'/(4r)
    Lp += [(z0 + math.log((1 - EPS1) / 2), 0), (math.log(1 - EPS1) + lg0_lo - LN2 - lr_b + math.log(FAC['low_S0']), 0)]
    Ln += [(math.log(kv / 4) - lb_a, 2), (math.log(kv / 4) - lb_a + lg0_hi - lr_a, 2),
           (z0 + math.log(kw * (X1 ** 2 + CZETA) / 4), 2), (math.log(kw * (X1 ** 2 + CZETA) / 4) + lg0_hi - lr_a, 2),
           (math.log(3 * kw / 4) + lr_b, 2), (math.log(3 * kw / 4) + lr_b + lg0_hi - np.log(2 + r_lo), 2),
           (z0 + lgam + math.log((1 + EPS1) / 4), 1), (lgam + math.log(2 * (1 + EPS1) / 4) + lg0_hi - lr_a, 1)]   # S0 part of den: upper sandwich factor <= 1/(1 - etabar) = 2 kept here, so that Upsilon = up_S2 only
    LnB1 = []
    # ---- H bracket at equal sigma (logs): b/(2(1+z)(1+w0)^2) + kF G r/(4(2+r)) - kq b/2 - sigma^2 gamma_P cz
    Hp = [lb_a - 2 * np.log1p(w0_max) - math.log(2 * (1 + Z))]
    if not G0:
        Hp += [np.log(kF_lo) + lGv + lr_a - np.log(2 + r_lo) - math.log(4) + math.log(FAC['low_S2'])]
    Hn = [lb_b + math.log(kq / 2), z0 + math.log(gP * CZ) + 2 * ls_hi]
    # ---- C (Phi/sigma constants)
    Cp += [(math.log(3 / (8 * (1 + Z))) + 2 * np.log1p(-b_hi) - np.log1p(w0_max), 1)]
    Cn += [(z0 + math.log(kw), 1), (math.log(kw) + lg0_hi - np.log(2 + r_lo), 1)]
    has = (1 - y_min) > 0
    Cn += [(np.where(has, lb_b + np.log(np.maximum(1 - y_min, 1e-300)), NEG) + lX1 - math.log(2 * math.sqrt(2) * (1 + Z)), 0)]
    Cn += [(z0 + math.log(gP * CZ * X1 / (2 * math.sqrt(2))), 2)]
    Cn += [(lb_b + math.log(kq * X1 / (4 * math.sqrt(2))), 0)]
    sb = 0.5 * lb_b
    Cn += [(z0 + math.log(CZ / 2 * 3 / 8 * QB * NU['nu0']), 1),
           (lb_b + math.log(CZ / 2 * (1 + QB) * cchi * NU['nu0'] * (1 + QB / 2)), 0),
           (z0 + math.log(CZ / 2 * QB * NU['tnu2'] / 8), 1),
           (lb_b + math.log(CZ / 2 * 0.5 * QB * (1 + QB) * cchi * NU['tnu2']), 0),
           (z0 + math.log(CZ / 2 * QB * NU['t2nu4'] / 8), 1),
           (lb_b + math.log(CZ / 2 * QB * NU['tnu2'] / (2 * (1 - NU['nu0']))), 1),
           (sb + math.log(CZ / 2 * QB * NU['t15nu3'] / (2 * math.sqrt(2))), 0)]
    Cn += [(z0 + lgam + lb_b + math.log(1 + NU['tnu2'] * QB / (2 * (1 - Z))), 0)]
    # ---- Lemma A*: a_5*, a_6* (replace the gamma_P forms of a_5, a_6)
    Lam = 4 * ell_b + 60.0; K1 = Lam / (1 - QB) + (U1 ** 2 + 0.5) / (1 - QB) ** 3
    c_a = (U1 + QB * math.exp(-0.5)) / ((1 - QB) * math.sqrt(2))
    c6G = math.sqrt(2 * c_a ** 2 + 0.25) / math.sqrt(2) + math.sqrt(0.75 * K1); c5 = math.sqrt(1.5 * K1)
    lT = math.log(4) - Lam + float(np.logaddexp(2 * math.log(Lam), 1.5 * math.log(Lam) + ell_b - math.log(1.25 * math.sqrt(5e11))))
    Cn += [(sb + math.log(CZ / 2 * NU['stnu1']) + 0.5 * lT, 0), (z0 + math.log(CZ / 2 * NU['nu0']) + 0.5 * (LN2 + lT), 1)]
    if lGt is not None:
        lt = lGt
        Cn += [(sb + math.log(CZ / 2 * NU['stnu1'] * c6G) + 0.5 * lt - LN2, 0), (z0 + math.log(CZ / 2 * NU['nu0'] * c5) + 0.5 * lt - LN2, 1)]
        if not G0:
            Cn += [(sb + math.log(CZ / 2 * NU['stnu1'] * c6G) + lGv - 0.5 * lt - LN2, 0), (lGv + math.log(CZ / 2 * NU['nu0'] * c5) - 0.5 * lt - LN2, 1)]
    # ---- S: S sqrt(A)/sigma = S/(2 sqrt(r sigma)) sqrt(a'); S = sigma^k * free -> power k - 1/2
    lR1 = np.log(math.sqrt(CZETA) + np.sqrt(2 * r_hi / math.pi)) - 0.5 * LN2
    Sn += [(np.where(has, lb_b + np.log(np.maximum(1 - y_min, 1e-300)), NEG) + lR1 - math.log(2 * (1 + Z)) - LN2 - 0.5 * lr_a, 0.5),
           (math.log(gP * CZ / 2) + lR1 - LN2 - 0.5 * lr_a, 2.5),
           (lb_b + math.log(kq / 4) + lR1 - LN2 - 0.5 * lr_a, 0.5)]
    CnB2, SnB2, CnB1, Mn = [], [], [], []
    if not G0:
        kr = lGv - math.log(4 * (1 - ETAB))
        Cp += [(lGv + math.log(FAC['low_S2']) + math.log((1 - NU['nu0']) / (4 * (1 + Z))) + 2 * np.log1p(-b_hi) - math.log(4) - lb_b - np.log1p(y_max / 2), 1)]
        Cn += [(z0 + lGv + math.log(1 / (4 * (1 + Z))), 1), (z0 + math.log(NU['tnu2'] * QB / (16 * (1 + Z))), 1)]
        lW = lGv - math.log(2 * (1 - ETAB))
        Cn += [(lW + math.log(CZ / 2 * NU['tnu2'] / 2) - LN2 - lb_a, 1),
               (lW + math.log(CZ / 2 * NU['stnu1']) - 0.5 * (LN2 + lb_a) + lX1 - 0.5 * LN2, 1)]
        lB21 = np.log(math.sqrt(CZETA / 2) + r_hi / np.sqrt(math.pi * (2 + r_lo)))
        Sn += [(lW + math.log(CZ / 2 * NU['stnu1']) - 0.5 * (LN2 + lb_a) + lB21 - LN2 - 0.5 * lr_a, 1.5)]
        CnB2 += [(z0 + kr + lX1 - 0.5 * LN2, 0)]
        SnB2 += [(kr + lB21 - LN2 - 0.5 * lr_a, 0.5)]
        lb10 = lse(math.log(math.sqrt(CZETA)) - 0.5 * np.log(2 * math.pi * (2 + r_lo)), lr_b - np.log(2 * (2 + r_lo)), z0 + 2 * lX1)
        CnB1 += [(kr + lb10, 1)]
        LnB1 += [(kr + np.logaddexp(z0 + math.log(CZETA), 2 * lr_b - np.log(2 * (2 + r_lo))) - math.log(4) - lr_a, 2)]
        Mn += [(kr + lX1 - 0.5 * np.log(2 * math.pi * (2 + r_lo)) + LN2 + 0.5 * lr_b, 0.5)]
        Cn += [(z0 + lgam + math.log((1 + NU['nu0']) * CZ / 2) + lGv + math.log(FAC['up_S2']), 0)]
    def net(pos, neg):
        """sum of (positive - negative) terms, netted at EQUAL sigma within each power p, then each net taken at its worst sigma end.
        Returns (log of total positive part, log of total negative part)."""
        P, N = [], []
        for p in sorted(set([t[1] for t in pos] + [t[1] for t in neg])):
            lp = lse(*[t[0] for t in pos if t[1] == p]) if any(t[1] == p for t in pos) else np.full(n, NEG)
            ln = lse(*[t[0] for t in neg if t[1] == p]) if any(t[1] == p for t in neg) else np.full(n, NEG)
            d = np.where(lp > ln, lp + np.log1p(-np.exp(np.minimum(ln - lp, 0.0))), NEG)        # positive net
            e_ = np.where(ln >= lp, ln + np.log1p(-np.exp(np.minimum(lp - ln, 0.0))), NEG)      # negative net
            e_ = np.where((ln >= lp) & ~np.isfinite(lp), ln, e_)
            P.append(d + p * ls_lo); N.append(e_ + p * ls_hi)
        return lse(*P), lse(*N)
    # ---- assemble in logarithms and scale by the AM-GM balance sqrt(alpha beta): a' = u t, u = sqrt(beta/alpha)
    #      branch B2 covers a' in [a_lo, u t*], branch B1 covers [u t*, a_hi]; the SAME t* is used on both sides.
    #      option keepF: keep the F-term with its remainder (B1/B2); or drop both (the F-term e sigma^2 <theta F (1 - cos psi)> is >= 0).
    lbep_all = lse(*Hp); lbep_noF = Hp[0]; lben = lse(*Hn)
    expos = lambda l: np.exp(np.minimum(l, 700.0))                    # positive parts: clipping from above only lowers them (conservative)
    exneg = lambda l: np.exp(np.minimum(l, 709.0)) + np.where(np.isfinite(l), 1e-300, 0.0)   # negative parts: never rounded down to 0
    final = np.full(n, -np.inf)
    for keepF in ((True, False) if not G0 else (False,)):
        res = {}
        for branch in (('B2', 'B1') if keepF else ('none',)):
            Lneg = Ln + (LnB1 if branch == 'B1' else [])
            Cneg = Cn + (CnB1 if branch == 'B1' else CnB2 if branch == 'B2' else [])
            Sneg = Sn + (SnB2 if branch == 'B2' else [])
            lalp, laln = net(Lp, Lneg); lkap, lkan = net(Cp, Cneg)
            lS = lse(*[t[0] + t[1] * ls_hi for t in Sneg]) if Sneg else np.full(n, NEG)
            lM = lse(*[t[0] + t[1] * ls_hi for t in Mn]) if branch == 'B1' else np.full(n, NEG)
            lbep = lbep_all if keepF else lbep_noF
            ok = (laln < lalp) & (lben < lbep)
            lal = lalp + np.log1p(-np.where(laln < lalp, np.exp(np.minimum(laln - lalp, 0.0)), 0.0))
            lbe = lbep + np.log1p(-np.where(lben < lbep, np.exp(np.minimum(lben - lbep, 0.0)), 0.0))
            lsc = 0.5 * (lal + lbe); lu = 0.5 * (lbe - lal)
            # f(t) = sc [ (1 - frac) t + 1/(t + sh) ] + kappa - S sqrt(u t) - M/sqrt(u t); normalise by Lst = the largest part
            lSu = lS + 0.5 * lu; lMu = lM - 0.5 * lu
            Lst = np.maximum.reduce([lsc, lkap, np.where(np.isfinite(lkan), lkan, -np.inf), np.where(np.isfinite(lSu), lSu, -np.inf)])
            Lst = np.where(np.isfinite(Lst), Lst, 0.0)
            w = np.exp(lsc - Lst)                                         # <= 1; underflow drops the non-negative t-part: conservative
            kpos = expos(lkap - Lst); kneg = exneg(lkan - Lst); kap = kpos - kneg
            Sp = exneg(lSu - Lst); Mp = exneg(lMu - Lst)
            sh = np.exp(np.clip(ls_hi - lu, -700.0, 700.0))
            p = np.exp(np.minimum(math.log(16) + lr_a + 2 * ls_lo - lu, 700.0)); p = np.where(math.log(16) + lr_a + 2 * ls_lo - lu < -700, 0.0, p)
            lq = math.log(4) + lr_b - ls_lo - lu
            q = np.where(lq > 700, np.inf, np.exp(np.maximum(lq, -700.0)))
            best = {}
            lp_ = np.where(p > 0, np.log(np.maximum(p, 1e-308)), -np.inf); lq_ = lq; lsh_ = ls_hi - lu
            for tstar in ((1e-9, 1e-3, 0.1, 1.0) if branch != 'none' else (None,)):
                ts = np.clip(np.full(n, tstar), p, q) if tstar is not None else q
                lts = np.clip(np.full(n, math.log(tstar)), lp_, lq_) if tstar is not None else lq_
                vb = np.full(n, -np.inf)
                for frac in (0.0, 1e-2, 0.3):
                    up = ts if branch == 'B2' else q
                    with np.errstate(all='ignore'):
                        Sterm = (np.where(w > 0, Sp * Sp / (4 * frac * np.maximum(w, 1e-300)), np.where(Sp > 1e-299, np.inf, 0.0)) if frac > 0 else Sp * np.sqrt(up))
                        Mterm = Mp / np.sqrt(ts) if branch == 'B1' else 0.0
                    K = kap - Sterm - Mterm
                    mag0 = kpos + kneg + Sterm + Mterm
                    if branch == 'B1':
                        v = interval_min(w * (1 - frac), w, K, ts, q, sh, mag0)
                    else:
                        v = interval_min(w * (1 - frac), w, K, p, ts if branch == 'B2' else q, sh, mag0)
                    vb = np.maximum(vb, v)
                # degenerate-interval option, all in logs: on [p_b, q_b], (1-frac) t + 1/(t + sh) >= 1/(q_b + sh); S sqrt(u t) <= S sqrt(u q_b); M/sqrt(u t) <= M/sqrt(u t*)
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
    return final


def interval_min(al, be, ka, p, q, s, mag0=None):
    """min over a in [p, q] of al a + be/(a + s) + ka (convex when be > 0, al > 0).  With mag0 (the sum of the magnitudes of the constant parts),
    returns the minimum divided by the total magnitude of the terms at the minimiser (rounding errors are ~1e-15 of that)."""
    conv = (be > 0) & (al > 0)
    au = np.where(conv, np.sqrt(np.maximum(be, 0) / np.maximum(al, 1e-300)) - s, 0.0)
    ac = np.clip(au, p, q)
    with np.errstate(all='ignore'):
        v_conv = al * ac + be / (ac + s) + ka
        v_nc = np.where(al >= 0, al * p, al * q) + np.where(be >= 0, be / (q + s), be / (p + s)) + ka
        v = np.where(conv, v_conv, v_nc)
        if mag0 is None: return v
        m = mag0 + np.where(conv, np.abs(al * ac) + np.abs(be / (ac + s)), np.abs(np.where(al >= 0, al * p, al * q)) + np.abs(be) / (p + s))
        return np.where(np.isfinite(v) & (m > 0), v / m, -np.inf)


def cell_ok(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, Dd):
    lPi = np.where(lmu_a >= math.log(2.0), CC.lPistar_vec(np.maximum(lmu_a, math.log(2.0))), lmu_b + 2 * math.log(QB)) + 1e-6
    lGc = lPi + 0.5 * (lr_b - np.log(2 + np.exp(lr_b)))
    ok = np.ones(len(lmu_a), dtype=bool)
    lGt = lGc - math.log(16.0)
    v0 = phi_cells(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, Dd, -np.inf, lGt)
    ok &= v0 > 1e-10
    # G = G_cell (array): evaluate term lists with a per-cell lG
    v1 = phi_cells(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, Dd, lGc, lGt)
    ok &= v1 > 1e-10
    return ok


def certify(ell_a, ell_b, Dd, lbmin_rel, max_cells=400_000_000, chunk=1_000_000):
    """cells in (log mu, log r, log b); b in [b_min, 1.5e-4] with b_min = sigma^2 * e^{lbmin_rel} (the edge below is NOT covered); r in [sigma^3, 1/4]."""
    lsig3 = -3 * ell_b; lbmin = -2 * ell_a + lbmin_rel; lbmax = math.log(0.003 / TAU_MIN)          # b = x/tau <= 0.003/tau_min
    mu_lo, mu_hi = 1e-4, max(4.0 * ell_b, 50.0)
    lm = np.linspace(math.log(mu_lo), math.log(mu_hi), 17); lr = np.linspace(lsig3, math.log(0.25), 25); lb = np.linspace(lbmin, lbmax, 25)
    A = np.array([(a, b, c, d, e, f) for a, b in zip(lm[:-1], lm[1:]) for c, d in zip(lr[:-1], lr[1:]) for e, f in zip(lb[:-1], lb[1:])])
    tails = [(-30.0, math.log(mu_lo), lsig3, math.log(0.25), lbmin, lbmax), (math.log(mu_hi), math.log(mu_hi) + 60, lsig3, math.log(0.25), lbmin, lbmax)]
    stack = [np.vstack([A, np.array(tails)])]; total = 0; minside = 1e9; bad_low_b = math.inf
    while stack:
        X = stack.pop()
        if len(X) > chunk: stack.extend(np.array_split(X, int(math.ceil(len(X) / chunk)))); continue
        total += len(X)
        if total > max_cells:
            return {'ell_a': ell_a, 'ell_b': ell_b, 'D': Dd, 'closed': False, 'cells': total, 'reason': 'cell budget',
                    'last_bad_log_b_over_sigma2': [float(X[:, 4].min() + 2 * ell_a), float(X[:, 5].max() + 2 * ell_a)],
                    'last_bad_log_r': [float(X[:, 2].min()), float(X[:, 3].max())], 'last_bad_log_mu': [float(X[:, 0].min()), float(X[:, 1].max())]}
        with np.errstate(all='ignore'):
            ok = cell_ok(X[:, 0], X[:, 1], X[:, 2], X[:, 3], X[:, 4], X[:, 5], ell_a, ell_b, Dd)
        bad = X[~ok]
        if not len(bad): continue
        w = np.stack([bad[:, 1] - bad[:, 0], (bad[:, 3] - bad[:, 2]) * 0.7, (bad[:, 5] - bad[:, 4]) * 0.7], axis=1)
        k = np.argmax(w, axis=1); ms = float(np.min(np.max(w, axis=1)))
        if ms < 1e-6:
            j = int(np.argmin(np.max(w, axis=1)))
            return {'ell_a': ell_a, 'ell_b': ell_b, 'D': Dd, 'closed': False, 'cells': total, 'reason': 'cells too small',
                    'stuck_cell': bad[j].tolist(), 'stuck_log_b_over_sigma2': float(bad[j, 4] + 2 * ell_a)}
        new = []
        for jj in range(3):
            sel = bad[k == jj]
            if not len(sel): continue
            mid = (sel[:, 2 * jj] + sel[:, 2 * jj + 1]) / 2
            L = sel.copy(); R = sel.copy(); L[:, 2 * jj + 1] = mid; R[:, 2 * jj] = mid; new += [L, R]
        stack.append(np.vstack(new))
    return {'ell_a': ell_a, 'ell_b': ell_b, 'D': Dd, 'lbmin_rel': lbmin_rel, 'tau_min': TAU_MIN, 'closed': True, 'cells_evaluated': total}


if __name__ == '__main__':
    ef, et, st = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]); Dd = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0
    lbr = float(sys.argv[5]) if len(sys.argv) > 5 else math.log(1e8)
    tmin = float(sys.argv[6]) if len(sys.argv) > 6 else 20.0
    set_tau_min(tmin)
    have = set()
    if os.path.exists(FN):
        for l in open(FN):
            if l.strip(): d = json.loads(l); have.add((round(d['ell_a'], 6), d['D'], round(d['lbmin_rel'], 6), d['tau_min']))
    ea = ef
    while ea < et - 1e-9:
        eb = min(et, ea + st)
        if (round(ea, 6), Dd, round(lbr, 6), tmin) not in have:
            t0 = time.time(); d = certify(ea, eb, Dd, lbr, max_cells=int(os.environ.get('CI_MAXCELLS', '400000000'))); d['lbmin_rel'] = lbr; d['tau_min'] = tmin; d['seconds'] = round(time.time() - t0, 1)
            open(FN, 'a').write(json.dumps(d) + '\n')
            print('ell [%.2f, %.2f] D %g tau >= %g b >= sigma^2 e^%.2f: %s cells %s %.0fs %s' % (ea, eb, Dd, tmin, lbr, 'CLOSED' if d['closed'] else 'OPEN (' + d.get('reason', '') + ')',
                  d.get('cells_evaluated', d.get('cells')), d['seconds'], d.get('stuck_log_b_over_sigma2', '')), flush=True)
        ea = eb
