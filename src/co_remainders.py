"""Stage CO 2.1-2.2: the catalogue of remainders of the LOCAL route (Lemma R_theta with the F term of Lemma CN-2), each bounded from above, and the
budget D measured at the minima of the local envelope.  Python 3.11, decimal arithmetic (ca_dec; sigma down to e^{-961} without underflow).

NOTATION (fixed w; E = e^{-w^2/2s^2}, h = 1 + E R; local average E_rho = theta_tau average at the focus; the line variable s ~ N(0, tau)):
  E0 = E_rho[R^2], E1 = E_rho[R'^2], mu = x E1/E0, P = 2 x E1, kappa = E_rho[(s^2/tau) R'^2]/E1, R_theta = E_rho[R] (|R_theta| <= eta),
  y = (1-b)^2 sigma^2/(2b), e = e^{-y/2}, beta^2 = tau y, k_E = int Gamma E = sqrt(r/(1+r)), se_E = s^2/(1+r), se2 = s^2/(2+r), g_s = sqrt(se2/L^2),
  W_theta = x g_s E1 = P g_s/2, D0_theta = E0 g_s se2/s^4, Z = 1 + k_E R_theta (normalisation; it multiplies Num and den alike, so the ratio is
  (Z Num)/(Z den) and only the b-pieces see Z), Psi_0 = 1 - (1-y) e e^{-L^2/4}, w = 3y/2 + L^2/4 (the linearised Psi of the x-free function).
THE LOCAL IDENTITY (exact: paper F eq. BLid with G[g] = E_rho[g], A[g] = E_rho[e^{i beta s/tau} g] (Lemma CN-2 (i)); for fixed w):
  sigma^2 int P s_v^2 phi dv = (B) + (F) + (X) + (L),
  (B) = b (E_rho[h] - (1-y) Re A[h e^{i psi}]),  (F) = sigma^2 (G[F] - Re A[F e^{i psi}]),  (X) = sigma^2 (1-b) Re A[h sigma^2 G e^{i psi}],
  (L) = sigma^2 Re A[D'] - sigma^2 G[h''],   and  Z Num = 2 sigma^2 (1/L^2 + Lambda_1 + J_theta,loc) - (sigma^2/2) Z E_p[s_w^2 phi] + (1/2) int Gamma [(B)+(F)+(X)+(L)],
  Z den = b Z + (W_theta,loc + T_w) + sigma^2 (1/L^2 + Lambda_1 + J_theta,loc),  T_w = x int Gamma sum n^2 q_n e_n cos phi_n (|T_w| <= 0.184 b k_E, L2),
  Lambda_1 = R_theta k_E (2+r)/((1+r) L^2) (L1), W_theta,loc in W_theta/(1 +- eta), J_theta,loc in D0_theta/(1 +- eta) (the main sandwiches).
MAIN PART (the x-free function of cn_envloc, sandwich factors 1):  Num_main = 2 sigma^2/L^2 + 2 sigma^2 D0_theta + b w/2 + W_theta se2/8 + W_theta kappa sigma^2/(8b),
  den_main = b + W_theta + sigma^2/L^2 + sigma^2 D0_theta.
REMAINDERS (derivation: report CO section 2; each N_i is a loss of (1/2)(...) or of 2 A_w, each Dn_i an excess of Z den):
  N1  [proof] L1 in 2 sigma^2 J_ww:                 2 sigma^2 eta k_E (2+r)/((1+r) L^2)
  N2  [proof] A_w, v'^2 term (ring-weighted, I8):   (sigma^4/(4b)) (1-b)^2 J2_up,  J2_up = (1/L^2)(1 + eta k_E r/(1+r)) + 2 eta k_E se_E/(L^2 s^2) + eta^2 g_s se2/(s^4 (1-eta))
  N3  [proof] A_w, w^2 term:                        (sigma^2/2) M4_up,  M4_up = 3(1 + eta k_E r^2/(1+r)^2) + 6 eta k_E r/(1+r)^2 + 3 se2 D0_theta/(1-eta)
  N4  [proof] A_w, chi^2 term:                      sigma^2 (3/2) [(xi^2 + zeta^2) J_up + (sigma^4/L^4) M4_up],  J_up = (1/L^2)(1 + eta k_E r/(1+r)) + 2 eta k_E se_E/(L^2 s^2) + D0_theta/(1-eta)
  N5  [proof] b-term, Z on b Psi_0 (I1):            (b/2) eta k_E Psi_0
  N6  [proof] b-term, the cos(beta s/tau) part (I3, I7, I10):  (b/2) eta k_E [(1-y)(e^{-se_E/4} - e e^{-L^2/4}) + y/2]
  N7  [proof] b-term, chi (|chi| <= chibar):         (b/2)(1+eta)(1-y) r_C
  N8  [proof] F term, chibar_e (F_rem of Lemma R with the local P):  (W_theta/(2(1-eta))) [(c0 sqrt(2 se2/pi) + c1 se2)/sqrt2 + c0^2 + c1^2 se2]
  N9  [proof] F term, first order of CN-2:          (W_theta/2) e^{-se2/4} xi_e sqrt(y kappa),  xi_e = xi/sqrt2
  N10 [proof] F term, second order of CN-2 (I9):    (1/2) 0.51 y^2 g_s U2 (1 + sqrt(pi/tau))      [x E1^{(2tau)} <= U2 theta_{2tau}(0)]
  N11 [proof] cross term, -sigma^2 (1-b)(W_theta,loc + T_w):  (sigma^2/2)(1-b)[W_theta/(1-eta) + 0.184 b k_E]
  N12 [proof] cross term, the local failure of the parity identity <h G> = -<F> (I5):  (sigma^3/2)(1-b) sqrt(x E1) k_E sqrt(2 se_E/pi)/(s^2 (1-eta))
  N13 [proof] cross term, gamma-tilde Phi_theta:    (sigma^2/2)(1-b)(1+eta) [int Gamma gamma-tilde (1 - e cos(w/sqrt2)) + int Gamma gamma-tilde chibar]
       (gamma-tilde(w) = the pointwise bound of |sigma^2 G| WITH its envelope factors E, (|w|/s^2) E; gamma_P of paper F is its sup over w)
  N14 [proof] last line, D - h' e^{iw/sqrt2} (I6):  (1/2) sqrt(1+y) [sqrt(P b) r_CE + sigma sqrt(b) (1+eta) int Gamma gamma-tilde],  r_CE = k_E (xi + zeta + sigma^2 sqrt(2 se_E/pi)/L^2)/sqrt2
       (A[d_v g] = (k'/tau) E_rho[e^{i beta s/tau}(s - i beta) g] and Cauchy-Schwarz in s: the local form of a_6; nu_1 -> sqrt((1+y)/tau), N_1 -> sqrt(E1))
  N15 [proof] last line, t_w (1 - cos(w/sqrt2)):   (1/2) 0.184 b k_E (1 - e^{-se_E/4})
  N16 [proof] last line, delta_2 (I4; a_3, a_8, a_10 of Lemma R in local form):  (1/2)(y/2) b Qbar k_E M2,  M2 = max_u u e^{-u/2}(1 + u cosh(beta sqrt(u/tau)))
  N17 [proof] linearisations of the x-free function: b (w - Psi_0)/2 + W_theta (se2/4 - (1 - e^{-se2/4}))/2 + (W_theta kappa/2)[(sigma^2/(4b) - (y/2) e^{-se2/4})]
  Dn1 [proof] Z on b in den:  b eta k_E;   Dn2 [proof] T_w in den:  0.184 b k_E;   Dn3 [proof] L1 in den:  sigma^2 eta k_E (2+r)/((1+r) L^2)
  Vanished (local forms of Lemma R's nu-terms): eps_1, z, 1 +- nu_0, omega_2 (exact local averages), a_1, a_2, a_5, a_9 (|theta_tau - 1| <= nu_0 is not used),
  a_4 (replaced by N9 + N10, Lemma CN-2), a_6 -> N14, a_3 + a_8 + a_10 -> N16 (+ N6), the (1/2) b Qbar (nu_0 + tau nu_2) Phi term -> N15 (+ N5).
D := ell [ (1 + Dden/ell)/(1 - Dnum/ell) - 1 ],  Dnum = ell sum N_i/Num_main,  Dden = ell sum Dn_i/den_main.
-> results/co_remainders.json
"""
import json, math, sys
from ca_dec import D, exp, log, sqrt, expm1, log1p, PI, E as EE
import numpy as np

QB = D('0.25'); SQ2 = sqrt(2)


def window_consts(ell):
    U1 = QB * sqrt(2 * D(ell) + D('3.1')) + D('1e-6'); U2 = QB * (2 * D(ell) + D('3.1')) + D('1e-4')
    gP = (U2 * (1 + QB) + U1 ** 2 + U1 * (1 + 2 * QB) / sqrt(EE)) / (1 - QB) ** 2
    cchi = (U1 / (1 - QB) + sqrt(2 / PI) + QB * exp(D('-0.5')) / (1 - QB)) / SQ2
    return U1, U2, gP, cchi


def M2_const(beta, tau):
    """max over u >= 0 of u e^{-u/2} (1 + u cosh(beta sqrt(u/tau))) (float; beta/sqrt(tau) = sqrt(y) tiny in the class)."""
    u = np.linspace(0, 400, 400001); c = float(beta / sqrt(tau))
    return float(np.max(u * np.exp(-u / 2) * (1 + u * np.cosh(c * np.sqrt(u))))) * (1 + 1e-6)


def catalogue(ell, mu, b, r, kap, P, x, A, eta):
    """all quantities at one state (decimal).  A = sigma^2/s^2.  Returns dict with main parts, remainders, D."""
    ell, mu, b, r, kap, P, x, A, eta = (D(t) for t in (ell, mu, b, r, kap, P, x, A, eta))
    sig = exp(-ell); s2 = sig ** 2 / A; L2 = s2 / r; s = sqrt(s2); L = sqrt(L2)
    se2 = s2 / (2 + r); gs = sqrt(se2 / L2); kE = sqrt(r / (1 + r)); seE = s2 / (1 + r)
    tau = x / b; y = (1 - b) ** 2 * sig ** 2 / (2 * b); e = exp(-y / 2); beta = sqrt(tau * y)
    E1 = P / (2 * x); E0 = P / (2 * mu)
    Wt = x * gs * E1; D0 = E0 * gs * se2 / s2 ** 2
    U1, U2, gP, cchi = window_consts(ell)
    xi = sig * U1 / (1 - QB); zeta = sig ** 2 * eta * exp(D('-0.5')) / (s * (1 - eta)); xie = xi / SQ2
    c0 = (xi + zeta) / SQ2; c1 = sig ** 2 / (SQ2 * L2)
    rC = (xi + zeta + sig ** 2 * sqrt(2 / PI) / L) / SQ2; rCE = kE * (xi + zeta + sig ** 2 * sqrt(2 * seE / PI) / L2) / SQ2
    eL = exp(-L2 / 4); eSE = exp(-seE / 4); ese2 = exp(-se2 / 4)
    Psi0 = -expm1(log1p(-y) - y / 2 - L2 / 4); w = 3 * y / 2 + L2 / 4
    # main part (x-free function, sandwich factors 1)
    Nm = {'2sig^2/L^2': 2 * sig ** 2 / L2, '2sig^2 D0': 2 * sig ** 2 * D0, 'b w/2': b * w / 2, 'W se2/8': Wt * se2 / 8, 'W kappa sig^2/(8b)': Wt * kap * sig ** 2 / (8 * b)}
    Dm = {'b': b, 'W': Wt, 'sig^2/L^2': sig ** 2 / L2, 'sig^2 D0': sig ** 2 * D0}
    Num_main = sum(Nm.values()); den_main = sum(Dm.values())
    # remainders
    J2 = (1 / L2) * (1 + eta * kE * r / (1 + r)) + 2 * eta * kE * seE / (L2 * s2) + eta ** 2 * gs * se2 / (s2 ** 2 * (1 - eta))
    Jup = (1 / L2) * (1 + eta * kE * r / (1 + r)) + 2 * eta * kE * seE / (L2 * s2) + D0 / (1 - eta)
    M4 = 3 * (1 + eta * kE * r ** 2 / (1 + r) ** 2) + 6 * eta * kE * r / (1 + r) ** 2 + 3 * se2 * D0 / (1 - eta)
    M2 = D(M2_const(beta, tau))
    Phi_up = (1 - e * eL) + eta * kE * (1 - e * eSE) + (1 + eta) * rC
    N = {}
    N['N1 L1'] = 2 * sig ** 2 * eta * kE * (2 + r) / ((1 + r) * L2)
    N['N2 A_w v2'] = sig ** 4 / (4 * b) * (1 - b) ** 2 * J2
    N['N3 A_w w2'] = sig ** 2 / 2 * M4
    N['N4 A_w chi2'] = sig ** 2 * D('1.5') * ((xi ** 2 + zeta ** 2) * Jup + sig ** 4 / L2 ** 2 * M4)
    N['N5 b Z'] = b / 2 * eta * kE * Psi0
    N['N6 b cos'] = b / 2 * eta * kE * ((1 - y) * (eSE - e * eL) + y / 2)
    N['N7 b chi'] = b / 2 * (1 + eta) * (1 - y) * rC
    N['N8 F chibar_e'] = Wt / (2 * (1 - eta)) * ((c0 * sqrt(2 * se2 / PI) + c1 * se2) / SQ2 + c0 ** 2 + c1 ** 2 * se2)
    N['N9 F first'] = Wt / 2 * ese2 * xie * sqrt(y * kap)
    N['N10 F second'] = D('0.5') * D('0.51') * y ** 2 * gs * U2 * (1 + sqrt(PI / tau))
    N['N11 X sig^2 W'] = sig ** 2 / 2 * (1 - b) * (Wt / (1 - eta) + D('0.184') * b * kE)
    N['N12 X parity'] = sig ** 3 / 2 * (1 - b) * sqrt(x * E1) * kE * sqrt(2 * seE / PI) / (s2 * (1 - eta))
    # gamma-tilde(w): the pointwise bound of |sigma^2 G| keeping its envelope factors (every term of sigma^2 G carries E or E_w = -(w/s^2) E):
    #   |sigma^2 G| <= x (E N2 + E^2 N1^2/(1-eta))/(1-eta) + sigma sqrt(x) (|w|/s^2) E N1 (1 + eta E/(1-eta))/(1-eta);  int Gamma of it (x N2 <= U2, sqrt x N1 <= U1):
    Ggt = U2 * kE + U1 ** 2 * gs / (1 - eta) ** 2 + sig * U1 * kE * sqrt(2 * seE / PI) * (1 + eta / (1 - eta)) / (s2 * (1 - eta))
    # int Gamma gamma-tilde (1 - e cos(w/sqrt2)):  E -> k_E (1 - e e^{-se_E/4}),  E^2 -> g_s (1 - e e^{-se2/4}),  (|w|/s^2) E -> (1-e) k_E sqrt(2 se_E/pi)/s^2 + e k_E sqrt(2/pi) se_E^{3/2}/(2 s^2)
    GgtPhi = (U2 * kE * (1 - e * eSE) + U1 ** 2 * gs * (1 - e * ese2) / (1 - eta) ** 2
              + sig * U1 * (1 + eta / (1 - eta)) / (s2 * (1 - eta)) * ((1 - e) * kE * sqrt(2 * seE / PI) + e * kE * sqrt(2 / PI) * seE * sqrt(seE) / 2))
    GgtChi = Ggt * (xi + zeta) / SQ2 + sig ** 2 * gP * sqrt(2 / PI) / (SQ2 * L)                       # int Gamma gamma-tilde chibar (|w| part with gamma_P)
    N['N13 X gamma_P'] = sig ** 2 / 2 * (1 - b) * (1 + eta) * (GgtPhi + GgtChi)
    # (L1) = sigma^2 Re A[d_v(D - h' e^{iw/sqrt2})], A[d_v g] = (k'/tau) E_rho[e^{i beta s/tau}(s - i beta) g]; Cauchy-Schwarz in s, then
    # |D - h' e^{iw/sqrt2}|^2 <= 2 k'^2 E^2 R'^2 chibar^2 + 2 h^2 |psi'|^2, |psi'| = |sigma^2 G|/sqrt2 <= gamma-tilde/sqrt2:
    N['N14 L D-term'] = D('0.5') * sqrt(1 + y) * (sqrt(P * b) * rCE + sig * sqrt(b) * (1 + eta) * Ggt)
    N['N15 L t_w'] = D('0.5') * D('0.184') * b * kE * (1 - eSE)
    N['N16 L delta2'] = D('0.5') * (y / 2) * b * QB * kE * M2
    N['N17 linearisation'] = b * (w - Psi0) / 2 + Wt * (se2 / 4 - (-expm1(-se2 / 4))) / 2 + Wt * kap / 2 * (sig ** 2 / (4 * b) - (y / 2) * ese2)
    Dn = {'Dn1 b Z': b * eta * kE, 'Dn2 T_w': D('0.184') * b * kE, 'Dn3 L1': sig ** 2 * eta * kE * (2 + r) / ((1 + r) * L2)}
    Dnum = ell * sum(N.values()) / Num_main; Dden = ell * sum(Dn.values()) / den_main
    Dtot = ell * ((1 + Dden / ell) / (1 - Dnum / ell) - 1) if Dnum < ell else D('Infinity')
    env_main = Num_main / (sig * den_main) * sqrt(ell)
    return {'inputs': {'ell': ell, 'mu': mu, 'b': b, 'r': r, 'kappa': kap, 'P': P, 'x': x, 'A': A, 'eta': eta, 'tau': tau, 's2': s2, 'L2': L2, 'y': y},
            'Num_main': Nm, 'den_main': Dm, 'N': N, 'Dn': Dn, 'Dnum': Dnum, 'Dden': Dden, 'D': Dtot, 'env_main': env_main,
            'rel_N': {k: ell * v / Num_main for k, v in N.items()}, 'rel_Dn': {k: ell * v / den_main for k, v in Dn.items()}}


def A_opt(mu, b, r, kap, P, ell):
    """the A of the x-free function's minimum: alpha A + beta sigma^2/A -> A* = sigma sqrt(beta/alpha) (cn_envloc).  Decimal (b, r may be e^{-900})."""
    mu, b, r, P = D(mu), D(b), D(r), D(P)
    G = P * sqrt(r / (2 + r)); kL = 2 * QB * sqrt(r / (1 + r))
    alpha = r * 2 * (1 - kL) + G / (mu * (2 + r)); beta = b / (8 * r) + G / (16 * (2 + r))
    return exp(-D(ell)) * sqrt(beta / alpha)


def brightness(mu, b, kap, kind):
    import cn_search as S, cn_lemmas as CL
    if kind == 'Psi2': return min(2 * mu * 0.0625, CL.Psi2(mu, b, kap))
    if kind == 'Psi2opt':
        import co_envmin as EM
        return min(2 * mu * 0.0625, EM.Psi2_opt(mu, b, kap))
    if kind == 'Psi': return min(2 * mu * 0.0625, S.psi_raw(mu, b))
    raise ValueError


def measure(ell, state, eta, kind='Psi2', nx=25):
    """D over x in [5e11 sigma^2, min(0.003, 20 b)] (tau < 20) and over all x <= 0.003 (the local route needs no tau restriction)."""
    mu, kap = state['mu'], state['kappa']; lb = state['log_b_over_sigma2'] - 2 * ell; lr = state['log_r']
    b = exp(D(lb)); r = exp(D(lr))
    P = brightness(mu, math.exp(lb) if lb > -700 else 0.0, kap, kind); A = A_opt(mu, b, r, kap, P, ell)
    lxlo = math.log(5e11) - 2 * ell; lxhi = math.log(0.003)
    rows = []
    for lx in np.linspace(lxlo, lxhi, nx):
        x = exp(D(lx)); c = catalogue(ell, mu, b, r, kap, P, x, A, eta)
        rows.append({'log_x': float(lx), 'log_tau': float(log(c['inputs']['tau'])), 'tau': float(c['inputs']['tau']), 'D': float(c['D']), 'Dnum': float(c['Dnum']), 'Dden': float(c['Dden']), 'env_main': float(c['env_main']),
                     'rel_N': {k: float(v) for k, v in c['rel_N'].items()}, 'rel_Dn': {k: float(v) for k, v in c['rel_Dn'].items()}})
    w_all = max(rows, key=lambda t: t['D']); w20 = max([t for t in rows if t['log_tau'] < math.log(20)] or rows, key=lambda t: t['D'])
    return {'state': {'mu': mu, 'log_b': lb, 'log_r': lr, 'kappa': kap, 'P': P, 'log_A': float(log(A)), 'eta': eta, 'kind': kind}, 'rows': rows, 'worst_all_x': w_all, 'worst_tau_lt_20': w20}


if __name__ == '__main__':
    env = json.load(open('results/cn_envloc.json'))
    try: env2 = json.load(open('results/co_envpsi2.json'))
    except FileNotFoundError: env2 = {}
    res = {}
    for key, row in env.items():
        ell = float(key) if key != '55.26' else 55.26204223
        ell = {'55.26': 55.26204223, '110.52': 110.52408446, '221.05': 221.04816893, '455.26': 455.26, '690.78': 690.77552790, '960.8': 960.8}[key]
        eta_i = 1 / ell
        out = {'T0': 1 + 2.1 / ell + 2.5e-4, 'T10': row['T10']}
        states = {'argmin_C=1.41421': row['argmin_C=1.41421'], 'argmin_C=50': row['argmin_C=50']}
        if key in env2: states['argmin_Psi2'] = env2[key]['argmin_Psi2']
        try: states['argmin_envmin'] = json.load(open('results/co_envmin.json'))[key]['argmin']        # the multi-start minimum (co_envmin.py): the reference state
        except (FileNotFoundError, KeyError): pass
        try: states['argmin_envmin_opt'] = json.load(open('results/co_envmin_opt.json'))[key]['argmin']    # the same with (K, L) of Psi_2 optimised continuously
        except (FileNotFoundError, KeyError): pass
        for sname, st in states.items():
            for eta_name, eta in (('eta=1/ell', eta_i), ('eta=0.111', 0.111), ('eta=0.5', 0.5)):
                m = measure(ell, st, eta, kind='Psi2opt' if 'opt' in sname else 'Psi2' if 'Psi2' in sname or 'C=50' in sname or 'envmin' in sname else 'Psi')
                out['%s %s' % (sname, eta_name)] = m
                print(key, sname, eta_name, 'D(worst, tau<20) %.3e  D(worst, all x) %.3e  env_main %.4f  top: %s' % (
                    m['worst_tau_lt_20']['D'], m['worst_all_x']['D'], m['worst_tau_lt_20']['env_main'],
                    sorted(m['worst_tau_lt_20']['rel_N'].items(), key=lambda t: -t[1])[:3]), flush=True)
        res[key] = out
    json.dump(res, open('results/co_remainders.json', 'w'), indent=1, default=float)
