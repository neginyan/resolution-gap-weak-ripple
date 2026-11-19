"""Paper H v1.3 (stage DW): Lemma BE at the upper end of sigma, sigma = e^{-55.26} = 1.00204e-24 (ell = 55.26).  Python 3.11 grammar.
Up to v1.2 the hypothesis of Lemma BE was stated as sigma in [1e-3000, 1e-24]; Theorem H1 starts at ell = 55.26, and 24 ln 10 = 55.2620..., so the
states with ell in [55.26, 24 ln 10) (sigma in (1e-24, e^{-55.26}]) were left out of the stated hypothesis.  v1.3 states the lemma for
sigma in [1e-3000, e^{-55.26}].  Every place where the proof uses an upper bound on sigma is re-checked here at
  sigma_up = e^{-ell_b} (1 + 1e-40),  ell_b = the binary64 value of 55.26 (55.25999999999999801...; it is stored slightly below 55.26, so this sigma is
  slightly larger than e^{-55.26}),
in 60-digit decimal arithmetic, with the formulas of ch_bedge.check_constants / check_remainders (copied below with sigma_hi as a parameter):
  (1) the constants of (BE-3): y_min, 1 - e, rho_W <= 0.0801, the coefficients and kappa >= 1.1e-5 (both cases L^2 >= 1 and L^2 < 1);
  (2) the sigma^2 term of rho_W (App. D 5: "below 1e-40"): 48 (1 + Qb) sigma^2 (U_1^2 + 1/2)/(1 - Qb)^3, U_1 at ell = 6907.76;
  (3) the tail piece T_Lambda, Lambda = 1e-3/sigma^2, ratio 22.3 T_Lambda/sigma^2 (worst at the largest sigma): below 10^{-4e44};
  (4) the consequence g sqrt(ell) >= (1 - epsbar) kappa sqrt(ell)/sigma >= (1 - epsbar) 8e19 at ell = 55.26;
  (5) the remainders a_i: worst at the SMALLEST sigma (1e-3000), unchanged (they decrease as sigma grows; recomputed at sigma_up as well).
Consistency: the copy with sigma_hi = 1e-24 reproduces results/ch_bedge.json to the last printed digit.
Power (deliberately wrong versions, each must fail): sigma_hi = 1e-20 against "below 1e-40"; sigma = 1e-23 against "8e19"; Lambda = 1e-1/sigma^2
against rho_W <= 0.0801.
-> results/dw_bedge_low.json"""
import json, os, sys, time
from decimal import Decimal as D, getcontext
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ch_bedge as CH

getcontext().prec = 60; getcontext().Emin = -10 ** 15; getcontext().Emax = 10 ** 15
mp = CH.mp; QB = CH.QB
ELL_B = D(55.26)                                     # exact binary64 value of 55.26
SIG_UP = (-ELL_B).exp() * (1 + D('1e-40'))


def constants(sig_hi, lam_c=D('1e-3')):
    Qb = mp.mpf(QB); bmax = 10 ** 4 * sig_hi ** 2
    y_min = (1 - bmax) ** 2 / (2 * 10 ** 4); one_m_e = y_min / (2 + y_min)
    U1, U2, gP, cchi = CH.U_consts(mp.mpf('6907.76'))
    sig_term = sig_hi ** 2 * (U1 ** 2 + mp.mpf('0.5')) / (1 - Qb) ** 3
    K1s2 = lam_c / (1 - Qb) + sig_term
    rhoW = 48 * (1 + Qb) * K1s2 / (1 - bmax) ** 2
    cW2 = 1 - mp.mpf('4.01e4') * sig_hi ** 2 - rhoW * (1 + mp.mpf('1e-100'))
    cW = cW2 / 2 * one_m_e
    cJ = 1 - mp.mpf(1) / (8 * (1 - mp.mpf('1e-100')))
    cb_coef = (1 / (1 + mp.mpf('1e-100')) - 1 / (2 * (1 - mp.mpf('1e-100')))) / 2
    Phi_min = 1 - mp.exp(-mp.mpf(1) / 4) - sig_hi * cchi
    Psi_min = min(Phi_min, mp.mpf('0.5')); cb = cb_coef * Psi_min
    rem_allow = mp.mpf('1e-6')
    case1 = min(cJ - rem_allow, cW, cb)
    case2 = min((cJ - rem_allow) / (1 + mp.mpf(10 ** 4) / mp.mpf('0.6')), cW)
    kappa = min(case1, case2)
    return {'y_min': y_min, 'one_minus_e_min': one_m_e, 'rho_W': rhoW, 'rho_W_sigma_term': 48 * (1 + Qb) * sig_term / (1 - bmax) ** 2,
            'coef_W_in_Num': cW, 'Phi_min_L2_ge_1': Phi_min, 'kappa_case_L2_ge_1': case1, 'kappa_case_L2_lt_1': case2, 'kappa': kappa}


def tail_log10(sig):
    Lam = mp.mpf('1e-3') / sig ** 2
    lpoly = mp.log10(4 * (Lam ** 2 + Lam ** D('1.5') / (mp.sqrt(mp.mpf('5e11')) * sig) / mp.mpf('1.25')))
    return mp.log10(mp.mpf('22.3')) + lpoly - Lam / D(10).ln() - 2 * mp.log10(sig)


def run(save=True):
    t0 = time.time(); out = {}
    rec = json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'results', 'ch_bedge.json')))['BE_3_constants']
    c24 = constants(mp.mpf('1e-24'))
    keys = ['y_min', 'one_minus_e_min', 'rho_W', 'coef_W_in_Num', 'Phi_min_L2_ge_1', 'kappa_case_L2_ge_1', 'kappa_case_L2_lt_1', 'kappa']
    out['consistency: the copy at sigma_hi = 1e-24 equals results/ch_bedge.json'] = all(float(c24[k]) == rec[k] for k in keys)
    cu = constants(SIG_UP)
    out['sigma_up'] = '%.12e' % SIG_UP; out['ell_b (binary64 of 55.26)'] = str(ELL_B)
    out['(1) constants at sigma_up'] = {k: float(v) for k, v in cu.items()}
    out['(1) rho_W <= 0.0801'] = bool(cu['rho_W'] <= D('0.0801'))
    out['(1) kappa >= 1.1e-5'] = bool(cu['kappa'] >= D('1.1e-5'))
    out['(1) change of kappa from sigma = 1e-24 (relative)'] = float(abs(cu['kappa'] / c24['kappa'] - 1))
    out['(2) sigma^2 term of rho_W at sigma_up'] = float(cu['rho_W_sigma_term'])
    out['(2) below 1e-40'] = bool(cu['rho_W_sigma_term'] < D('1e-40'))
    tl = tail_log10(SIG_UP)
    out['(3) log10 of the tail ratio at sigma_up'] = float(tl)
    out['(3) below 10^(-4e44)'] = bool(tl < D('-4e44'))
    g = D('1.1e-5') * ELL_B.sqrt() / SIG_UP
    out['(4) kappa sqrt(ell)/sigma at ell = 55.26 (without 1 - epsbar)'] = float(g)
    out['(4) >= 8e19'] = bool(g >= D('8e19'))
    with_remainders = CH.check_remainders()
    out['(5) remainders (worst at sigma = 1e-3000; ch_bedge.check_remainders): largest log10 ratio'] = with_remainders['all_below_log10']
    out['(5) ok'] = bool(with_remainders['ok'])
    # POWER
    p = {}
    p['sigma_hi = 1e-20: the sigma^2 term of rho_W is NOT below 1e-40'] = bool(constants(mp.mpf('1e-20'))['rho_W_sigma_term'] >= D('1e-40'))
    p['sigma = 1e-23: kappa sqrt(55.26)/sigma is NOT >= 8e19'] = bool(D('1.1e-5') * ELL_B.sqrt() / D('1e-23') < D('8e19'))
    p['Lambda = 1e-1/sigma^2: rho_W is NOT <= 0.0801'] = bool(constants(SIG_UP, lam_c=D('1e-1'))['rho_W'] > D('0.0801'))
    out['power (each wrong version broken)'] = p
    out['ok'] = bool(out['consistency: the copy at sigma_hi = 1e-24 equals results/ch_bedge.json'] and out['(1) rho_W <= 0.0801'] and out['(1) kappa >= 1.1e-5']
                     and out['(2) below 1e-40'] and out['(3) below 10^(-4e44)'] and out['(4) >= 8e19'] and out['(5) ok'] and all(p.values()))
    out['seconds'] = time.time() - t0
    if save:
        json.dump(out, open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'results', 'dw_bedge_low.json'), 'w'), indent=1)
    return out


if __name__ == '__main__':
    print(json.dumps(run(), indent=1))
