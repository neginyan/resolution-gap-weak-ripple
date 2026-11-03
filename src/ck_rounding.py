"""Stage CK 2.2 (stopped here): is the acceptance margin of the cells above the floating-point rounding?  Python 3.11.  [numerical]
(1) Main-inequality cells (cc_cells.cell_ok, used by certificate (i)): acceptance is lpos >= lneg + 1e-12 on LOGARITHMS whose size reaches 3 ell
    (r, b down to sigma^3) and more in the tails.  We compare the margin 1e-12 with one ulp of the compared numbers, and recompute lpos - lneg in
    long double (64-bit mantissa) on random cells to measure the actual double-precision error.
(2) Link cells (ce_link, certificate (ii)): acceptance is (minimum)/(magnitude) > 1e-10 after netting equal sigma-powers; we record the size of the
    logarithms there (per-term relative error ~ (#operations) * u * |log|).
-> results/ck_rounding.json"""
import json, math
import numpy as np
import cc_cells as CC, cc_core as C

LD = np.longdouble
def cc_margin(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, delta, T):
    """lpos - lneg of cc_cells.cell_ok at both G end points, in dtype T (same formulas)."""
    f = lambda v: np.asarray(v, dtype=T)
    lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b = map(f, (lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b))
    LN2 = np.log(T(2)); gam = (1 + T(delta)) / np.sqrt(T(ell_a)); lgam = np.log(gam)
    lsig_lo = T(-ell_b); sig_hi = np.exp(T(-ell_a)); r_hi = np.exp(lr_b)
    la0 = lr_a + np.log(2 - gam * sig_hi); la1 = np.log(1 - gam * sig_hi / 2) - lmu_b - np.log(2 + r_hi)
    lk0 = lb_a - LN2 * 3 - lr_b; lk1 = -np.log(T(16)) - np.log(2 + r_hi)
    lPi = np.array([C.log_Pistar(max(float(m), 2.0)) for m in np.exp(lmu_a.astype(float))], dtype=T)
    lPi = np.where(lmu_a >= np.log(T(2)), lPi, lmu_b + 2 * np.log(T(C.QB))) + T(1e-6)
    lG = lPi + T(0.5) * (lr_b - np.log(2 + r_hi))
    lae = lambda a, b: np.maximum(a, b) + np.log1p(np.exp(-np.abs(a - b)))
    out = []
    lT1 = LN2 + T(0.5) * (lae(la0, la1 + lG) + lae(lk0, lk1 + lG)); lT3 = lG + lsig_lo - lb_b - np.log(T(16))
    lpos = lae(lae(lT1, lsig_lo + np.log(T(3) / 8) + 0 * lT1), lT3); lneg = lae(lgam + lb_b, lgam + lG - LN2)
    out.append(lpos - lneg)
    lT1 = LN2 + T(0.5) * (la0 + lk0); lpos = lae(lT1, lsig_lo + np.log(T(3) / 8) + 0 * lT1); lneg = lgam + lb_b
    out.append(lpos - lneg)
    return out, np.maximum(np.abs(lpos), np.abs(lneg))

res = {'cc_cells': [], 'ce_link': []}
rng = np.random.default_rng(3)
for ell in (55.26, 960.8, 6897.8):
    n = 20000; e_a, e_b = ell, ell + 10
    lmu = rng.uniform(math.log(1e-4), math.log(4 * e_b), n); dm = rng.uniform(1e-6, 0.3, n)
    lr = rng.uniform(-3 * e_b, math.log(0.25), n); dr = rng.uniform(1e-6, 1.0, n)
    lb = rng.uniform(-3 * e_b, 0.0, n); db = rng.uniform(1e-6, 1.0, n)
    args = (lmu, lmu + dm, lr - dr, lr, lb - db, lb, e_a, e_b, 10.0 / e_a + 2.1 / e_a + 2.5e-4)
    (d1, d0), mag = cc_margin(*args, np.float64)
    (D1, D0), magL = cc_margin(*args, LD)
    err = float(max(np.max(np.abs(d1 - D1.astype(float))), np.max(np.abs(d0 - D0.astype(float)))))
    res['cc_cells'].append({'ell': ell, 'max_abs_log_compared': float(np.max(mag)), 'ulp_of_largest_log': float(np.spacing(np.max(mag))),
                            'margin_used': 1e-12, 'max_double_minus_longdouble_in_margin': err})
    # ce_link: the logarithms are of size p*ell (p <= 3) and lr >= -3 ell; per-term error ~ 10 ops * u * |log|
    L = 3 * (ell + 10) + 700
    res['ce_link'].append({'ell': ell, 'max_abs_log': L, 'per_term_relative_error_est': 10 * 1.11e-16 * L, 'threshold_used': 1e-10})
json.dump(res, open('results/ck_rounding.json', 'w'), indent=1)
print(json.dumps(res, indent=1))
