"""Stage CL: the acceptance rule of the main-inequality cells (cc_cells.cell_ok, used by cc_cells.py and cd_cells.py) made to COVER ROUNDING.
Python 3.11.

ROUNDING MODEL (as cl_link.py): IEEE 754 binary64, round to nearest, u = 2^-53; log, exp, log1p within 2 ulp; cell ends exact.
DERIVED BOUND.  cell_ok compares two logarithms, lpos and lneg.  Each is obtained from the exact cell ends through a chain of at most 25 floating
operations (additions, products by 1/2, log, log1p, exp, max), each with error at most u times the magnitude of its operands (2u for the
elementary functions).  The operands are the cell ends (each used at most twice), ell_b, constants below 100, and lG = log G_cell, which can be
huge in the mu tails; but lG enters only inside a log-add-exp, where it either dominates (then |lG| <= |lpos| + |lneg| + B) or is dropped
exactly.  Hence, with B = 2(|lr_a| + |lr_b| + |lb_a| + |lb_b| + |lmu_a| + |lmu_b|) + ell_b + 100,
    |lpos^ - lpos| + |lneg^ - lneg| <= delta := 60 u (|lpos^| + |lneg^| + B)      (30 u per side, 25 operations plus the elementary functions).
RULE: accept iff lpos^ >= lneg^ + max(1e-12, S delta) with the safety factor S = 10.  (The old rule had the margin 1e-12 only.)"""
import math
import numpy as np
import cc_cells as CC, cc_core as C

U = 2.0 ** -53; SAFE = 10.0
MIN_SLACK = [math.inf]        # smallest (lpos - lneg)/(S delta) over accepted cells, recorded for the report
lae2 = CC.lae2; lPistar_vec = CC.lPistar_vec; LN2 = math.log(2.0)


def cell_ok(lmu_a, lmu_b, lr_a, lr_b, lb_a, lb_b, ell_a, ell_b, delta):
    gam = (1 + delta) / math.sqrt(ell_a); lgam = math.log(gam)
    lsig_lo = -ell_b; sig_hi = math.exp(-ell_a)
    r_hi = np.exp(lr_b)
    la0 = lr_a + math.log(2 - gam * sig_hi)
    la1 = math.log(1 - gam * sig_hi / 2) - lmu_b - np.log(2 + r_hi)
    lk0 = lb_a - LN2 * 3 - lr_b
    lk1 = -math.log(16) - np.log(2 + r_hi)
    lPi = np.where(lmu_a >= math.log(2.0), lPistar_vec(np.maximum(lmu_a, math.log(2.0))), lmu_b + 2 * math.log(C.QB)) + 1e-6
    lG = lPi + 0.5 * (lr_b - np.log(2 + r_hi))
    B = 2 * (np.abs(lr_a) + np.abs(lr_b) + np.abs(lb_a) + np.abs(lb_b) + np.abs(lmu_a) + np.abs(lmu_b)) + ell_b + 100
    out = np.ones(len(lmu_a), dtype=bool); sl = []
    for Gz in (False, True):
        if Gz:
            lT1 = LN2 + 0.5 * (lae2(la0, la1 + lG) + lae2(lk0, lk1 + lG))
            lT3 = lG + lsig_lo - lb_b - math.log(16)
            lpos = lae2(lae2(lT1, np.full_like(lT1, lsig_lo + math.log(3 / 8))), lT3)
            lneg = lae2(lgam + lb_b, lgam + lG - LN2)
        else:
            lT1 = LN2 + 0.5 * (la0 + lk0)
            lpos = lae2(lT1, np.full_like(lT1, lsig_lo + math.log(3 / 8)))
            lneg = lgam + lb_b
        dl = 60 * U * (np.abs(lpos) + np.abs(lneg) + B)
        need = np.maximum(1e-12, SAFE * dl); okk = lpos >= lneg + need
        sl.append(np.where(okk, (lpos - lneg) / need, np.inf))
        out &= okk
    if out.any(): MIN_SLACK[0] = min(MIN_SLACK[0], float(np.min(np.minimum(sl[0], sl[1])[out])))
    return out
