"""Stage CL, task 2.1.3: do the tau cells (App. A, factor 0.99977) and the Lemma BE constants (ch_bedge.py) have an acceptance rule of the
same kind as certificates (i)-(iii), i.e. a double-precision comparison whose margin might be below the rounding error?

tau cells (cb_tau.tau_cells, tau_c = 20): the factor F(tau, beta) = sqrt(k1 k2)/c0 is evaluated in double at the left end of each of 400 cells;
the claimed bound is min F >= 0.99977.  Here we recompute F at every left end in 50-digit decimal (exact model of the same formula; exp is
decimal's correctly rounded exp) and compare: |F_double - F_decimal| must be far below the margin min F - 0.99977.
Derived bound for the double value: F is built from <= 80 exponentials (relative error <= 4u each, values <= 1e-4), one sum, eight
+,-,*,/ on O(1) numbers and one sqrt; Higham: |dF| <= (80*4u*nu0 + 12u + 4u) * F <= 20u.  The safety rule is margin >= 10 * 20u.

Lemma BE (ch_bedge.py): the constants (y_min, 1 - e, rho_W, coefficients, kappa) are computed in 60-digit decimal; the remainders are
compared in log10 at sizes ~ -1e7, i.e. with margins of millions of decades.  No double-precision acceptance rule of the cell type is used
for the stated constant kappa >= 1.1e-5; we record the margins.  The in-code double checks (formula vs finite differences,
envelope ratio 0.767, mean-square ratio 0.244, 0.830) have margins >= 0.17, i.e. > 1e14 u.
"""
import json, math
import numpy as np
from decimal import Decimal as D, getcontext

getcontext().prec = 50
U = 2.0 ** -53


def F_double(tau, beta, Qb=0.25):
    n = np.arange(1, 80); en = np.exp(-n * n * tau / 2 + n * beta)
    nu0 = 2 * en.sum(); z = Qb * nu0 / 2; eps1 = nu0 * (1 + 2 * Qb) / (1 - z)
    k1 = 1 - eps1; k2 = (1 - nu0) / (1 + z); c0 = (1 + nu0) / (1 - z)
    return math.sqrt(k1 * k2) / c0


def F_dec(tau, beta, Qb=D('0.25')):
    tau = D(repr(float(tau))); beta = D(repr(float(beta)))
    nu0 = 2 * sum((-(D(n) * n) * tau / 2 + n * beta).exp() for n in range(1, 80))
    z = Qb * nu0 / 2; eps1 = nu0 * (1 + 2 * Qb) / (1 - z)
    k1 = 1 - eps1; k2 = (1 - nu0) / (1 + z); c0 = (1 + nu0) / (1 - z)
    return (k1 * k2).sqrt() / c0


def main():
    out = {}
    edges = np.geomspace(20.0, 80.0, 401)
    worst_d = 9.0; worst_D = D(9); maxerr = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        fd = F_double(float(a), 1e-6 * float(b)); fD = F_dec(float(a), 1e-6 * float(b))
        worst_d = min(worst_d, float(fd)); worst_D = min(worst_D, fD)
        maxerr = max(maxerr, abs(D(repr(float(fd))) - fD) / fD)
    margin = float(worst_D - D('0.99977'))
    out['tau_cells'] = {'cells': 400, 'min_F_double': worst_d, 'min_F_decimal50': float(worst_D), 'claimed': 0.99977,
                        'margin': margin, 'max_rel_err_double_vs_decimal': float(maxerr), 'derived_bound_20u': 20 * U,
                        'margin_over_bound': margin / (20 * U), 'ok': bool(margin >= 10 * 20 * U and maxerr <= 20 * U)}
    be = json.load(open('results/ch_bedge.json'))
    c = be['BE_3_constants']; r = be['BE_3_remainders_log10']
    out['lemma_BE'] = {'arithmetic': '60-digit decimal (ported from mpmath)', 'kappa': c['kappa'], 'kappa_claimed': c['kappa_claimed'],
                       'kappa_rel_margin': c['kappa'] / c['kappa_claimed'] - 1, 'remainders_max_log10': r['all_below_log10'],
                       'double_checks': {'envelope_ratio': be['BE_a_formula_and_envelope']['max_ratio_to_envelope_bound'],
                                         'mean_square_ratio': be['BE_2_mean_square']['max_MG_over_bound'],
                                         'x2S4_ratio': be['BE_2_mean_square']['max_x2S4_over_Lambda_split']},
                       'rule_of_cell_type': False,
                       'ok': bool(be['all_ok'] and c['kappa'] / c['kappa_claimed'] - 1 > 1e-3)}
    out['all_ok'] = out['tau_cells']['ok'] and out['lemma_BE']['ok']
    json.dump(out, open('results/cl_tau_be.json', 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
