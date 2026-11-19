"""Paper H v1.3: the raw values behind the numbers printed in Sec. IX A (the same selection as verification/hp_numbers.py, lines on
co_envmin_opt, co_envmin, co_remainders, co_rtheta, co_envpsi2), read from a results/ folder.  Python 3.11 grammar.
  python src/dv_section9.py RESULTS_DIR OUT.json   writes the values (used once on the v1.2.0 results -> results/v12_section9_values.json);
  values(dir) is used by verification/run_v13.py to compare v1.3 with v1.2 in the printed digits."""
import json, os, sys


def values(R):
    L = lambda f: json.load(open(os.path.join(R, f)))
    EO = L('co_envmin_opt.json'); EM = L('co_envmin.json'); CO = L('co_remainders.json'); RT = L('co_rtheta.json'); E2 = L('co_envpsi2.json')
    neg = [x for x in RT['states'] if x['lemma_R_period'] < 0]
    v = {}
    for k in EO: v['envmin_opt %s' % k] = (EO[k]['min'], '%.4f'); v['T0 %s' % k] = (EO[k]['T0'], '%.4f')
    for k in EM: v['envmin grid %s' % k] = (EM[k]['min'], '%.4f')
    v['D at 55.26 (printed as ceil to 2 digits)'] = (CO['55.26']['argmin_envmin_opt eta=1/ell']['worst_all_x']['D'], 'ceil2')
    v['rtheta states'] = (len(RT['states']), '%d')
    v['rtheta min ratio (Lemma R negative)'] = (min(x['ratio_theta'] for x in neg), 'pct')
    v['rtheta max ratio (Lemma R negative)'] = (max(x['ratio_theta'] for x in neg), 'pct')
    v['envpsi2 Cinf 55.26'] = (E2['55.26']['Cinf_kmin=1e-12'], '%.2f'); v['envpsi2 Psi2 55.26'] = (E2['55.26']['Psi2'], '%.4f')
    return v


def printed(val, fmt):
    import math
    if fmt == 'pct': return '%d' % round(100 * val)
    if fmt == 'ceil2':
        e = math.floor(math.log10(val)); m = math.ceil(val / 10 ** e * 10 - 1e-9) / 10; return '%.1fe%d' % (m, e)
    return fmt % val


if __name__ == '__main__':
    v = values(sys.argv[1]); json.dump({k: [a, f] for k, (a, f) in v.items()}, open(sys.argv[2], 'w'), indent=1)
    for k, (a, f) in v.items(): print(k, a, printed(a, f))
