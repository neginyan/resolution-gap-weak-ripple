"""Paper H (draft): every number of the manuscript re-derived from results/ and checked against the TeX sources (the bw_numbers.py pattern of paper F).
Bounds are rounded outward (floor for lower bounds, ceil for upper bounds) where the text states a bound.  Python 3.11.
Usage: python hp_numbers.py [ROOT] [TEXDIR]  (ROOT holds results/, TEXDIR holds sections/)."""
import json, math, os, sys, glob

def dn(v, d): return math.floor(v * 10 ** d) / 10 ** d
def up(v, d): return math.ceil(v * 10 ** d) / 10 ** d
def sci(v, d=1):
    e = int(math.floor(math.log10(abs(v)))); m = v / 10 ** e
    return '%s%.*f\\times10^{%d}' % ('-' if v < 0 else '', d, abs(m), e)

def numbers(root):
    R = lambda f: os.path.join(root, 'results', f)
    jl = lambda f: [json.loads(l) for l in open(R(f)) if l.strip()]
    out = []
    add = lambda txt, val, src: out.append({'text': txt, 'value': val, 'source': src})
    env = {json.loads(l)['sig_exp']: json.loads(l)['gsl_env2'] for l in open(R('ca_envelope.jsonl')) if l.strip()}
    add('$g\\sqrt\\ell=%.4f$ at $\\sigma=10^{-24}$' % round(env[24], 4), env[24], 'ca_envelope (numerical)')
    for e in (48, 96, 300, 3000):
        add('$%.4f$ at $10^{-%d}$' % (round(env[e], 4), e), env[e], 'ca_envelope (numerical)')
    MAIN = jl('cl_main_cells.jsonl')          # certificate (i) with the rule of App. E (stage CL)
    CC = sorted([r for r in MAIN if r['cert'] == 'cc'], key=lambda r: r['ell_a'])
    assert all(r['closed'] for r in CC) and sum(r['audit_violations'] for r in CC) == 0
    add('$%d$ columns $\\ell\\in[%.3f,%.1f]$ close, with $%d$ audit points' % (len(CC), CC[0]['ell_a'], CC[-1]['ell_b'], sum(r['audit_points'] for r in CC)), len(CC), 'cc_cells (v1.1: the exact start 55.262)')
    D = sorted([r for r in MAIN if r['cert'] == 'cd'], key=lambda r: r['ell_a'])
    assert all(r['closed'] for r in D) and sum(r['audit_violations'] for r in D) == 0
    add('all $%d$ columns $\\ell\\in[%.2f,%.1f]$ close' % (len(D), D[0]['ell_a'], D[-1]['ell_b']), len(D), 'cd_cells')
    add('with $%d$ audit points' % sum(r['audit_points'] for r in D), sum(r['audit_points'] for r in D), 'cd_cells audit')
    add('up to $%s$ boxes' % sci(up(max(r['cells_evaluated'] for r in D) / 1e7, 1) * 1e7), max(r['cells_evaluated'] for r in D), 'cd_cells max')
    L = sorted(jl('cl_link_cells.jsonl'), key=lambda r: r['ell_a'])      # certificate (ii) with the rule of App. E (stage CL)
    assert all(r['closed'] for r in L) and all(abs(a['ell_b'] - b['ell_a']) < 1e-6 for a, b in zip(L, L[1:]))
    reach = L[-1]['ell_b']
    add('the $%d$ columns $\\ell\\in[%.2f,%.1f]$ close' % (len(L), L[0]['ell_a'], reach), reach, 'ce_link_cells (reach)')
    add('$\\ell\\in[55.26,%.1f]$' % reach, reach, 'Theorem H1 range')
    W = jl('cf_windows.jsonl'); assert len(W) == 17 and all(w['closed'] for w in W)
    add('All $%d$ windows close' % len(W), len(W), 'cf_windows')
    F = json.load(open(R('ce_a4.json')))['forced_tau']
    g = lambda e, k: [r for r in F if r['sig_exp'] == e and r['tau'] == 20 and r['state'] == 'optimum L2,s2'][0][k]
    for e in (48, 300):
        add('$%s$ & $%s$' % (sci(g(e, 'TA_sqrt_ell_old')), sci(g(e, 'TA_sqrt_ell_new'))), g(e, 'TA_sqrt_ell_new'), 'ce_a4 forced tau (numerical)')
    add('$T_B\\sqrt\\ell\\ge%.3f$' % dn(min(r['TB_sqrt_ell_new'] for r in F), 3), min(r['TB_sqrt_ell_new'] for r in F), 'ce_a4 T_B')
    CW = json.load(open(R('ce_counter_windows_link.json')))
    add('between $%.4f$ and $%.4f$' % (dn(min(r['product_with_T10'] for r in CW), 4), up(max(r['product_with_T10'] for r in CW), 4)), min(r['product_with_T10'] for r in CW), 'ce_counter_windows_link')
    K = json.load(open(R('cc_counter.json')))
    for r in K:
        if r['eta_target'] == 0: continue
        if (r['sig_exp'], r['eta_target']) in ((48, 0.03), (96, 0.1), (300, 0.03), (300, 0.1), (300, 0.2)):
            add('$10^{-%d}$ & $%.3f$ & $%s$ & $%.1f$ & $%.3f$' % (r['sig_exp'], r['eta_upper'], ('%.2f' % r['mu']) if r['sig_exp'] < 300 else ('%.1f' % r['mu']), r['log10_eta_e_mu2_lower'], r['gsl_hi']), r['gsl_hi'], 'cc_counter')
    mx = max(r['log10_eta_e_mu2_lower'] for r in K if r['log10_eta_e_mu2_lower'] is not None)
    add('$\\log_{10}(\\eta e^{\\mu/2})=%.1f$' % mx, mx, 'cc_counter max')
    hi = [r['gsl_hi'] for r in K if r['eta_target'] > 0]
    add('$%.2f$ to $%.2f$' % (min(hi), max(hi)), min(hi), 'cc_counter hi range (numerical values, plain rounding)')
    M = json.load(open(R('ce_split_map.json')))['search']
    add('$%s$ and $%s$' % (sci(dn(M['S_gap_windows']['min_margin'], 4)), sci(dn(M['Sprime_nogap_windows']['min_margin'], 3))), M['S_gap_windows']['min_margin'], 'ce_split_map search')
    T = json.load(open(R('cb_tau_cells.json')))
    add('at least $%.5f$' % dn(T['tau_c=20']['min_tau_factor_F'], 5), T['tau_c=20']['min_tau_factor_F'], 'cb_tau_cells [20,80]')
    add('$1-%s$' % sci(up(T['tau_c=40']['one_minus_F'], 9)), T['tau_c=40']['one_minus_F'], 'cb_tau_cells [40,80]')
    # ---- stage CH/CI/CJ: Lemma BE, the narrow focus, single A* columns
    B = json.load(open(R('ch_bedge.json'))); c = B['BE_3_constants']; assert B['all_ok']
    add('\\kappa=%s' % sci(dn(c['kappa'] / 1e-5, 1) * 1e-5), c['kappa'], 'ch_bedge kappa (claimed)')
    add('(1-\\bar\\varepsilon)\\,%d\\times10^{%d}' % (int(c['g_sqrt_ell_lower_at_1e-24_without_epsbar'] // 1e19), 19), c['g_sqrt_ell_lower_at_1e-24_without_epsbar'], 'ch_bedge g sqrt(ell) at 1e-24 (floor)')
    assert c['rho_W'] <= 0.0801; add('\\le0.0801', c['rho_W'], 'ch_bedge rho_W (<= 0.0801: 0.08 plus a sigma^2 term)')
    add('\\ge%s$' % sci(dn(c['one_minus_e_min'] / 1e-5, 2) * 1e-5, 2), c['one_minus_e_min'], 'ch_bedge 1-e (floor)')
    add('\\ge%s$, so' % sci(dn(c['y_min'] * 0.9999 / 1e-5, 2) * 1e-5, 2), c['y_min'], 'ch_bedge y_min (floor, (1-b)^2 < 1)')
    add('\\ge%.4f$' % dn(c['Phi_min_L2_ge_1'], 4), c['Phi_min_L2_ge_1'], 'ch_bedge Phi at L^2 = 1 (floor)')
    add('\\ge%s\\,\\mathrm{den}' % sci(dn(min(c['kappa_case_L2_ge_1'], c['kappa_case_L2_lt_1']) / 1e-5, 2) * 1e-5, 2), c['kappa'], 'ch_bedge kappa both cases (floor)')
    rl = B['BE_3_remainders_log10']['all_below_log10']; assert rl < -1e7
    add('below $10^{-10^7}$ times', rl, 'ch_bedge remainders')
    TT = json.load(open(R('ch_tau.json'))); row = {round(x['ell'], 2): x['tau_ceiling_D10'] for x in TT['ceiling']['rows']}
    add('$\\tau\\approx%.1f$ ($\\ell=55.26$), $%.1f$ ($\\ell=690.8$) and $%.1f$ ($\\ell=6907.8$)' % (row[55.26], row[690.78], row[6907.76]), row[690.78], 'ch_tau ceiling (derivation)')
    bt = {x['tau']: x for x in TT['beats']['rows']}[0.01]
    add('$[%.3f,%.3f]$ at $\\tau=0.01$' % (bt['local_over_period_min'], bt['local_over_period_max']), bt['local_over_period_max'], 'ch_tau beats (numerical)')
    LO = json.load(open(R('ci_local.json')))
    add('The value $%.2f$ came from' % up(LO['worst']['rho_max'], 2), LO['worst']['rho_max'], 'ci_local: the v1.0 ratio, now the source of the error (ceil)')
    NP = json.load(open(R('cj_nop2.json')))
    add('by factors $%.2f$ to $%d$' % (dn(NP['min_ratio'], 2), int(NP['max_ratio'])), NP['min_ratio'], 'cj_nop2 (numerical)')
    PR = jl('ci_link_cells.jsonl'); best = {}
    for p_ in PR:
        if p_['closed']: best[p_['ell_a']] = min(best.get(p_['ell_a'], 99), p_['tau_min'])
    add('close for $\\tau\\ge%g$ at $\\ell\\in[55.26,65.26]$, for $\\tau\\ge%g$ at $\\ell\\in[215.26,225.26]$, and for $\\tau\\ge%g$ at $\\ell\\in[455.26,465.26]$' % (best[55.26], best[215.26], best[455.26]), best[55.26], 'ci_link single columns')
    assert best[695.26] == best[855.26] == best[455.26]
    add('needs $%s$ boxes' % sci(up(max(r['cells_evaluated'] for r in L) / 1e7, 1) * 1e7), max(r['cells_evaluated'] for r in L), 'ce_link max column')
    # ---- stage CL: the rounding-covering rule (App. E)
    ms = min(r['min_slack_over_S_delta'] for r in MAIN)
    add('over all accepted boxes is $%.2f$' % dn(ms, 2), ms, 'cl_main_cells min (lpos-lneg)/(S delta) (floor)')
    WC = jl('cl_windows.jsonl'); assert len(WC) == 17 and all(w['closed'] and w['cells_evaluated'] == w['old_cells'] for w in WC)
    nch = sum(1 for r in MAIN + L if r['cells_evaluated'] != r['old_cells']); mx = max(r['cells_evaluated'] - r['old_cells'] for r in MAIN + L)
    assert min(r['cells_evaluated'] - r['old_cells'] for r in MAIN + L) >= 0
    add('in all but $%d$ columns, where at most %s more boxes' % (nch, {1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six'}.get(mx, str(mx))), nch, 'cl_* changed columns')
    add('In (i) the smallest ratio of the accepted margin to $S\\delta$ is $%.2f$' % dn(ms, 2), ms, 'cl_main_cells min slack (floor)')
    TB = json.load(open(R('cl_tau_be.json'))); assert TB['all_ok']
    add('Its margin $%s$ above $0.99977$' % sci(dn(TB['tau_cells']['margin'] / 1e-6, 1) * 1e-6), TB['tau_cells']['margin'], 'cl_tau_be margin (floor)')
    add('more than $%d\\times10^{9}$ times this bound' % int(TB['tau_cells']['margin_over_bound'] // 1e9), TB['tau_cells']['margin_over_bound'], 'cl_tau_be ratio (floor)')
    add('agrees to within $%s$' % sci(up(TB['tau_cells']['max_rel_err_double_vs_decimal'] / 1e-16, 1) * 1e-16), TB['tau_cells']['max_rel_err_double_vs_decimal'], 'cl_tau_be decimal agreement (ceil)')
    # ---- v1.1: the covering (stages CP, CQ)
    T2 = jl('cq_tails_ii.jsonl'); T1 = jl('cq_tails_i_w0.jsonl'); WG = jl('cq_wedge.jsonl'); CV = json.load(open(R('cq_cover.json')))
    assert len(T2) == len(L) and all(t['closed'] for t in T2) and len(WG) == len(L) and all(w['closed'] for w in WG) and len(T1) == len(MAIN) and all(t['closed'] for t in T1)
    s3 = lambda v: '%.2f\\times10^{%d}' % (v / 10 ** int(math.floor(math.log10(v))), int(math.floor(math.log10(v))))
    add('$%d$ of $10^5$ states for (ii)' % CV['POWER ii: v1.0 table alone']['missed'], CV['POWER ii: v1.0 table alone']['missed'], 'cq_cover v1.0 table')
    assert CV['ii repaired (body + tails + wedge + BE)']['missed'] == 0 and CV['ii repaired (body + tails + wedge + BE)']['states'] == 100000
    assert all(CV['i %s repaired (body + mu tails + b strip)' % k]['missed'] == 0 and CV['i %s repaired (body + mu tails + b strip)' % k]['states'] == 50000 for k in ('cc', 'cd'))
    npw = sum(1 for k in CV if k.startswith('POWER')); assert all(CV[k]['missed'] > 0 for k in CV if k.startswith('POWER'))
    add('each of %s deliberately damaged tables' % {11: 'eleven'}[npw], npw, 'cq_cover damaged tables')
    add('close on all $%d+%d$ columns ($%s$ boxes)' % (len(CC), len(D), s3(sum(t['cells_evaluated'] for t in T1))), sum(t['cells_evaluated'] for t in T1), 'cq_tails_i')
    add('($%s$ boxes, $%d$ CPU minutes)' % (s3(sum(t['cells_evaluated'] for t in T2)), round(sum(t['seconds'] for t in T2) / 60)), sum(t['cells_evaluated'] for t in T2), 'cq_tails_ii')
    add('($%s$ boxes, $%.1f$ CPU hours; at most $%s$ boxes in a column)' % (s3(sum(w['cells_evaluated'] for w in WG)), sum(w['seconds'] for w in WG) / 3600, s3(max(w['cells_evaluated'] for w in WG))), sum(w['cells_evaluated'] for w in WG), 'cq_wedge')
    add('is skipped ($%s$ boxes)' % s3(sum(w['skipped_boxes'] for w in WG)), sum(w['skipped_boxes'] for w in WG), 'cq_wedge skipped')
    add('$\\mu\\in[e^{-30},\\mu_{\\max}e^{60}]$ ($%s$ boxes)' % s3(sum(r['cells_evaluated'] for r in L)), sum(r['cells_evaluated'] for r in L), 'cl_link_cells total')
    QC = json.load(open(R('cq_checks.json')))
    add('bit for bit on $%d$ boxes' % QC['F1_identity']['boxes'], QC['F1_identity']['boxes'], 'cq_checks identity')
    assert QC['F1_identity']['phi_b_maxdiff'] == 0 and QC['F1_identity']['phi_bl_maxdiff'] == 0 and all(QC['wrong_broken'].values())
    add('$%s$ states in the tails of (ii), $%s$ in those of (i), $%s$ in the wedge' % tuple('%.1f\\times10^{%d}' % (v / 10 ** int(math.log10(v)), int(math.log10(v))) for v in (QC['F3_ii']['points'], QC['F3_i']['points'], QC['F3_w']['points'])), QC['F3_w']['points'], 'cq_checks F3')
    add('A $60$-digit evaluation on $%d$ points of $[20,10^8]$' % QC['K6']['grid'], QC['K6']['grid'], 'cq_checks K6')
    # ---- v1.1: the narrow focus (Sec. IX A; stages CM, CN, CO)
    CM = json.load(open(R('cm_counter.json')))
    add('$P_\\theta/\\Pi^*(\\mu_\\theta)=%.1f$ at $x=0.003$ and $%.1f$ at $x=10^{-6}$' % (CM['headline'][0]['rho'], CM['best']['rho']), CM['best']['rho'], 'cm_counter (numerical)')
    EO = json.load(open(R('co_envmin_opt.json'))); EM = json.load(open(R('co_envmin.json')))
    lab = {'55.26': '$\\ell=55.26$', '110.52': '$110.5$', '221.05': '$221$', '455.26': '$455$', '690.78': '$690.8$', '960.8': '$960.8$'}
    for k in EO:
        add('$%.4f$ vs $%.4f$ (%s)' % (EO[k]['min'], EO[k]['T0'], lab[k]), EO[k]['min'], 'co_envmin_opt (numerical)')
    add(', '.join('$%.4f$' % EM[k]['min'] for k in ('55.26', '110.52', '221.05', '455.26', '690.78')) + ' and $%.4f$' % EM['960.8']['min'], EM['55.26']['min'], 'co_envmin grid (K, L) (numerical)')
    CO = json.load(open(R('co_remainders.json'))); Dv = CO['55.26']['argmin_envmin_opt eta=1/ell']['worst_all_x']['D']
    add('$D\\le%s$ at $\\ell=55.26$' % sci(up(Dv / 1e-5, 1) * 1e-5), Dv, 'co_remainders (ceil)')
    RT = json.load(open(R('co_rtheta.json'))); neg = [x for x in RT['states'] if x['lemma_R_period'] < 0]
    add('on $%d$ states with $\\tau$ from $0.3$ to $20$' % len(RT['states']), len(RT['states']), 'co_rtheta')
    add('it keeps $%d$--$%d\\,\\%%$ of the exact value' % (round(100 * min(x['ratio_theta'] for x in neg)), round(100 * max(x['ratio_theta'] for x in neg))), min(x['ratio_theta'] for x in neg), 'co_rtheta (rounded)')
    E2 = json.load(open(R('co_envpsi2.json'))); WT = json.load(open(R('cm_wterm.json')))
    add('falls to $%.2f$ ($\\kappa\\ge10^{-12}$) and to the $W$-free value $%.2f$' % (E2['55.26']['Cinf_kmin=1e-12'], WT['55.26']['x0']), E2['55.26']['Cinf_kmin=1e-12'], 'co_envpsi2, cm_wterm (numerical)')
    add('($%.4f$ instead of $%.4f$)' % (E2['55.26']['Psi2'], EM['55.26']['min']), E2['55.26']['Psi2'], 'co_envpsi2 grid vs co_envmin')
    CN = json.load(open(R('cn_search.json'))); C1 = json.load(open(R('cn_search_round1.json')))
    add('with ratio $%.2f$' % C1['worst']['C1']['ratio'], C1['worst']['C1']['ratio'], 'cn_search round 1')
    add('with largest ratio $%.3f$' % CN['worst']['C1']['ratio'], CN['worst']['C1']['ratio'], 'cn_search C1')
    add('by a factor $%s$' % sci(CN['worst']['W1']['ratio']), CN['worst']['W1']['ratio'], 'cn_search W1')
    add('replaced by $\\kappa$, by $%.1f$' % CN['worst']['W2']['ratio'], CN['worst']['W2']['ratio'], 'cn_search W2')
    add('the constant $0.3$, by $%.1f$' % CN['worst']['W3']['ratio'], CN['worst']['W3']['ratio'], 'cn_search W3')
    add('divided by $10$, by $%.1f$' % CN['worst']['W4']['ratio'], CN['worst']['W4']['ratio'], 'cn_search W4')
    return out

def check(root, texdir):
    tex = ''.join(open(f).read() for f in sorted(glob.glob(os.path.join(texdir, 'sections', '*.tex'))))
    tex = tex.replace('\n', ' ')
    rows = numbers(root); bad = [r for r in rows if r['text'] not in tex]
    for r in rows: print(('OK   ' if r not in bad else 'MISS ') + r['text'] + '   <- ' + r['source'])
    return rows, bad

if __name__ == '__main__':
    root = sys.argv[1] if len(sys.argv) > 1 else '.'; texdir = sys.argv[2] if len(sys.argv) > 2 else 'paperH/draft'
    rows, bad = check(root, texdir); print('%d numbers, %d missing' % (len(rows), len(bad))); sys.exit(1 if bad else 0)
