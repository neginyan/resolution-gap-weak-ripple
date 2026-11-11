"""Re-run the LAST column of certificate (ii) (Lemma R* boxes, wedge, mu tails) with ell_b moved OUTWARD to the next binary64 number above 960.8,
so that the decimal end point 960.8 of Theorem H1 is inside the boxes.  lam = 1 (published certificate + G -> 0+) or lam = 2/3 (second certificate)."""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault('DA_OUT', os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'results'))
import da_scaling as DS
lam = float(sys.argv[1]); out = os.path.join(os.environ['DA_OUT'], 'top_end_lam%s.jsonl' % ('1' if lam == 1 else '2_3'))
for kind in ('tails', 'cl', 'wedge'):
    pub = [json.loads(l) for l in open(os.path.join(DS.RES, DS.PUB[kind])) if l.strip()]
    p = dict(max(pub, key=lambda r: r['ell_a'])); old = p['ell_b']; p['ell_b'] = math.nextafter(960.8, math.inf)
    assert old == 960.8 and p['ell_b'] > 960.8
    d = DS.run_column(kind, p, lam, seed=int(round(p['ell_a'] * 1000)) + {'cl': 1, 'wedge': 7, 'tails': 3}[kind])
    d['ell_b_repr'] = repr(p['ell_b']); open(out, 'a').write(json.dumps(d) + '\n')
    print(kind, lam, d['closed'], d['cells_evaluated'], d['published_cells'], d['min_margin_on_accepted'], d['mid_negative'], d['seconds'], flush=True)
