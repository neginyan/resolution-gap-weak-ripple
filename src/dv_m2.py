"""Paper H v1.3: a proved upper bound of the constant M2 of Lemma R_theta (remainder N16, Sec. IX A).  Python 3.11 grammar.
  M2(y) := sup_{u >= 0} u e^{-u/2} (1 + u cosh sqrt(y u)),   y = (1 - b)^2 sigma^2/(2b) = beta^2/tau.
v1.2 took M2 as the maximum on a grid of u in [0, 400] times (1 + 1e-6) (src/co_remainders.py, src/co_rtheta.py): a numerical value, not a proved bound.
v1.3 replaces it by M2_bound(y):
  (a) y <= 1e-20 (every state of src/co_remainders.py: y <= 10^{-21.16}):  M2(y) <= (1 + 1e-18)(4u* + 2) e^{-u*/2} <= 2.737605,  u* = (3 + sqrt17)/2.
      PROOF (by hand; App. D of the paper v1.3).  At y = 0, d/du [u(1+u)e^{-u/2}] = e^{-u/2}(2 + 3u - u^2)/2 has one zero u* on u > 0 (positive before, negative
      after), and u*(1 + u*) = 4u* + 2 since u*^2 = 3u* + 2.  For u <= 100: c = sqrt(yu) <= 1e-9, cosh c <= 1 + (c^2/2) cosh c <= 1 + 1e-18, and the function is
      linear in cosh c with non-negative coefficients, so it is at most (1 + 1e-18) times its value at y = 0.  For u >= 100: c <= 1e-10 sqrt(u) <= 0.01 u,
      cosh c <= e^{0.01u}, so it is at most u(1 + u) e^{-0.49u}, which decreases for u >= 100 (log-derivative <= 1/u + 1/(1+u) - 0.49 <= 2/u - 0.49 < 0) and is below 6e-18 at u = 100.
      (4u* + 2) e^{-u*/2} = 2.737604646...
  (b) y > 1e-20 (the states of src/co_rtheta.py, y in [8.3e-4, 0.31]): interval arithmetic (src/dv_interval.py, 40 digits, outward widening, with the
      decimal.Inexact trap of the default context on while it runs).  u in [0, U_T], U_T = max(100, 64 y): pieces [u0, u1] of width 0.001 (u < 20) and 0.01
      (u >= 20), on which the function is at most u1 e^{-u0/2} (1 + u1 cosh sqrt(y u1));  u >= U_T: sqrt(y u) <= u/8, so the function is at most
      u(1 + u) e^{-3u/8}, decreasing for u >= 6, bounded by its value at U_T.
  POWER (run()): the claim 2.7376 (below the true value at y = 0) fails; the bound of (b) at y = 1e-4 is below the bound 2.74484 of the follow-up work and
  above the grid value; (b) with the pieces' ends exchanged (a wrong bound) falls below the grid value somewhere.
-> results/dv_m2.json"""
import json, math, sys, os, functools, decimal, time
from decimal import Decimal as Dc
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np

M2_HAND = 2.737605          # (a): (1 + 1e-18) (4u* + 2) e^{-u*/2} rounded up
Y_HAND = 1e-20


def trapped(fn):
    @functools.wraps(fn)
    def w(*a, **k):
        ctx = decimal.getcontext(); old = ctx.traps[decimal.Inexact]; ctx.traps[decimal.Inexact] = True
        try: return fn(*a, **k)
        finally: ctx.traps[decimal.Inexact] = old
    return w


@trapped
def _pieces(y_str, wrong=False):
    from dv_interval import Iv, exp, sqrt
    def cosh(x): return (exp(x) + exp(-x)) / 2
    Y = Iv(Dc(y_str)); UT = Dc(max(100, math.ceil(64 * float(y_str)) + 1))          # an integer; the tail condition is checked below
    best = Dc(0); u0 = Dc(0)
    while u0 < UT:
        h = Dc('0.001') if u0 < 20 else Dc('0.01'); u1 = u0 + h
        ue = Iv(u0) if wrong else Iv(u1); ud = Iv(u1) if wrong else Iv(u0)      # wrong: increasing part at the left end, decreasing part at the right end
        v = (ue * exp(-ud / 2) * (1 + ue * cosh(sqrt(ue * Y)))).hi
        if v > best: best = v
        u0 = u1
    U = Iv(UT); tail = (U * (1 + U) * exp(-Iv(3) * U / 8)).hi
    cond = (sqrt(Y * U) - U / 8).hi <= 0
    return max(best, tail), cond


def M2_bound(y):
    """a proved upper bound of M2(y) (float, rounded up)."""
    y = float(y)
    if y <= Y_HAND: return M2_HAND
    v, cond = _pieces(repr(y))
    assert cond
    return math.nextafter(float(v), math.inf)


def grid(y, umax=400.0, n=400001):
    u = np.linspace(0, umax, n); return float(np.max(u * np.exp(-u / 2) * (1 + u * np.cosh(math.sqrt(y) * np.sqrt(u)))))


def run(save=True):
    t0 = time.time(); out = {}
    from decimal import getcontext
    with decimal.localcontext() as c:
        c.prec = 60; us = (3 + Dc(17).sqrt()) / 2; val = (4 * us + 2) * (-us / 2).exp()
    out['(a) (4u*+2)e^{-u*/2}, 60 digits'] = str(val)[:22]
    out['(a) (1 + 1e-18) times it <= 2.737605'] = bool(val * (1 + Dc('1e-18')) <= Dc('2.737605'))
    out['(a) grid at y = 0 and y = 1e-20 (must be <= 2.737605)'] = [grid(0.0), grid(1e-20)]
    rows = {}
    for y in (1e-20 * 1.0000001, 1e-6, 1e-4, 8.3e-4, 0.01, 0.1, 0.31, 1.0):
        b = M2_bound(y); g = grid(y); rows['%.3g' % y] = {'bound': b, 'grid (v1.2 without the factor 1 + 1e-6)': g, 'bound/grid - 1': b / g - 1}
    out['(b) interval bound vs grid'] = rows
    okb = all(v['bound'] >= v['grid (v1.2 without the factor 1 + 1e-6)'] and v['bound/grid - 1'] < 3e-3 for v in rows.values())
    pw = {'claim 2.7376 at y = 0 (must fail: grid exceeds it)': grid(0.0) > 2.7376,
          'bound at y = 1e-4 <= 2.74484 (the bound of the follow-up work)': rows['0.0001']['bound'] <= 2.74484,
          'wrong pieces (ends exchanged) fall below the grid somewhere': any(float(_pieces(repr(y), wrong=True)[0]) < grid(y) for y in (1e-4, 0.1))}
    out['power'] = pw
    out['ok'] = bool(out['(a) (1 + 1e-18) times it <= 2.737605'] and max(out['(a) grid at y = 0 and y = 1e-20 (must be <= 2.737605)']) <= M2_HAND and okb and all(pw.values()))
    out['seconds'] = time.time() - t0
    print(json.dumps(out, indent=1))
    if save: json.dump(out, open(os.path.join(os.path.dirname(HERE), 'results', 'dv_m2.json'), 'w'), indent=1)
    return out


if __name__ == '__main__':
    sys.exit(0 if run()['ok'] else 1)
