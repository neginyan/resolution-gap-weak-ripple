"""Paper H v1.3: a copy of the interval arithmetic of the follow-up work (dl_interval.py of the paper I code, same author, MIT), used only by src/dv_m2.py.
AUDIT GUIDE (stage DN; for an independent reader).  STAGE DO: two lines fixed -- Iv.__neg__ now uses copy_negate (exact), and the widening in lq is done
  inside the 40-digit context; before DO both rounded in the default context (28 digits, to nearest), so subtractions were not outward (DN found it with
  dn_ivtest.py).  No other decimal operation runs outside the context CTX: checked by running every interval proof with Inexact trapped in the default
  context (do_check.py).  The pre-DO version is verification/pre_DO/dl_interval.py.
  WHAT: a small interval arithmetic in Python's decimal (40 digits, outward widening 10^-38 per operation) used by the computer-assisted proof of the
  constants of Theorem DJ-T.  WHERE IN THE DOCUMENT: DO_証明まとめ.md section 3.9 (the proof) and section 3.9.1 (why the enclosure is correct: decimal's
  rounding guarantees, the widening, infinite ends, ends at 0, sign-changing products and quotients).  TESTS: dn_ivtest.py (correct enclosures on random
  cases checked against 80-digit values; three deliberately wrong versions caught).
Stage DL: a small interval arithmetic in decimal, with outward rounding, for the computer-assisted proof of the constants of Theorem DJ-T.  Python 3.11 grammar.
ENCLOSURE.  An interval is [lo, hi] with lo <= hi, Decimal endpoints, +-Infinity allowed.  Every operation (+, -, x, /, sqrt, exp, ln, power) is evaluated at
the endpoints (the elementary functions are monotone) with decimal arithmetic of PREC digits, rounded to nearest (an error of at most half a unit in the last
place, i.e. a relative error <= 0.5 x 10^{1-PREC}, also for exp and ln, which decimal rounds correctly), and the result is then WIDENED outward by the
relative amount WID = 10^{2-PREC} > that error.  So the computed interval contains the exact image (the fundamental property of interval arithmetic; the
products take the min and max over the four endpoint products, with 0 x inf := 0, which is correct for intervals whose endpoint is the exact value 0).
SPECIAL FUNCTIONS (monotone in a way proved by one derivative each; DL_証明まとめ.md, section B):
  pe(k, x, c):  x^k e^{-c x}  (x > 0, c > 0 a number, k >= 0): increasing on (0, k/c], decreasing on [k/c, inf).
  pec(k, c, x): x^k e^{-c x}  with c an interval of positive numbers: decreasing in c; in x as pe.
  H(x, L, a, b, c): a x L - b ln L - c ln x  (x, L > 0; a, b, c > 0): increasing in L where a x L >= b, increasing in x where a x L >= c.
  lol(L): ln(L)/L, decreasing for L >= e.
  lq(U): ln(1 + U)/U (= 1 at U = 0), decreasing on (-1, inf)."""
from decimal import Decimal as Dc, Context, localcontext, ROUND_HALF_EVEN
PREC = 40
CTX = Context(prec=PREC, Emax=10 ** 17, Emin=-10 ** 17, rounding=ROUND_HALF_EVEN)
for _s in list(CTX.traps): CTX.traps[_s] = False
CTX.traps[__import__('decimal').InvalidOperation] = True
WID = Dc(10) ** (2 - PREC)
INF = Dc('Infinity')


def _dn(x):
    if x.is_infinite() or x == 0: return x
    return x - abs(x) * WID


def _up(x):
    if x.is_infinite() or x == 0: return x
    return x + abs(x) * WID


def D(v):
    if isinstance(v, Dc): return v
    if isinstance(v, int): return Dc(v)
    if isinstance(v, str): return Dc(v)
    return Dc(repr(float(v)))      # a constant written in the source as a decimal number (0.184, 3.1, ...): its decimal value, exactly


class Iv:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        with localcontext(CTX):
            lo = D(lo); hi = lo if hi is None else D(hi)
            if lo > hi: raise ValueError('empty interval %s %s' % (lo, hi))
            self.lo, self.hi = lo, hi

    def __repr__(self): return '[%s, %s]' % (format(self.lo, '.6e') if self.lo.is_finite() else self.lo, format(self.hi, '.6e') if self.hi.is_finite() else self.hi)
    def f(self): return (float(self.lo), float(self.hi))

    @staticmethod
    def c(v): return v if isinstance(v, Iv) else Iv(v)

    def __add__(self, o):
        o = Iv.c(o)
        with localcontext(CTX): return Iv(_dn(self.lo + o.lo), _up(self.hi + o.hi))
    __radd__ = __add__

    def __neg__(self):
        # DO: copy_negate is exact ("unaffected by the context ... no rounding is performed"); the unary minus of decimal rounds in the CURRENT context,
        # which made every subtraction round its subtrahend to the default 28 digits (nearest) before stage DO
        return Iv(self.hi.copy_negate(), self.lo.copy_negate())
    def __sub__(self, o): return self + (-Iv.c(o))
    def __rsub__(self, o): return Iv.c(o) + (-self)

    def __mul__(self, o):
        o = Iv.c(o)
        with localcontext(CTX):
            def m(a, b):
                if (a == 0 and b.is_infinite()) or (b == 0 and a.is_infinite()): return Dc(0)
                return a * b
            ps = [m(a, b) for a in (self.lo, self.hi) for b in (o.lo, o.hi)]
            return Iv(_dn(min(ps)), _up(max(ps)))
    __rmul__ = __mul__

    def inv(self):
        with localcontext(CTX):
            if self.lo <= 0 <= self.hi: raise ZeroDivisionError('interval contains 0: %r' % self)
            a = Dc(0) if self.hi.is_infinite() else 1 / self.hi
            b = Dc(0) if self.lo.is_infinite() else 1 / self.lo
            if self.hi < 0: a, b = (Dc(0) if self.lo.is_infinite() else 1 / self.lo), 1 / self.hi
            return Iv(_dn(min(a, b)), _up(max(a, b)))

    def __truediv__(self, o): return self * Iv.c(o).inv()
    def __rtruediv__(self, o): return Iv.c(o) * self.inv()

    def __pow__(self, p):
        if isinstance(p, int) and p >= 0:
            r = Iv(1)
            for _ in range(p): r = r * self
            if p % 2 == 0 and self.lo < 0 < self.hi: r = Iv(0, r.hi)
            return r
        return exp(log(self) * p)

    def ge0(self): return self.lo >= 0
    def le(self, v): return self.hi <= D(v)
    def gee(self, v): return self.lo >= D(v)


def _mono(fn, x, inc=True):
    with localcontext(CTX):
        a, b = fn(x.lo), fn(x.hi)
        if not inc: a, b = b, a
        return Iv(_dn(a), _up(b))


def _exp(v):
    if v == INF: return INF
    if v == -INF: return Dc(0)
    return v.exp()


def _ln(v):
    if v == 0: return -INF
    if v == INF: return INF
    return v.ln()


def exp(x): return _mono(_exp, Iv.c(x))


def log(x):
    x = Iv.c(x)
    if x.lo < 0: raise ValueError('log of %r' % x)
    return _mono(_ln, x)


def sqrt(x):
    x = Iv.c(x)
    if x.lo < 0: raise ValueError('sqrt of %r' % x)
    return _mono(lambda v: v if v.is_infinite() else v.sqrt(), x)


def log1p(x):
    x = Iv.c(x)
    return log(1 + x) if x.lo > Dc('1e-10') or x.hi < -Dc('1e-10') else _mono(lambda v: _l1p(v), x)


def _l1p(v):
    if abs(v) < Dc('1e-10'):
        return v - v * v / 2 + v ** 3 / 3 - v ** 4 / 4            # the remainder < v^5 < 1e-50 (relative to v: < 1e-40, inside the widening)
    with localcontext(CTX) as c2:
        c2.prec = PREC + 30; r = (1 + v).ln()                    # 1 + v exact to 70 digits: ln(1 + v) to relative 1e-58 for |v| >= 1e-10
    return +r


def imin(a, b):
    a, b = Iv.c(a), Iv.c(b); return Iv(min(a.lo, b.lo), min(a.hi, b.hi))


def imax(a, b):
    a, b = Iv.c(a), Iv.c(b); return Iv(max(a.lo, b.lo), max(a.hi, b.hi))


def _pe_iv(k, x, c):
    """the interval value of x^k e^{-c x} at the point x (computed in interval arithmetic)."""
    if x.is_infinite(): return Iv(0)
    if x == 0: return Iv(1) if k == 0 else Iv(0)
    return exp(log(Iv(x)) * k - Iv(x) * c)


def pe(k, x, c=1):
    """x^k e^{-c x}, k >= 0, c > 0 (numbers): increasing on (0, k/c], decreasing on [k/c, inf)."""
    x = Iv.c(x); k = D(k); c = D(c)
    with localcontext(CTX):
        m = k / c
        if x.lo >= m: return Iv(_pe_iv(k, x.hi, c).lo, _pe_iv(k, x.lo, c).hi)
        if x.hi <= m: return Iv(_pe_iv(k, x.lo, c).lo, _pe_iv(k, x.hi, c).hi)
        mm = Iv(m) * Iv(1)                                                                  # the maximum at k/c (k/c rounded: the max of a smooth function is flat there)
        top = _pe_iv(k, m, c).hi * (1 + WID * 1000)
        return Iv(min(_pe_iv(k, x.lo, c).lo, _pe_iv(k, x.hi, c).lo), top)


def pec(k, c, x):
    """x^k e^{-c x} with c an interval of positive numbers: decreasing in c, so in [pe(k, x, c.hi), pe(k, x, c.lo)]."""
    c = Iv.c(c); a = pe(k, x, c.hi); b = pe(k, x, c.lo)
    return Iv(a.lo, b.hi)


def H(x, L, a, b, c):
    """a x L - b ln L - c ln x: increasing in L where a x L >= b and in x where a x L >= c; else the plain interval form."""
    x, L = Iv.c(x), Iv.c(L); a, b, c = D(a), D(b), D(c)
    with localcontext(CTX):
        ok = x.lo > 0 and L.lo > 0 and a * x.lo * L.lo >= max(b, c) * (1 + WID * 10)
    if ok:
        def h(xx, LL):
            if LL.is_infinite() or xx.is_infinite(): return Iv(INF)
            return Iv(a) * Iv(xx) * Iv(LL) - b * log(Iv(LL)) - c * log(Iv(xx))
        return Iv(h(x.lo, L.lo).lo, h(x.hi, L.hi).hi)
    return a * x * L - b * log(L) - c * log(x)


def lol(L):
    """ln(L)/L for L >= e (decreasing)."""
    L = Iv.c(L)
    assert L.lo >= Dc('2.7183')
    f = lambda v: Iv(0) if v.is_infinite() else log(Iv(v)) / Iv(v)
    return Iv(f(L.hi).lo, f(L.lo).hi)


def lq(U):
    """ln(1 + U)/U on (-1, inf), decreasing, = 1 at 0."""
    U = Iv.c(U)
    def f(v):
        if v.is_infinite(): return Iv(0)
        with localcontext(CTX) as c2:
            if abs(v) < Dc('1e-6'):
                c2.prec = PREC + 10
                s_ = sum(((-v) ** k_ / (k_ + 1) for k_ in range(9)), Dc(0))                 # remainder < |v|^9 < 1e-54
                return Iv(s_ - abs(s_) * WID - Dc('1e-50'), s_ + abs(s_) * WID + Dc('1e-50'))
            c2.prec = PREC + 30; r = (1 + v).ln() / v                                         # 1 + v exact to 70 digits: relative error < 1e-60
        with localcontext(CTX):                                                               # DO: the widening inside the 40-digit context (was outside: 28 digits)
            return Iv(r - abs(r) * WID, r + abs(r) * WID)
    return Iv(f(U.hi).lo, f(U.lo).hi)
