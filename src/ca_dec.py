"""Tiny arbitrary-exponent arithmetic on Python's decimal module (mpmath is not available here): sigma down to 1e-500 and below without underflow.
Python 3.11."""
from decimal import Decimal, getcontext, localcontext, Underflow, Subnormal, Inexact, Rounded, Clamped
C = getcontext(); C.prec = 40; C.Emax = 999999999999999999; C.Emin = -999999999999999999
for t in (Underflow, Subnormal, Inexact, Rounded, Clamped): C.traps[t] = False


def D(x):
    if isinstance(x, Decimal): return x
    if isinstance(x, int): return Decimal(x)
    if isinstance(x, str): return Decimal(x)
    return Decimal(repr(float(x)))                    # floats and numpy floats


def exp(x): return D(x).exp()
def log(x): return D(x).ln()
def log10(x): return D(x).log10()
def sqrt(x): return D(x).sqrt()


def expm1(x):
    x = D(x)
    if abs(x) > D('1e-5'): return x.exp() - 1
    return x + x * x / 2 + x ** 3 / 6 + x ** 4 / 24


PI = Decimal('3.141592653589793238462643383279502884197')
E = Decimal(1).exp()
INF = Decimal('Infinity')


def log1p(x):
    x = D(x)
    if abs(x) > D('1e-5'): return (1 + x).ln()
    return x - x * x / 2 + x ** 3 / 3 - x ** 4 / 4
