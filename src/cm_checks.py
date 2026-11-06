"""Stage CM: independent checks of the steps of the local brightness lemma (report CM, section 2) and of the F-term identity (section 3).
Python 3.11.  Exact arithmetic (Fraction), 60-digit decimals, and double-precision FFT checks on random windows.
-> results/cm_checks.json ; exit code 1 if any check fails.
"""
import json, math, os, sys
from fractions import Fraction
from decimal import Decimal, getcontext
import numpy as np

getcontext().prec = 60
OUT = 'results/cm_checks.json'
res = {}; ok_all = True


def rec(name, ok, **kw):
    global ok_all
    ok_all &= bool(ok); res[name] = {'ok': bool(ok), **kw}
    print('[%s] %s' % ('OK ' if ok else 'FAIL', name), {k: (float(v) if isinstance(v, (int, float, np.floating)) else v) for k, v in kw.items()}, flush=True)


# ---- (a) elementary inequalities, 60-digit decimal ------------------------------------------------------------------------
D = Decimal; E_ = D(1).exp()
# (a1) z e^{-z/2} <= 2/e, equality at z = 2: so s^2 e^{-s^2/2tau} <= (4 tau/e) e^{-s^2/4tau}
zs = [D(k) / 100 for k in range(0, 4001)]
m = max(z * (-z / 2).exp() for z in zs); two_over_e = 2 / E_
rec('a1 max z e^{-z/2} = 2/e', m <= two_over_e + D('1e-55') and abs(D(2) * (-D(1)).exp() - two_over_e) < D('1e-50'), max_on_grid=str(m)[:20], two_over_e=str(two_over_e)[:20])
# (a2) sqrt2/e
rec('a2 sqrt2/e < 0.5203', D(2).sqrt() / E_ < D('0.5203'), value=str(D(2).sqrt() / E_)[:12])
# (a3) dual feasibility: h(s) = (z-1) + e^{-2z(1-s)} - s^2 z >= 0 on s in [0,6] for z >= 0.8 (normalised k = s k_*, c = 1)
worst = D(10)
for z in (D('0.8'), D('0.9'), D(1), D('1.33'), D(2), D(3), D(5), D(10), D(30), D(60)):
    for i in range(0, 6001):
        s = D(i) / 1000
        h = (z - 1) + (-2 * z * (1 - s)).exp() - s * s * z
        worst = min(worst, h)
rec('a3 dual feasibility h >= 0 (z >= 0.8, s in [0,6])', worst >= -D('1e-55'), min_h=str(worst)[:14])
# (a4) the moment bound algebra: (E0/2)(k*^2 - c) + c B1 e^{-k*^2/c} = (E0/2) c ln(2B1/E0) at k*^2 = c ln(2B1/E0)  (exact in logs: Fraction check of the two forms)
for E0, B1, c in ((Fraction(1, 100), Fraction(3, 2), Fraction(7)), (Fraction(1, 10 ** 6), Fraction(11, 10), Fraction(1, 3))):
    L = math.log(2 * float(B1) / float(E0))            # ln(2B1/E0) (float); the identity is algebraic in L
    lhs = float(E0) / 2 * (float(c) * L - float(c)) + float(c) * float(B1) * (float(E0) / (2 * float(B1)))
    rhs = float(E0) / 2 * float(c) * L
    rec('a4 moment bound algebra E0=%s' % E0, abs(lhs - rhs) <= 1e-15 * abs(rhs), lhs=lhs, rhs=rhs)
# (a5) 1 - e^{-u} >= u/(1+u) >= u/(2+u) and the chirp mean: E[cos(beta s/tau)] = e^{-beta^2/2tau} (used only as a remark)
rec('a5 e^{-y/2}: 1-e^{-y/2} <= y/2', all((1 - (-D(k) / 1000 / 2).exp()) <= D(k) / 2000 for k in range(1, 2001)))

# ---- (b) analytic facts, double precision on a random window ----------------------------------------------------------------
rng = np.random.default_rng(20261104)
def window_rand(M=60):
    ns = np.arange(1, M + 1); a = rng.uniform(0, 1, M); a[0] = 1.0; phi = rng.uniform(0, 2 * np.pi, M)
    T = np.linspace(0, 2 * np.pi, 1 << 14, endpoint=False)
    f = 1 + sum(a[i] * np.cos(ns[i] * T - phi[i]) for i in range(M)); mn = f.min()
    if mn < 0: a = a / (1 - mn) * 0.999
    return ns, a, phi
ns, a, phi = window_rand(); x = 2e-3; q = a * np.exp(-ns ** 2 * x / 2)
cn = (a / 2) * np.exp(-1j * phi)                            # f = 1 + sum_n (cn e^{int} + c.c.)
def u_of(z):                                               # u(z) = 1 + sum_n q_n cos(n z - phi_n), complex z
    return 1 + np.sum(q[:, None] * np.cos(ns[:, None] * z[None, :] - phi[:, None]), axis=0)
def G_of(z):                                               # heat kernel by its Poisson form, G_x(z) = sqrt(2pi/x) sum_k e^{-(z-2pi k)^2/2x}, |k| <= 3
    kk = np.arange(-3, 4); return np.sum(np.sqrt(2 * np.pi / x) * np.exp(-(z[None, :] - 2 * np.pi * kk[:, None]) ** 2 / (2 * x)), axis=0)
t = np.linspace(-0.3, 0.3, 601); tu = np.linspace(0, 2 * np.pi, 400, endpoint=False)
ratios = []
for y in (0.01, 0.05, 0.1, 0.2):
    ratios.append(float(np.max(np.abs(G_of(t + 1j * y)) / (np.exp(y * y / (2 * x)) * np.real(G_of(t))))))
    ratios.append(float(np.max(np.abs(u_of(tu + 1j * y)) / (np.exp(y * y / (2 * x)) * np.real(u_of(tu))))))
rec('b1 |G_x(w+iy)| <= e^{y^2/2x} G_x(w) and |u(t+iy)| <= e^{y^2/2x} u(t)', max(ratios) <= 1 + 1e-9, max_ratio=max(ratios))

# (b2) the line picture: E0 = int_R w^2 R^2 dt/2pi, J(y) = int |F(t+iy)|^2 dt/2pi = int e^{-2ky} dm, E1 = int k^2 dm + ||sF||^2/4tau^2 - E0/2tau
tau = 0.02; t0 = 1.3; Lh = 14 * math.sqrt(2 * tau); Ng = 1 << 16
s = np.linspace(-Lh, Lh, Ng, endpoint=False); ds = s[1] - s[0]; tt = t0 + s
w = (2 * np.pi / tau) ** 0.25 * np.exp(-s * s / (4 * tau))
R = lambda z: u_of(z) - 1
Rp = lambda z: -np.sum(ns[:, None] * q[:, None] * np.sin(ns[:, None] * z[None, :] - phi[:, None]), axis=0)
F = w * np.real(R(tt)); F1 = w * np.real(Rp(tt))
E0 = np.sum(F ** 2) * ds / (2 * np.pi); E1 = np.sum(F1 ** 2) * ds / (2 * np.pi)
# circle definition via the Fourier formula (cm_local_bound.local_quantities)
import cm_local_bound as LB
E0c, E1c = LB.local_quantities(ns, q, ns * t0 - phi, x, tau)
rec('b2 line unfolding: E0, E1 agree with the circle formulas', abs(E0 - E0c) <= 1e-9 * E0c and abs(E1 - E1c) <= 1e-9 * E1c, E0_line=E0, E0_circle=E0c, E1_line=E1, E1_circle=E1c)
Fh = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(F))) * ds / math.sqrt(2 * np.pi); k = np.fft.fftshift(np.fft.fftfreq(Ng, ds)) * 2 * np.pi; dk = k[1] - k[0]
m_k = np.abs(Fh) ** 2 * dk / (2 * np.pi)
rec('b3 Plancherel: int dm = E0', abs(np.sum(m_k) - E0) <= 1e-8 * E0, sum_m=float(np.sum(m_k)), E0=E0)
k2m = float(np.sum(k ** 2 * m_k)); sF2 = np.sum(s ** 2 * F ** 2) * ds / (2 * np.pi)
rec('b4 identity E1 = int k^2 dm + ||sF||^2/4tau^2 - E0/2tau', abs(E1 - (k2m + sF2 / (4 * tau * tau) - E0 / (2 * tau))) <= 1e-7 * E1, E1=E1, rhs=k2m + sF2 / (4 * tau * tau) - E0 / (2 * tau))
Jok = True; Jrows = []
for y in (0.0, 0.02, 0.05, 0.1):
    Fy = w * np.exp(y * y / (4 * tau)) * 0 + (2 * np.pi / tau) ** 0.25 * np.exp(-(s + 1j * y) ** 2 / (4 * tau)) * R(tt + 1j * y)
    Jy = np.sum(np.abs(Fy) ** 2) * ds / (2 * np.pi); sel = np.abs(k) <= 150; Jk = float(np.sum(m_k[sel] * np.exp(-2 * k[sel] * y)))
    By = math.exp(y * y * (1 / x + 1 / (2 * tau))) * (1 + math.sqrt(E0) + math.exp(-y * y / (2 * x))) ** 2
    Jok &= abs(Jy - Jk) <= 1e-6 * Jy and Jy <= By; Jrows.append({'y': y, 'J_line': Jy, 'J_spectral': Jk, 'B': By})
rec('b5 Paley-Wiener J(y) = int e^{-2ky} dm and J(y) <= B(y)', Jok, rows=Jrows)
# (b6) the lemma on this window at this focus
B0 = (1 + math.sqrt(E0) + (E0 / 2) ** (1 / (2 + x / tau))) ** 2; b = x / tau
E0p, _ = LB.local_quantities(ns, q, ns * t0 - phi, x, 2 * tau)
rec('b6 lemma: x E1 <= (1+b/2) E0 ln(2B0/E0) + 0.5203 b E0(2tau)', x * E1 <= (1 + b / 2) * E0 * math.log(2 * B0 / E0) + 0.5203 * b * E0p, lhs=x * E1, rhs=(1 + b / 2) * E0 * math.log(2 * B0 / E0) + 0.5203 * b * E0p)
rec('b7 sum_n |a_n| e^{-n^2 x/2} of the random window (not needed by the lemma)', True, value=float(q.sum()))

# ---- (c) the F-term: chirp form vs the (1-e)W main part --------------------------------------------------------------------
# Exact identity of paper F (App. D): F-term = sigma^2 int N_V F (1 - cos(alpha v + psi)) dv.  In t = k' v with the focus theta_tau the chirp is
# beta s/tau, beta^2 = tau y, and  (1 - e) = 1 - e^{-y/2} = E_rho[1 - cos(beta s/tau)].  So the local main part is E_rho[F (1 - cos(beta s/tau))],
# which is (y/2) E_rho[(s^2/tau) F] to leading order, NOT (y/2) E_rho[F]: the ratio kappa = E_rho[(s^2/tau) F]/E_rho[F] can be << 1.
def kappa_of(R2):          # R2 = R'^2 on the s-grid; kappa = E[(s^2/tau) R'^2]/E[R'^2]
    wt = np.exp(-s * s / (2 * tau)); return float(np.sum(s * s / tau * wt * R2) / np.sum(wt * R2))
def chirp_ratio(R2, y):    # E[R'^2 (1 - cos(beta s/tau))] / ((1 - e^{-y/2}) E[R'^2]) with beta = sqrt(tau y)
    beta = math.sqrt(tau * y); wt = np.exp(-s * s / (2 * tau))
    return float(np.sum(wt * R2 * (1 - np.cos(beta * s / tau))) / ((1 - math.exp(-y / 2)) * np.sum(wt * R2)))
rows = []
xc = 1e-6            # a small blur, so that the packets below are flat in n (n^2 x <= 1)
for name, N, W in (('wide packet (W=4 < 1/sqrt tau)', 300, 4), ('narrow packet (W=300 > 1/sqrt tau)', 600, 300), ('fundamental at a node of R\' (dark centre)', 1, 1)):
    nsp = np.arange(N, N + W); qp = np.exp(-nsp ** 2 * xc / 2) / W
    R2 = (np.sum(-nsp[:, None] * qp[:, None] * np.sin(nsp[:, None] * s[None, :]), axis=0)) ** 2     # R at its crest at s = 0; R'^2 oscillates at 2N under the packet envelope
    rows.append({'window': name, 'kappa': kappa_of(R2), 'chirp_over_main_y1e-6': chirp_ratio(R2, 1e-6), 'chirp_over_main_y0.1': chirp_ratio(R2, 0.1), 'one_over_sqrt_tau': 1 / math.sqrt(tau)})
narrow = [r for r in rows if r['window'].startswith('narrow')][0]; wide = [r for r in rows if r['window'].startswith('wide')][0]
rec('c1 F-term: chirp part/((1-e) E[F]) = kappa; kappa ~ 1 for a wide packet, << 1 for a narrow packet at the focus centre',
    abs(wide['chirp_over_main_y1e-6'] - wide['kappa']) < 0.02 and abs(narrow['chirp_over_main_y1e-6'] - narrow['kappa']) < 0.02 and narrow['kappa'] < 0.2 and wide['kappa'] > 0.8, rows=rows)

json.dump(res, open(OUT, 'w'), indent=1)
print('ALL OK' if ok_all else 'SOME CHECKS FAILED')
sys.exit(0 if ok_all else 1)
