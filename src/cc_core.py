"""Stage CC core: the x-free reduced lower bound and its minimisation.
Main reduced function (Lemma R's leading terms, mid, eta = 0): with A = sigma^2/s^2, r = s^2/L^2 in (0, 1/4], G = P g_s, g_s = sqrt(r/(2 + r)),
    N = 2 A r + G A/(mu(2 + r)) + 3 sigma^2/8 + b sigma^2/(8 r A) + G sigma^2/(16 b) + G sigma^2/(16 A (2 + r)),
    D = b + G/2 + A r + G A/(2 mu (2 + r)),
and g >= gamma iff Phi = N - gamma sigma D >= 0.  Minimising over A > 0 exactly (alpha A + beta/A >= 2 sqrt(alpha beta)) gives Phi/sigma >= phi,
    phi(mu, r, b, G) = 2 sqrt(alpha (b/(8r) + G/(16(2 + r)))) + 3 sigma/8 + G sigma/(16 b) - gamma b - gamma G/2,
    alpha = r (2 - gamma sigma) + G (1 - gamma sigma/2)/(mu (2 + r)).
phi is concave in G (geometric mean of two affine functions plus affine terms), so over G in [0, G_max] its minimum is at an end point;
G_max = Pi*(mu) g_s(r) with the x-free bound of Lemma Z3' (x <= 0.003):  Pi*(mu) = Qbar^2 mu (1 + a t0)/(1 + t0)^2, a = c_s(0.003) e^{-mu/2},
t0 = (e/2)(mu - 0.003) e^{(mu - 0.003)/2}  (and Pi* <= mu Qbar^2).  At G = 0, phi >= sqrt(b (1 - gamma sigma/2)) - gamma b >= 0 (b <= 1, gamma < 1).
Python 3.11."""
import math
QB = 0.25
CS = math.exp(0.0015) * (1 + math.sqrt(math.e) + 2 / math.e)


def log_Pistar(mu):
    if mu <= 0.003: return math.log(max(mu, 1e-300)) + 2 * math.log(QB)
    lt0 = math.log(math.e / 2) + math.log(mu - 0.003) + (mu - 0.003) / 2
    la = math.log(CS) - mu / 2
    lae = lambda *v: max(v) + math.log(sum(math.exp(t - max(v)) for t in v))
    val = 2 * math.log(QB) + math.log(mu) + lae(0.0, la + lt0) - 2 * lae(0.0, lt0)
    return min(val, math.log(mu) + 2 * math.log(QB))


def phi_over_G(lmu, lr, lb, lsig, gamma, lG=None):
    """phi/G at G = G_max (in logarithms; G may be ~ e^{-mu}).  Returns phi/G (dimensionless)."""
    mu = math.exp(lmu); r = math.exp(lr); b = math.exp(lb); sig = math.exp(lsig)
    lgs = 0.5 * (lr - math.log(2 + r))
    if lG is None: lG = log_Pistar(mu) + lgs
    # alpha/G = r(2 - gamma sig)/G + (1 - gamma sig/2)/(mu(2+r));  (b/(8r) + G/(16(2+r)))/G
    a1 = math.exp(lr - lG) * (2 - gamma * sig) + (1 - gamma * sig / 2) / (mu * (2 + r))
    c1 = math.exp(lb - lr - lG) / 8 + 1 / (16 * (2 + r))
    return 2 * math.sqrt(a1 * c1) + 3 * math.exp(lsig - lG) / 8 + math.exp(lsig - lb) / 16 - gamma * math.exp(lb - lG) - gamma / 2
