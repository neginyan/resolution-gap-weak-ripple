"""Stage CI 2.1: the three local corrections of Lemma R_theta, checked numerically against their proofs.  Python 3.11.
(L1) focus term of the w-information:  F_w = sum_n q_n e_n cos(phi_n) int Gamma E (2/L^2 - w^2/L^4);  claim |F_w| <= 2 Qbar sqrt(r/(1+r))/L^2.
     Proof: int Gamma E = k_E = sqrt(s^2/(s^2+L^2)) = sqrt(r/(1+r)); int Gamma E w^2 = k_E s_E^2, s_E^2 = L^2 s^2/(L^2+s^2) < L^2, so the integral is in
     [k_E/L^2, 2 k_E/L^2]; |sum q_n e_n cos phi_n| <= sum q_n <= Qbar.
(L2) curvature term of W:  |x sum n^2 q_n e_n cos phi_n| <= (2/e) Qbar b,  because x n^2 e_n = b u e^{-u/2} (u = n^2 tau) and max u e^{-u/2} = 2/e.
(L3) normalisation: Z = 1 + sum q_n e_n cos phi_n int Gamma E  in  [1 - Qbar k_E, 1 + Qbar k_E].
Random windows (q0_n, phases), random tau in [1e-6, 20], random (s^2, L^2) with r = s^2/L^2 <= 1/4.  -> results/ci_rtheta.json"""
import json, math
import numpy as np

QB = 0.25


def main(ntrial=3000, seed=20261103):
    rng = np.random.default_rng(seed)
    w1 = w2 = w3 = 0.0
    for _ in range(ntrial):
        nh = int(rng.integers(1, 400)); n = np.arange(1, nh + 1)
        q = rng.uniform(0, 1, nh) * (rng.uniform(0, 1, nh) < 0.3 + 0.7 * rng.uniform()); q = q / max(q.sum(), 1e-300) * QB * rng.uniform(0.2, 1)
        ph = rng.uniform(0, 2 * math.pi, nh)
        tau = 10 ** rng.uniform(-6, math.log10(20)); x = 10 ** rng.uniform(-9, math.log10(0.003)); b = x / tau
        L2 = 10 ** rng.uniform(-3, 3); s2 = L2 * 10 ** rng.uniform(-6, math.log10(0.25)); r = s2 / L2
        e = np.exp(-n * n * tau / 2); S = float(np.sum(q * e * np.cos(ph)))
        w = np.linspace(-12 * math.sqrt(L2), 12 * math.sqrt(L2), 20001); dw = w[1] - w[0]
        Gam = np.exp(-w ** 2 / (2 * L2)) / math.sqrt(2 * math.pi * L2); E = np.exp(-w ** 2 / (2 * s2))
        I = float(np.sum(Gam * E * (2 / L2 - w ** 2 / L2 ** 2)) * dw); kE = math.sqrt(r / (1 + r))
        w1 = max(w1, abs(S * I) / (2 * QB * kE / L2))
        w2 = max(w2, abs(x * float(np.sum(n * n * q * e * np.cos(ph)))) / (2 / math.e * QB * b))
        kE_num = float(np.sum(Gam * E) * dw)
        w3 = max(w3, abs(S * kE_num) / (QB * kE))
    out = {'trials': ntrial, 'L1_max_ratio': w1, 'L2_max_ratio': w2, 'L3_max_ratio': w3, 'ok': bool(max(w1, w2, w3) <= 1 + 1e-9)}
    json.dump(out, open('results/ci_rtheta.json', 'w'), indent=1)
    print(out)


if __name__ == '__main__':
    main()
