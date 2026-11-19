# resolution-gap-weak-ripple

**Beyond few harmonics: the lower side of the logarithmic law of the resolution gap for windows with a weak ripple**

Code, stored results, computer-assisted certificates and a verification script for the paper by Taishi Namba (2026).

**Paper:** T. Namba, *Beyond few harmonics: the lower side of the logarithmic law of the resolution gap for windows with a weak ripple*, preprint, Zenodo (2026), [doi:10.5281/zenodo.23932698](https://doi.org/10.5281/zenodo.23932698) (all versions; v1.0 is [10.5281/zenodo.23932699](https://doi.org/10.5281/zenodo.23932699); v1.1 is [10.5281/zenodo.23993985](https://doi.org/10.5281/zenodo.23993985); v1.2 is [10.5281/zenodo.24114426](https://doi.org/10.5281/zenodo.24114426); the DOI of v1.3 is added on publication).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23932577.svg)](https://doi.org/10.5281/zenodo.23932577)
[![verify](https://github.com/neginyan/resolution-gap-weak-ripple/actions/workflows/verify.yml/badge.svg)](https://github.com/neginyan/resolution-gap-weak-ripple/actions/workflows/verify.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## The result in brief

This is the follow-up to [paper F](https://doi.org/10.5281/zenodo.23856796) (code: [resolution-gap-two-sided](https://github.com/neginyan/resolution-gap-two-sided)). Paper F proved the logarithmic law of the resolution gap of the arrow of time, with the sharp constant, for periodic windows with at most twenty harmonics or with damping exponent K ≥ 0.003. Here the **lower side** is extended to windows with **any number of harmonics** whose ripple is weak in the sup norm or poor in low-frequency energy.

- **Theorem H1.** For every state whose focus spans at least about 0.7 ripple periods (τ = x/b ≥ 20), and for ℓ = ln(1/σ) ∈ [55.26, 960.8] (σ from about 10⁻²⁴ down to about 10⁻⁴¹⁷), g√ℓ ≥ (1 − ε̄)(λ/Υ)(1 + 12.1/ℓ + 2.5×10⁻⁴). The coefficient one follows, up to 1 − ε̄, for every window with η ≤ 1/ℓ, for a family of intermediate windows, and for seventeen thin-block windows certified one by one. States with b < 10⁴σ² obey a far stronger bound (Lemma BE). Since v1.1 the boxes of the computer-assisted proofs cover every state of the theorem; this is shown line by line and checked on 10⁵ random states. v1.2 corrects two steps of the proof and closes both by re-running the certificates; v1.3 makes a constant of the local route (§IX A) a proved bound and corrects three smaller points of the open problems and the code (below). Theorem H1 is unchanged since v1.0 in its statement and range.
- **Three tools.** A remainder of the curvature step made relative to the ripple (Lemma A4*); a split sandwich in which a low-frequency block enters through its energy, not its sup (Lemma S); a mean-square bound of the curvature cross term that closes the edge of very wide focus (Lemma BE).
- **A conjecture of our own is false.** The working hypothesis ηe^{μ/2} ≤ C fails on explicit windows (up to 10³² at σ = 10⁻³⁰⁰); the coefficient one survives on all of them.
- **The road to an unconditional theorem.** The paper lists what is still open, from the nearest to the farthest: states with a narrow focus (τ < 20), ℓ > 960.8 and the limit σ → 0, windows outside the covered families, the class constant ε̄, and the upper side. The narrow focus is now described with proved lemmas and a numerical envelope, but it is not closed; see "What v1.1 fixes" below.

Every statement in the paper carries its status: proof, computer-assisted proof, derivation, numerical, or conjecture (Table I of the paper).

## What v1.3 fixes

v1.3 closes the known issues listed below (2026-11-11, 2026-11-15, two of 2026-11-18) and one more point found while preparing it. **The statement and the range of Theorem H1 are unchanged, and no certificate of Theorem H1 is re-run.** One fix (item 6) concerns a lemma of its proof, Lemma BE: its hypothesis is now stated for the whole range of the theorem; the proof is the same. `verification/run_v13.py` checks this by undoing the four code changes and running `verification/run_v12.py`, unchanged, on the result. Every number printed in §IX A of the paper is unchanged in the printed digits.

1. **The constant M₂ of Lemma R_θ (remainder N16) is a proved bound** (known issue of 2026-11-18). `M2_const` in `src/co_remainders.py` and `src/co_rtheta.py` now calls `src/dv_m2.py`:
   - for y ≤ 10⁻²⁰, M₂(y) ≤ (1 + 10⁻¹⁸)(4u\* + 2)e^{−u\*/2} ≤ 2.737605 with u\* = (3 + √17)/2, proved by hand (App. D of the paper v1.3, new subsection "The constant M₂ of Lemma R_θ");
   - for larger y, a bound by interval arithmetic (`src/dv_interval.py`, 40 digits with outward rounding, run with the `decimal.Inexact` trap on), on pieces of u and a hand tail.
   The range of y = (1 − b)²σ²/(2b) used in §IX A is now stated: every state of the remainder budget (`co_remainders.py`) has y ≤ 10⁻²¹·¹⁶, so the hand bound applies; the 46 quadrature states (`co_rtheta.py`) have y from 8.3 × 10⁻⁴ to 0.30, where the interval bound applies. The interval bound exceeds the old grid value by about 5 × 10⁻⁴ relative. Re-running `co_remainders.py` and `co_rtheta.py` changes no printed number (D ≤ 1.2 × 10⁻⁵ at ℓ = 55.26; 46 states; 27–94 %).
2. **The coefficient √2/(16e) in `src/cn_lemmas.py`** (known issue of 2026-11-15): `0.0325` is replaced by `math.sqrt(2) / (16 * math.e)` = 0.0325163…, the value of Lemma CM-1. The five §IX scripts (`cn_lemmas.py`, `co_envpsi2.py`, `co_envmin.py`, `co_envmin.py opt`, `co_remainders.py`) and `co_rtheta.py` have been re-run; every printed number is unchanged (`verification/run_v13.py` compares them with the v1.2 values kept in `results/v12_section9_values.json`). `cn_lemmas.json` and `co_envmin.json` are bit-for-bit those of v1.2. `co_envmin_opt.json` was re-run with the default 40 random starts (the v1.2 record used 10): the six minima agree to 4 × 10⁻¹⁰. `co_envpsi2.json` now keeps the full double-precision values (the v1.2 record kept four digits).
3. **The setting of Lemma CM-1.** The paper stated the lemma for "every focus"; the header of `src/cm_local_bound.py` said "b ≤ 1". No step of the proof uses a bound on b. The text now says "every focus (t₀, τ), whatever b = x/τ > 0", and the header is corrected. The bound Ψ₂ applies the lemma at a narrower focus, where b′ = b/(LKκ) may exceed one; restricting Ψ₂ to b′ ≤ 1 leaves the §IX A envelope values unchanged (a check of the follow-up work, numerical).
4. **Reproducing the 46 states.** `src/co_rtheta.py` now has the default argument 40, which gives the 46 states of the paper (40 random states and 6 wave packets); v1.2's default 30 gave 36 states.
5. **The top end ℓ = 960.8 in binary64** (note of 2026-11-11): the paper now says in §VI how the stored value 960.79999999999995452… is covered (the re-runs of `src/top_end.py`, made for v1.2, with the top end moved to the next binary64 number; `run_v13.py` checks the six stored runs).

6. **The stated range of Lemma BE** (known issue of 2026-11-18, found in an external review): the lemma is now stated for σ ∈ [10⁻³⁰⁰⁰, e^{−55.26}], that is ℓ ∈ [55.26, 6907.8], the range of Theorem H1. Up to v1.2 it was stated for σ ≤ 10⁻²⁴, which leaves out ℓ ∈ [55.26, 24 ln 10) (24 ln 10 = 55.262…). The proof uses an upper bound on σ in three places (ρ_W, the tail piece, and the last line of the lemma); `src/dw_bedge_low.py` re-checks them at σ = e^{−55.26} and at the binary64 value of 55.26: ρ_W ≤ 0.0801, κ ≥ 1.1 × 10⁻⁵ (the relative change of κ from σ = 10⁻²⁴ is 8 × 10⁻⁴⁶), the σ² term of ρ_W is 1.2 × 10⁻⁴³, the tail ratio is below 10^{−4.3×10⁴⁴}, and κ√ℓ/σ ≥ 8.16 × 10¹⁹. The text of §VII, App. D 5, the abstract ("σ from about 10⁻²⁴") and §VI are adjusted.

The new and changed files: `src/dw_bedge_low.py`, `results/dw_bedge_low.json` (new in this item), `src/dv_m2.py`, `src/dv_interval.py`, `src/dv_section9.py` (new); `src/co_remainders.py`, `src/co_rtheta.py`, `src/cn_lemmas.py`, `src/cm_local_bound.py` (changed; the exact changes are listed in `results/v13_code_changes.json`); `results/dv_m2.json`, `results/v12_section9_values.json`, `results/v13_code_changes.json` (new); the re-run results `results/co_rtheta.json` (with its log `co_rtheta_v13.log`), `co_remainders.json`, `cn_lemmas.json`, `co_envpsi2.json`, `co_envmin.json`, `co_envmin_opt.json`; `verification/run_v13.py` (new) and `verification/run_all.py` (calls `run_v13.py` instead of `run_v12.py`); `paper/sections/` (v1.3 text).

## What v1.2 fixes

Two steps of the proof of Theorem H1 were incomplete as written in v1.0 and v1.1. Both were found by the author's own audit while working on the narrow focus, were listed here as known issues as soon as they were found (2026-11-07 and 2026-11-10), and are closed in v1.2 by re-running the certificates. **The statement and the range of Theorem H1 are unchanged.** The code of v1.1 is not changed (only the text of one docstring in `src/cl_link.py`); the new runs wrap the published acceptance functions from outside.

1. **The end point G → 0⁺.** The cells of Lemma R\* are accepted at the two ends of the brightness G, because the bound is concave in G. One term of these tables, τν₂Q̄/(16(1+z)), does not contain G and enters only for G > 0, so the lower end is the limit G → 0⁺, not G = 0 alone. Every certificate of Theorem H1 has been re-run with this limit (`src/ct_v11_G0plus.py` → `results/ct_v11_*_G0plus*.jsonl`, `results/ct_summary.json`): the Lemma R\* boxes, the wedge and the μ tails of certificate (ii) (92 columns each), the 17 windows of certificate (iii) and the 987 tail columns of certificate (i). Every box passes, every cell count equals the published one, and 2.25 × 10⁸ evaluations at intermediate G found no negative value.
2. **The scaling step.** The certificates are run with the sandwich factors of the main terms set to one, and the factors are restored afterwards. With a lower factor λ < 1 the bound actually needed is T(λ) = P_bg + λ(P₀ + P₂) − N − λγσ(d_bg + d₂), and T(λ) = λT(1) + (1 − λ)(P_bg − N); v1.0 and v1.1 had used −N ≥ −λN, which fails for λ < 1. T is affine in λ, and every window of Theorem H1 has λ ≥ 4/5 (η ≤ Q̄ = 1/4) and an upper factor Υ ≥ 1, so T(1) ≥ 0 (the published certificate) and T(2/3) ≥ 0 suffice. The second certificate, certificate (ii) with the factors (2/3, 2/3, 1) and the target (2/3)γ, closes on all 92 + 92 + 92 columns, with G → 0⁺ in every box (`src/da_scaling.py` → `results/da_cl_w*.jsonl`, `da_wedge_w*.jsonl`, `da_tails_w*.jsonl`, `da_scaling.json`): 2.57 × 10⁹ cells, smallest normalised margin 2.4 × 10⁻¹⁰, and no negative value among 2.85 × 10⁸ evaluations at intermediate G. With the factors 2/3 but the target not lowered the first column does not close (`results/da_power.json`), so the run is not vacuous; the wrapped acceptance with the factors one is the published one bit for bit (`results/da_equiv.json`). Certificates (i) and (iii) do not use this step.

`verification/run_v12.py` (called by `run_all.py` up to v1.2; from v1.3 on, by `run_v13.py` on a copy with the v1.3 code changes undone) checks that the v1.1 code is unchanged, the stored re-runs, the algebra of the scaling step, and re-runs three columns of the second certificate.

## What v1.1 fixes

**Coverage of the computer-assisted proofs (affected the proof of Theorem H1 in v1.0).** An audit after v1.0 found three places where the union of the boxes did not cover the set of states that the certificates talk about. In none of them is there a sign that the inequality fails, and all three are closed in v1.1; every box of v1.0 is kept.

1. **The wedge of certificate (ii).** The v1.0 boxes (`src/cl_link.py`, 92 columns of width 10 in ℓ = ln(1/σ)) start at b = 10⁴σ² measured with the largest σ of each column, while Theorem H1 and Lemma BE split the states with the state's own σ. The wedge 10⁴σ² ≤ b < 10⁴σ_top² was covered by neither. v1.1 covers it by boxes in (μ, r, β, ℓ) with β = b/σ² and ℓ itself subdivided (`src/cq_core.py`, `results/cq_wedge.jsonl`): all 92 columns close, 4.44×10⁸ boxes. With β boxes but ℓ fixed per column the first column does not close (`results/cp_link_cells.jsonl`), because in the wedge the main balance and the largest cost both grow like σ.
2. **The μ tails of certificates (i) and (ii).** The v1.0 boxes covered μ ∈ [e⁻³⁰, μ_max·e⁶⁰] only, while μ = xS₂/S₀ lies anywhere in [x, ∞). v1.1 adds tail boxes reaching 0 and ∞ (`results/cq_tails_ii.jsonl`, `results/cq_tails_i_w0.jsonl`); the lower tail uses that, for fixed A, the bound is affine in G and in G/μ, so its minimum over the tail is at a corner.
3. **The smallest b of certificate (i).** Its claim is b ∈ (0, 1], but the boxes stopped at b = e^{−2100−3ℓ}, and at σ³ for the two finite μ tails. v1.1 adds a strip down to b = 0. Theorem H1 does not use certificate (i).

A random audit of 10⁵ admissible states (5×10⁴ for each claim of certificate (i)), weighted towards every edge, finds no uncovered state in the corrected tables, while it finds the gaps in the v1.0 tables and in eleven deliberately damaged tables (`src/cq_cover.py`). The new tail and wedge functions are generated from the v1.0 term table by a fixed list of textual replacements (`src/cq_make.py`) and agree with the v1.0 functions bit for bit on 42 000 boxes; ten deliberately wrong versions of the new steps are all broken (`src/cp_checks.py`, `src/cq_checks.py`). The monotonicity of the τ factors, checked on [20, 400] in v1.0, is now proved for every τ ≥ 20 (App. A). **The statement and the range of Theorem H1 are unchanged.**

**The narrow focus (τ < 20, Sec. IX A).** Two statements of v1.0 were wrong and are corrected:

1. The conjectured local Z3′, P_θ ≤ C Π*(μ_θ) with one universal constant C, is false: a coherent packet of high harmonics at a node of the fundamental gives P_θ/Π*(μ_θ) = 26.9 at x = 0.003 and 71.4 at x = 10⁻⁶. The value 2.13 of v1.0 came from a search whose packets had small raw coefficients. The correct form, P_θ ≤ c μ_θ e^{−μ_θ}, is proved from positivity alone (Lemma CM-1).
2. The local form of Lemma R does not keep the algebraic form of Lemma R: the term (1 − e)W carries the spread κ of the slope energy in the focus (Lemma CN-2).

v1.1 adds the proved Lemmas CM-1, CN-1 (two foci, with the joint brightness bounds Ψ₂ and Ψ₂*), CN-2 and R_θ (the local route with all remainders explicit), and a numerical envelope that exceeds the target T₀ at every ℓ of the table. The states with τ < 20 are **not** closed: the cell proof of this envelope is the next step.

## Known issues of v1.0 (kept for the record)

The three coverage gaps and the two narrow-focus statements above were listed here as known issues of v1.0 on 2026-11-04 and 2026-11-05, and are fixed in v1.1. In v1.0, Theorem H1 was proved only for the states with b ≥ 10⁴σ_top² in each column and e⁻³⁰ ≤ μ ≤ μ_max·e⁶⁰, and for all states with b < 10⁴σ² (Lemma BE).

## Known issues of v1.1 (fixed in v1.2, kept for the record)

Two known issues of v1.1 were listed here on 2026-11-07 (the end point G → 0⁺: a missing check, not a failed one) and on 2026-11-10 (the scaling step: the condition P_bg ≥ N, or T(λ) ≥ 0, was checked only at 8 × 10⁵ random states). Until v1.2, Theorem H1 was to be read as proved up to these checks. Both are closed in v1.2; see "What v1.2 fixes" above.

## Known issues of v1.2 (fixed in v1.3, kept for the record)

The entries below were listed here on 2026-11-11, 2026-11-15 and 2026-11-18. They are fixed in v1.3; see "What v1.3 fixes" above. The text is kept as it was.

### Note on v1.2 (2026-11-11): the binary64 value of the top end ℓ = 960.8

A review of v1.2 asked whether the decimal end points of Theorem H1 are inside the boxes once they are converted to binary64. The lower end is: 55.26 is stored as 55.25999999999999801…, slightly below it. The upper end was not: 960.8 is stored as 960.79999999999995452…, so the last column of certificate (ii) stopped 4.5 × 10⁻¹⁴ below the stated end. This was a missing check at the edge, not a failed one. The last column of each part of certificate (ii) (Lemma R\* boxes, wedge, μ tails) has been re-run with its top end moved outward to the next binary64 number, 960.8000000000001, both as published (factors one, with G → 0⁺) and as the second certificate of the scaling step (factors 2/3): all six runs close, each with exactly the number of cells of the same column with the old end point (published, and second certificate) (`src/top_end.py` → `results/top_end_lam1.jsonl`, `results/top_end_lam2_3.jsonl`). So the stated range ℓ ∈ [55.26, 960.8] is covered. The paper will say this in its next version.

### Known issue of v1.2 (2026-11-15): a rounded constant in `src/cn_lemmas.py`

In `src/cn_lemmas.py` (the numerical bound Ψ₂ from Lemma CN-1, used by the §8 scripts `co_envmin.py`, `co_envpsi2.py`, `co_remainders.py` and `co_rtheta.py`), the coefficient of the b-term is written as `0.0325`. The value that follows from Lemma CM-1 is (√2/e)/16 = 0.0325163…, so the code rounds it down by 0.05 %.

- **Theorem H1 is not affected.** Its proof and its checks (`verification/run_v12.py`) do not use this file. The text of the paper states Lemma CM-1 with the exact coefficient √2/e.
- **The numbers in §8 (open problems, marked as numerical) do not change.** The author re-ran the five §8 scripts on a copy of this repository with the coefficient replaced by √2/(16e): every number that `verification/hp_numbers.py` checks against the text is unchanged in the printed digits, and the raw values are bit-for-bit identical. (At the states that matter in §8, b is very small, so the b-term does not reach the double-precision sum.)
- The constant will be corrected in the next version of the paper and repository.

Note on reproducing §8: the saved result of `co_rtheta.py` ("46 states" in the text) was made with `python src/co_rtheta.py 40` (40 states + 6 wave packets). Running it with the default argument gives 36 states; the range of ratios (27–94 %) is the same.

### Known issue of v1.2 (2026-11-18): a constant of Lemma R_θ taken from a grid

In `src/co_remainders.py`, the remainder N16 of Lemma R_θ (§8, the local route) uses the constant M₂ = max over u ≥ 0 of u e^{−u/2}(1 + u cosh √(yu)). The function `M2_const` takes this maximum on a grid of u in [0, 400] and multiplies it by (1 + 10⁻⁶). This is a numerical value, not a proved bound, although the lemma is marked as proved.

- **Theorem H1 is not affected.** Lemma R_θ belongs to the local route of §8 (open problems), which Theorem H1 does not use.
- The grid value is very close to the true maximum: for y = 0 the maximum is attained at u = (3 + √17)/2 and equals 2.73760…, and in the follow-up work an interval-arithmetic bound M₂ ≤ 2.74484 has been proved for every y ≤ 10⁻⁴.
- In the next version of the paper and repository, the grid value will be replaced by a proved bound, and the range of y used in §8 will be stated.

### Known issue of v1.0–v1.2 (2026-11-18): the stated range of Lemma BE

Lemma BE (§VII, the small-b edge) is stated for σ ∈ [10⁻³⁰⁰⁰, 10⁻²⁴]. Theorem H1 is stated for ℓ ≥ 55.26, and 24 ln 10 = 55.26204…, so the states with ℓ ∈ [55.26, 24 ln 10), that is σ ∈ (10⁻²⁴, e^{−55.26} = 1.00204 × 10⁻²⁴], are outside the stated hypothesis of the lemma that the proof of Theorem H1 uses for b ≤ 10⁴σ².

- **This is not a counterexample, and the statement of Theorem H1 does not change.** The proof of Lemma BE uses an upper bound on σ in only three places, and each has an enormous margin: the lemma gives a lower bound of order 8 × 10¹⁹ for g√ℓ, where about 1 is needed.
- The fix is to state the hypothesis as σ ≤ e^{−55.26} and to re-check those three places at that end (`src/dw_bedge_low.py`); this is done in v1.3.
- Found by an external review of the follow-up work.

## How the computer-assisted proofs are accepted

The cells are evaluated in double precision (IEEE 754 binary64, round to nearest). A cell is accepted only when its margin is at least **ten times a rounding bound derived in App. E** of the paper from the operations of the cell. For the full bound, every term is moved in the safe direction before equal powers of σ are netted, so cancellation is covered. The derived bounds are checked against long-double re-evaluation (`src/ck_rounding.py`). The bound assumes that `exp`, `log`, `log1p`, `expm1` and `sqrt` are within two ulp; `src/ulp_check.py` measures the library actually used against a 60-digit evaluation (largest error 0.71 ulp with NumPy 2.5.3 on x86-64 Linux). The older records (`results/cc_cells.jsonl`, `cd_cells.jsonl`, `ce_link_cells.jsonl`, `cf_windows.jsonl`) are kept for the comparison of box counts.

## Verify

```bash
pip install -r requirements.txt
python verification/run_all.py              # under 30 minutes: every number of the paper, a sample of every certificate, the random covering audit
python verification/run_all.py --recompute  # every column and window (many CPU hours)
```

`run_all.py` re-derives every number of the manuscript from `results/` and checks that it appears in `paper/sections/`, checks every stored certificate, recomputes a sample of each one with the rounding-covering rule (including the μ tails and the wedge), and runs the random covering audit with its damaged tables. No GPU is used.

## What is where

| result | status | script → results |
|---|---|---|
| Lemma R (lower bound for every state) | proof | `src/ca_reduction.py` |
| Lemma A4* (the remainder a₄ made relative) | proof | `src/ce_a4.py` → `results/ce_a4.json` |
| Lemma H, τ cells (factor ≥ 0.99977 on [20, 80]) | proof / computer-assisted | `src/cb_tau.py` → `results/cb_tau_cells.json` |
| Lemma S, S′, Corollary S (split sandwich) | proof | `src/cd_split.py`, `src/ce_split_map.py` |
| main reduced inequality, ℓ ≤ 6907.8 | computer-assisted | `src/cc_cells.py`, `src/cd_cells.py` with the rule of `src/cl_cells.py` → `results/cl_main_cells.jsonl`; μ tails and b → 0: `src/cq_core.py` → `results/cq_tails_i_w0.jsonl` |
| full Lemma R* bound, τ ≥ 20, b ≥ 10⁴σ², ℓ ≤ 960.8 | computer-assisted | `src/cl_link.py` → `results/cl_link_cells.jsonl`; μ tails and wedge: `src/cq_core.py` (generated `src/cq_phi.py`) → `results/cq_tails_ii.jsonl`, `results/cq_wedge.jsonl` |
| covering of the states by the boxes and Lemma BE | proof (line by line) + random audit | `src/cq_cover.py` → `results/cq_cover.json`; checks of the new steps `src/cp_checks.py`, `src/cq_checks.py` |
| Lemma BE (every state with b < 10⁴σ²) | proof (constants checked) | `src/ch_bedge.py` → `results/ch_bedge.json` |
| v1.2: end point G → 0⁺ in every certificate of Theorem H1 | computer-assisted (re-run) | `src/ct_v11_G0plus.py` → `results/ct_v11_*_G0plus*.jsonl`, `results/ct_summary.json` |
| v1.2: second certificate of the scaling step (factors 2/3, target (2/3)γ) | computer-assisted | `src/da_scaling.py` → `results/da_*_w*.jsonl`, `results/da_scaling.json`, `da_equiv.json`, `da_power.json` |
| 17 thin-block windows | computer-assisted | `src/cf_windows.py` with `src/cl_link.py` → `results/cl_windows.jsonl` |
| counterexamples to conjecture C | numerical (exact sums) | `src/cc_counter.py` → `results/cc_counter.json` |
| local corrections (L1)–(L3), Lemma A* | proof | `src/ci_rtheta.py` → `results/ci_rtheta.json` |
| local Z3′ with one constant (v1.0 conjecture) | false (counterexample) | `src/cm_counter.py` → `results/cm_counter.json` |
| v1.3: the constant M₂ of Lemma R_θ (N16) | proof (y ≤ 10⁻²⁰) / computer-assisted (interval arithmetic) | `src/dv_m2.py`, `src/dv_interval.py` → `results/dv_m2.json` |
| v1.3: Lemma BE at the upper end σ = e^{−55.26} | proof; constants checked in 60-digit decimal | `src/dw_bedge_low.py` → `results/dw_bedge_low.json` |
| v1.3: the top end ℓ = 960.8 in binary64 (runs made for v1.2) | computer-assisted (re-run of the last columns) | `src/top_end.py` → `results/top_end_lam1.jsonl`, `results/top_end_lam2_3.jsonl` |
| Lemmas CM-1, CN-1 (Ψ₂, Ψ₂*), CN-2, R_θ | proof (+ numerical checks) | `src/cm_local_bound.py`, `src/cm_checks.py`, `src/cn_lemmas.py`, `src/co_rtheta.py`, `src/co_remainders.py` |
| joint brightness bound search | numerical (search with power checks) | `src/cn_search.py` → `results/cn_search.json` |
| local envelope with Ψ₂* above T₀ | numerical | `src/co_envmin.py` → `results/co_envmin_opt.json` |
| single columns with Lemma A* (τ ≥ 8 to 12) | computer-assisted (single columns) | `src/ci_link.py` → `results/ci_link_cells.jsonl` |
| rounding bounds (App. E), τ cells and Lemma BE in decimal | derivation / checks | `src/cl_cells.py`, `src/cl_link.py`, `src/ck_rounding.py`, `src/cl_tau_be.py` |
| App. D line by line, in exact fractions and 60-digit decimal | checks | `src/ck_details.py` → `results/ck_details.json` |

| folder | content |
|---|---|
| `src/` | lemmas, certificates and checks (Python 3, NumPy, SciPy) |
| `results/` | stored results and certificates |
| `verification/` | `run_all.py`, `hp_numbers.py` (the numbers of the paper), `run_v12.py` (the checks of v1.2) and `run_v13.py` (the checks of v1.3; it runs `run_v12.py` on a copy with the v1.3 code changes undone) |
| `paper/` | LaTeX source of the paper v1.3 (REVTeX 4.2) and the PDFs of v1.0, v1.1 and v1.2 (the PDF of v1.3 is added on publication) |

## Related work by the author

- T. Namba, *How fine can a hidden ripple be? A two-sided logarithmic law for the resolution gap of the arrow of time*, Zenodo (2026), [doi:10.5281/zenodo.23856796](https://doi.org/10.5281/zenodo.23856796) (paper F, v1.5; code [resolution-gap-two-sided](https://github.com/neginyan/resolution-gap-two-sided)). Paper F and its repository are not changed by this work.

## Use of AI assistants

The study was carried out with AI assistants: Anthropic Claude (planning, proofs, computations, code, verification and drafting), Google Gemini and OpenAI ChatGPT (review of the manuscript). Their roles are described in the acknowledgments of the paper. The author directed the research and takes full responsibility for it.

## Citation

Please cite the paper ([doi:10.5281/zenodo.23932698](https://doi.org/10.5281/zenodo.23932698)) and this repository; see `CITATION.cff`.
The code is archived at Zenodo: [doi:10.5281/zenodo.23932577](https://doi.org/10.5281/zenodo.23932577) (all versions; v1.0.0 is [10.5281/zenodo.23932578](https://doi.org/10.5281/zenodo.23932578); v1.1.0 is [10.5281/zenodo.23993961](https://doi.org/10.5281/zenodo.23993961); v1.2.0 is [10.5281/zenodo.24115114](https://doi.org/10.5281/zenodo.24115114); the DOI of v1.3.0 is added on publication).

## License

MIT © 2026 Taishi Namba
