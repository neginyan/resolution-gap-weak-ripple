# resolution-gap-weak-ripple

**Beyond few harmonics: the lower side of the logarithmic law of the resolution gap for windows with a weak ripple**

Code, stored results, computer-assisted certificates and a verification script for the paper by Taishi Namba (2026).

**Paper:** T. Namba, *Beyond few harmonics: the lower side of the logarithmic law of the resolution gap for windows with a weak ripple*, preprint, Zenodo (2026), [doi:10.5281/zenodo.23932698](https://doi.org/10.5281/zenodo.23932698) (all versions; v1.0 is [10.5281/zenodo.23932699](https://doi.org/10.5281/zenodo.23932699); v1.1 is [10.5281/zenodo.23993985](https://doi.org/10.5281/zenodo.23993985)).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23932577.svg)](https://doi.org/10.5281/zenodo.23932577)
[![verify](https://github.com/neginyan/resolution-gap-weak-ripple/actions/workflows/verify.yml/badge.svg)](https://github.com/neginyan/resolution-gap-weak-ripple/actions/workflows/verify.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## The result in brief

This is the follow-up to [paper F](https://doi.org/10.5281/zenodo.23856796) (code: [resolution-gap-two-sided](https://github.com/neginyan/resolution-gap-two-sided)). Paper F proved the logarithmic law of the resolution gap of the arrow of time, with the sharp constant, for periodic windows with at most twenty harmonics or with damping exponent K ≥ 0.003. Here the **lower side** is extended to windows with **any number of harmonics** whose ripple is weak in the sup norm or poor in low-frequency energy.

- **Theorem H1.** For every state whose focus spans at least about 0.7 ripple periods (τ = x/b ≥ 20), and for ℓ = ln(1/σ) ∈ [55.26, 960.8] (σ from 10⁻²⁴ down to about 10⁻⁴¹⁷), g√ℓ ≥ (1 − ε̄)(λ/Υ)(1 + 12.1/ℓ + 2.5×10⁻⁴). The coefficient one follows, up to 1 − ε̄, for every window with η ≤ 1/ℓ, for a family of intermediate windows, and for seventeen thin-block windows certified one by one. States with b < 10⁴σ² obey a far stronger bound (Lemma BE). In v1.1 the boxes of the computer-assisted proofs cover every state of the theorem; this is shown line by line and checked on 10⁵ random states.
- **Three tools.** A remainder of the curvature step made relative to the ripple (Lemma A4*); a split sandwich in which a low-frequency block enters through its energy, not its sup (Lemma S); a mean-square bound of the curvature cross term that closes the edge of very wide focus (Lemma BE).
- **A conjecture of our own is false.** The working hypothesis ηe^{μ/2} ≤ C fails on explicit windows (up to 10³² at σ = 10⁻³⁰⁰); the coefficient one survives on all of them.
- **The road to an unconditional theorem.** The paper lists what is still open, from the nearest to the farthest: states with a narrow focus (τ < 20), ℓ > 960.8 and the limit σ → 0, windows outside the covered families, the class constant ε̄, and the upper side. The narrow focus is now described with proved lemmas and a numerical envelope, but it is not closed; see "What v1.1 fixes" below.

Every statement in the paper carries its status: proof, computer-assisted proof, derivation, numerical, or conjecture (Table I of the paper).

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
| 17 thin-block windows | computer-assisted | `src/cf_windows.py` with `src/cl_link.py` → `results/cl_windows.jsonl` |
| counterexamples to conjecture C | numerical (exact sums) | `src/cc_counter.py` → `results/cc_counter.json` |
| local corrections (L1)–(L3), Lemma A* | proof | `src/ci_rtheta.py` → `results/ci_rtheta.json` |
| local Z3′ with one constant (v1.0 conjecture) | false (counterexample) | `src/cm_counter.py` → `results/cm_counter.json` |
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
| `verification/` | `run_all.py` and `hp_numbers.py` (the numbers of the paper) |
| `paper/` | LaTeX source of the paper v1.1 (REVTeX 4.2) and the PDFs of v1.0 and v1.1 |

## Related work by the author

- T. Namba, *How fine can a hidden ripple be? A two-sided logarithmic law for the resolution gap of the arrow of time*, Zenodo (2026), [doi:10.5281/zenodo.23856796](https://doi.org/10.5281/zenodo.23856796) (paper F, v1.5; code [resolution-gap-two-sided](https://github.com/neginyan/resolution-gap-two-sided)). Paper F and its repository are not changed by this work.

## Use of AI assistants

The study was carried out with AI assistants: Anthropic Claude (planning, proofs, computations, code, verification and drafting), Google Gemini and OpenAI ChatGPT (review of the manuscript). Their roles are described in the acknowledgments of the paper. The author directed the research and takes full responsibility for it.

## Citation

Please cite the paper ([doi:10.5281/zenodo.23932698](https://doi.org/10.5281/zenodo.23932698)) and this repository; see `CITATION.cff`.
The code is archived at Zenodo: [doi:10.5281/zenodo.23932577](https://doi.org/10.5281/zenodo.23932577) (all versions; v1.0.0 is [10.5281/zenodo.23932578](https://doi.org/10.5281/zenodo.23932578); v1.1.0 is [10.5281/zenodo.23993961](https://doi.org/10.5281/zenodo.23993961)).

## License

MIT © 2026 Taishi Namba
