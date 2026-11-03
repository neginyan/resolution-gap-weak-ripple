# resolution-gap-weak-ripple

**Beyond few harmonics: the lower side of the logarithmic law of the resolution gap for windows with a weak ripple**

Code, stored results, computer-assisted certificates and a verification script for the paper by Taishi Namba (2026).

[![verify](https://github.com/neginyan/resolution-gap-weak-ripple/actions/workflows/verify.yml/badge.svg)](https://github.com/neginyan/resolution-gap-weak-ripple/actions/workflows/verify.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## The result in brief

This is the follow-up to [paper F](https://doi.org/10.5281/zenodo.23856796) (code: [resolution-gap-two-sided](https://github.com/neginyan/resolution-gap-two-sided)). Paper F proved the logarithmic law of the resolution gap of the arrow of time, with the sharp constant, for periodic windows with at most twenty harmonics or with damping exponent K ≥ 0.003. Here the **lower side** is extended to windows with **any number of harmonics** whose ripple is weak in the sup norm or poor in low-frequency energy.

- **Theorem H1.** For every state whose focus spans at least about 0.7 ripple periods (τ = x/b ≥ 20), and for ℓ = ln(1/σ) ∈ [55.26, 960.8] (σ from 10⁻²⁴ down to about 10⁻⁴¹⁷), g√ℓ ≥ (1 − ε̄)(λ/Υ)(1 + 12.1/ℓ + 2.5×10⁻⁴). The coefficient one follows, up to 1 − ε̄, for every window with η ≤ 1/ℓ, for a family of intermediate windows, and for seventeen thin-block windows certified one by one. States with b < 10⁴σ² obey a far stronger bound (Lemma BE).
- **Three tools.** A remainder of the curvature step made relative to the ripple (Lemma A4*); a split sandwich in which a low-frequency block enters through its energy, not its sup (Lemma S); a mean-square bound of the curvature cross term that closes the edge of very wide focus (Lemma BE).
- **A conjecture of our own is false.** The working hypothesis ηe^{μ/2} ≤ C fails on explicit windows (up to 10³² at σ = 10⁻³⁰⁰); the coefficient one survives on all of them.
- **The road to an unconditional theorem.** The paper lists what is still open, from the nearest to the farthest: states with a narrow focus (τ < 20; everything reduces to one local inequality, which holds within a factor 2.13 on the windows tried while the bound tolerates 10), ℓ > 960.8 and the limit σ → 0, windows outside the covered families, the class constant ε̄, and the upper side.

Every statement in the paper carries its status: proof, computer-assisted proof, derivation, numerical, or conjecture (Table I of the paper).

## How the computer-assisted proofs are accepted

The cells are evaluated in double precision (IEEE 754 binary64, round to nearest). A cell is accepted only when its margin is at least **ten times a rounding bound derived in App. E** of the paper from the operations of the cell. For the full bound, every term is moved in the safe direction before equal powers of σ are netted, so cancellation is covered. The derived bounds are checked against long-double re-evaluation (`src/ck_rounding.py`). The older records (`results/cc_cells.jsonl`, `cd_cells.jsonl`, `ce_link_cells.jsonl`, `cf_windows.jsonl`) are kept for the comparison of box counts.

## Verify

```bash
pip install -r requirements.txt
python verification/run_all.py              # about 20 minutes: every number of the paper + a sample of every certificate
python verification/run_all.py --recompute  # every column and window (many CPU hours)
```

`run_all.py` re-derives every number of the manuscript from `results/` and checks that it appears in `paper/sections/`, checks every stored certificate, and recomputes a sample of each one with the rounding-covering rule. No GPU is used.

## What is where

| result | status | script → results |
|---|---|---|
| Lemma R (lower bound for every state) | proof | `src/ca_reduction.py` |
| Lemma A4* (the remainder a₄ made relative) | proof | `src/ce_a4.py` → `results/ce_a4.json` |
| Lemma H, τ cells (factor ≥ 0.99977 on [20, 80]) | proof / computer-assisted | `src/cb_tau.py` → `results/cb_tau_cells.json` |
| Lemma S, S′, Corollary S (split sandwich) | proof | `src/cd_split.py`, `src/ce_split_map.py` |
| main reduced inequality, ℓ ≤ 6907.8 | computer-assisted | `src/cc_cells.py`, `src/cd_cells.py` with the rule of `src/cl_cells.py` → `results/cl_main_cells.jsonl` |
| full Lemma R* bound, τ ≥ 20, b ≥ 10⁴σ², ℓ ≤ 960.8 | computer-assisted | `src/cl_link.py` → `results/cl_link_cells.jsonl` |
| Lemma BE (every state with b < 10⁴σ²) | proof (constants checked) | `src/ch_bedge.py` → `results/ch_bedge.json` |
| 17 thin-block windows | computer-assisted | `src/cf_windows.py` with `src/cl_link.py` → `results/cl_windows.jsonl` |
| counterexamples to conjecture C | numerical (exact sums) | `src/cc_counter.py` → `results/cc_counter.json` |
| local corrections (L1)–(L3), Lemma A* | proof | `src/ci_rtheta.py` → `results/ci_rtheta.json` |
| local Z3′ for τ < 20 | conjecture / numerical | `src/ci_local.py`, `src/ci_envloc.py`, `src/ci_nop2.py` |
| single columns with Lemma A* (τ ≥ 8 to 12) | computer-assisted (single columns) | `src/ci_link.py` → `results/ci_link_cells.jsonl` |
| rounding bounds (App. E), τ cells and Lemma BE in decimal | derivation / checks | `src/cl_cells.py`, `src/cl_link.py`, `src/ck_rounding.py`, `src/cl_tau_be.py` |
| App. D line by line, in exact fractions and 60-digit decimal | checks | `src/ck_details.py` → `results/ck_details.json` |

| folder | content |
|---|---|
| `src/` | lemmas, certificates and checks (Python 3, NumPy, SciPy) |
| `results/` | stored results and certificates |
| `verification/` | `run_all.py` and `hp_numbers.py` (the numbers of the paper) |
| `paper/` | LaTeX source of the paper (REVTeX 4.2) and the PDF of v1.0 |

## Related work by the author

- T. Namba, *How fine can a hidden ripple be? A two-sided logarithmic law for the resolution gap of the arrow of time*, Zenodo (2026), [doi:10.5281/zenodo.23856796](https://doi.org/10.5281/zenodo.23856796) (paper F, v1.5; code [resolution-gap-two-sided](https://github.com/neginyan/resolution-gap-two-sided)). Paper F and its repository are not changed by this work.

## Use of AI assistants

The study was carried out with AI assistants: Anthropic Claude (planning, proofs, computations, code, verification and drafting), Google Gemini and OpenAI ChatGPT (review of the manuscript). Their roles are described in the acknowledgments of the paper. The author directed the research and takes full responsibility for it.

## License

MIT © 2026 Taishi Namba
