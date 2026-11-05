# resolution-gap-weak-ripple

**Beyond few harmonics: the lower side of the logarithmic law of the resolution gap for windows with a weak ripple**

Code, stored results, computer-assisted certificates and a verification script for the paper by Taishi Namba (2026).

**Paper:** T. Namba, *Beyond few harmonics: the lower side of the logarithmic law of the resolution gap for windows with a weak ripple*, preprint, Zenodo (2026), [doi:10.5281/zenodo.23932698](https://doi.org/10.5281/zenodo.23932698) (all versions; v1.0 is [10.5281/zenodo.23932699](https://doi.org/10.5281/zenodo.23932699)).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23932577.svg)](https://doi.org/10.5281/zenodo.23932577)
[![verify](https://github.com/neginyan/resolution-gap-weak-ripple/actions/workflows/verify.yml/badge.svg)](https://github.com/neginyan/resolution-gap-weak-ripple/actions/workflows/verify.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## The result in brief

This is the follow-up to [paper F](https://doi.org/10.5281/zenodo.23856796) (code: [resolution-gap-two-sided](https://github.com/neginyan/resolution-gap-two-sided)). Paper F proved the logarithmic law of the resolution gap of the arrow of time, with the sharp constant, for periodic windows with at most twenty harmonics or with damping exponent K ≥ 0.003. Here the **lower side** is extended to windows with **any number of harmonics** whose ripple is weak in the sup norm or poor in low-frequency energy.

- **Theorem H1.** For every state whose focus spans at least about 0.7 ripple periods (τ = x/b ≥ 20), and for ℓ = ln(1/σ) ∈ [55.26, 960.8] (σ from 10⁻²⁴ down to about 10⁻⁴¹⁷), g√ℓ ≥ (1 − ε̄)(λ/Υ)(1 + 12.1/ℓ + 2.5×10⁻⁴). The coefficient one follows, up to 1 − ε̄, for every window with η ≤ 1/ℓ, for a family of intermediate windows, and for seventeen thin-block windows certified one by one. States with b < 10⁴σ² obey a far stronger bound (Lemma BE). *(The computer-assisted part has a coverage gap found after v1.0; see the first known issue below.)*
- **Three tools.** A remainder of the curvature step made relative to the ripple (Lemma A4*); a split sandwich in which a low-frequency block enters through its energy, not its sup (Lemma S); a mean-square bound of the curvature cross term that closes the edge of very wide focus (Lemma BE).
- **A conjecture of our own is false.** The working hypothesis ηe^{μ/2} ≤ C fails on explicit windows (up to 10³² at σ = 10⁻³⁰⁰); the coefficient one survives on all of them.
- **The road to an unconditional theorem.** The paper lists what is still open, from the nearest to the farthest: states with a narrow focus (τ < 20), ℓ > 960.8 and the limit σ → 0, windows outside the covered families, the class constant ε̄, and the upper side. Two statements about the narrow focus in Sec. IX A of v1.0 are wrong, and certificate (ii) has a coverage gap; see the known issues below.

Every statement in the paper carries its status: proof, computer-assisted proof, derivation, numerical, or conjecture (Table I of the paper).

## Known issue (2026-11-05, after v1.0): coverage gaps in certificates (i) and (ii) — fixed in the lab, to appear in v1.1

**This one affects the proof of Theorem H1.** The boxes of certificate (ii) (`src/cl_link.py`, 92 columns of width 10 in ℓ = ln(1/σ)) start at b = 10⁴σ² measured with the *largest* σ of each column, while Theorem H1 and Lemma BE split the states at b = 10⁴σ² with the state's *own* σ. Inside each column, the states with 10⁴σ² ≤ b < 10⁴σ_top² (σ_top = e^{−ℓ_a}, the top of the column) are therefore covered by neither the certificate nor Lemma BE. For these states Theorem H1 is, at present, not proved.

There is no sign that the inequality fails there. Covering the wedge needs a different box layout: with columns of width 1 or 2 whose lower edge in b is taken at the bottom σ of the column, the column at ℓ = 55.26 closes; at ℓ ≈ 960 a column of width 1 did not close within the budget, because the margin there is thin. The planned fix for v1.1 is to use β = b/σ² instead of b as a box variable, so that the edge b = 10⁴σ² becomes a face of the boxes, and to subdivide ℓ inside the wedge, where the balance of the main terms is sensitive to σ (a first test with β boxes and columns of width 10 did not close; columns of width 1 do).

A second coverage gap was found in the same audit. The boxes of certificates (i) and (ii) cover μ ∈ [e⁻³⁰, μ_max·e⁶⁰] only (μ_max = max(4ℓ_b, 50)), while the states of the theorem have μ = xS₂/S₀ anywhere in [x, ∞). The states with smaller or larger μ are not covered by the published boxes. Point evaluations there are all positive (smallest normalised value 0.59 on 12 000 random states), and v1.1 will add tail boxes reaching 0 and ∞.

A third, smaller gap of the same kind was found while closing these: the boxes of certificate (i) did not cover b < σ³ in its two finite μ tails, nor b < e^{−2100−3ℓ} for any μ, although its claim is b ∈ (0, 1]. Theorem H1 does not use certificate (i).

**Status (2026-11-06):** all three gaps are closed in the lab results that will become v1.1. The wedge is covered by boxes in (μ, r, β, ℓ) with β = b/σ² and ℓ itself subdivided (all 92 columns close), the μ tails by boxes reaching 0 and ∞ for certificates (i) and (ii), and the b tail of certificate (i) by a band down to b = 0. A random audit of 10⁵ admissible states finds no uncovered state in the corrected tables, while it does find the gaps in the v1.0 tables and in eleven deliberately damaged tables. The heaviest new column was recomputed independently with the same box count. Until v1.1 is published, Theorem H1 should be read as proved for the states with b ≥ 10⁴σ_top² in each column and e⁻³⁰ ≤ μ ≤ μ_max·e⁶⁰, and for all states with b < 10⁴σ² (Lemma BE). The seventeen-window certificate (iii), Lemma BE and paper F are not affected. The random audit will be part of the v1.1 verification.

## Known issue (2026-11-04, after v1.0)

Two statements about the states with a narrow focus (τ < 20) in Sec. IX A of v1.0, repeated in one sentence of the abstract, are wrong. They will be corrected in v1.1.
**Theorem H1, the three certificates, Lemma BE and Apps. D and E are not affected.** They concern the states with τ ≥ 20, where the remainder Θ − θ_τ is bounded by Lemma A4*, and they use Lemma Z3′ in its proved period-average form.

1. **The conjectured local Z3′, P_θ ≤ C Π*(μ_θ) with one universal constant C, is false.** A coherent packet of high harmonics placed at a node of the fundamental (it satisfies (P2), (P3) and positivity) gives P_θ/Π*(μ_θ) = 27 at x = 0.003 and 71 at x = 10⁻⁶, and the ratio grows with μ_θ. The value 2.13 quoted in v1.0 came from a search whose packets had small raw coefficients, and the statement that the bound tolerates C = 10 loses its meaning. The correct form is P_θ ≤ c μ_θ e^{−μ_θ}; a bound of this form follows from the positivity of the profile alone (lab stage CM, to appear in v1.1).
2. **The sentence that the local form of Lemma R "keeps the algebraic form of Lemma R" does not hold for the term (1 − e)W.** For a narrow focus this term is weighted by how the slope energy of the ripple is spread inside the focus, and it is smaller when the ripple is concentrated at the centre of the focus. So the narrow-focus states do not reduce to one local inequality, as the abstract of v1.0 says.

The states with τ < 20 remain open, as stated in Theorem H1.

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

## Citation

Please cite the paper ([doi:10.5281/zenodo.23932698](https://doi.org/10.5281/zenodo.23932698)) and this repository; see `CITATION.cff`.
The code is archived at Zenodo: [doi:10.5281/zenodo.23932577](https://doi.org/10.5281/zenodo.23932577) (all versions; v1.0.0 is [10.5281/zenodo.23932578](https://doi.org/10.5281/zenodo.23932578)).

## License

MIT © 2026 Taishi Namba
