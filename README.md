# Simulating analogue Hawking radiation with cross-review between two AIs

**An acoustic black hole, a spurious result, and how its origin was identified**

Jordi Torres Ribas · v1.5 · 18 November 2026

A numerical solver for analogue Hawking radiation in a one-dimensional Bose-Einstein
condensate flowing faster than sound (an "acoustic black hole"). Halfway through the
project the solver produced a striking result — the radiation temperature collapsing to
zero when the short-scale physics was modified — that turned out to be a modelling error.
This repository contains the solver, the validations, the investigation of that spurious
result, and every intermediate result, including the negative ones.

- **Write-up for a general/technical audience:** [`docs/WRITEUP_EN.md`](docs/WRITEUP_EN.md)
  (Spanish: [`docs/WRITEUP.md`](docs/WRITEUP.md))
- **Scientific document of the frozen benchmark:** [`docs/FP_BENCHMARK.md`](docs/FP_BENCHMARK.md)

Most working documents are in Spanish; this README and the English write-up are the entry
points.

## Main results

| Result | Value | Reproduce |
|---|---|---|
| Benchmark v1.0 vs Finazzi & Parentani, PRD 83, 084010 (2011) | T_up/(κ/2π) = 0.957 → 0.985 for w = 2 → 8, approaching 1 as the horizon gets smoother. With FP's broadened-horizon κ̄ the agreement is ≤ 0.75 %, but only for their window width d = 1.58 (fitted to their profiles); with d = 1 or 2 the maximum deviation is 2.4 % or 3.4 % (`docs/FP_BENCHMARK.md` §5b) | `python3 benchmarks/fp_v1/FP_check_reproduction.py` (~30 s) |
| External validation vs Larré, Recati, Carusotto & Pavloff, PRA 85, 013621 (2012) | leading coefficient: relative error 2.4–2.9·10⁻⁵; secondary: < 10⁻³ | `python3 residual_anatomy/gates/gate_E1_lrcp.py` (~1 min) |
| E1-bis: same validation with the ω → 0 extrapolation controlled, criteria pre-registered in commit `001a7ac` before the run | PASS: leading coefficient (A8) to 5.4·10⁻⁹ (m_u = 0.6) and 4.2·10⁻⁸ (m_u = 0.7); secondary (A14) to 7.9·10⁻⁷ and 3.1·10⁻⁶ (`residual_anatomy/E1_RESULT.md`) | `python3 residual_anatomy/gates/gate_E1bis_lrcp.py --selftest` (~15 s; the full gate refuses to overwrite its results) |
| The "collapse" was an artefact (inconsistent background) | with a consistent background T_eff/T_H varies smoothly 0.979 → 0.892 for λ = 0.25 → 4. By the exact equivalence (λ, w) ≡ (1, w/√λ) (next row), this is the effect of the horizon width w/√λ = 4 → 1 in healing lengths, **not** of short-scale (UV) physics | `python3 residual_anatomy/diag_T2_qp_consistent.py` (~1 min) |
| Goldstone (symmetry) test identifying the defect | measured / predicted residual = 1.0000 | `python3 residual_anatomy/diag_goldstone.py` (~7 min) |
| With a consistent background, λ is an exact rescaling: (λ, w) ≡ (1, w/√λ) | the λ-solver and the v1.0 solver on the rescaled profile agree to ≤ 3.2·10⁻⁶ (rounding floor of the plain sparse solve ~7·10⁻⁶) in 90 points; with an equilibrated, iteratively refined solve of both sides, to ≤ 3.5·10⁻⁹ (`residual_anatomy/diag_conditioning_equiv.py`, ~2 min). So λ does **not** isolate short-scale physics | `python3 residual_anatomy/diag_scaling_equivalence.py` (~30 s) |

**Scope.** No new physics is claimed: the robustness of the radiation against short-scale
modifications was already established (Macher, Finazzi, Parentani and others). The
artefact is shown to be one within the tested scope (one flow profile, one range of λ).
The external validation covers the basic solver (λ = 1). See §7 of the write-up.

> **Audit notice (2026-11-17).** `docs/CORRELATIONS.md` §4, §5, §7 and §8 (robustness of the
> density-density "tongue" against λ, and the w × λ grid) were computed with the defective λ
> operator and **have no physical validity**; their conclusions are withdrawn. §1–3 (λ = 1) are
> not affected. Details at the top of that document.

## Quick start

```bash
# Python 3.9.6 (see .python-version); pinned dependencies:
python3 -m pip install -r requirements.txt

# Always run from the repository root:
python3 benchmarks/fp_v1/FP_check_reproduction.py

# All validations with pass/fail exit codes (~2 min):
python3 tests/run_all.py
```

`gate_E1_lrcp.py` refuses to overwrite existing results: move `out/e1_lrcp_*` aside first
to re-run it. The other scripts overwrite their outputs in `out/`.

## Layout

| Folder | Contents |
|---|---|
| `solver/` | Core solver (frozen v1.0): GPE background, finite-difference BdG operator, S-matrix, exact modal matching, dispersive (λ) operator, symbolic κ |
| `tests/` | Solver validations (soliton, Thomas-Fermi, operator checks) |
| `benchmarks/fp_v1/` | Finazzi-Parentani benchmark v1.0 and its controls |
| `reference/` | Frozen reference table for the benchmark |
| `methods/` | Validation chain of the methods, including failed attempts (see `docs/REPRODUCE.md` §4) |
| `correlations/` | Density-density correlation line ("Rumbo #1", `docs/CORRELATIONS.md`) |
| `residual_anatomy/` | λ / Padé investigation, the "collapse", its diagnosis and the external validation; `gates/GATES_REGISTRY.md` indexes ~20 pre-registered gates |
| `out/` | All results: CSV, manifests (with SHA-256), logs, figures — including superseded versions kept on purpose |
| `docs/` | Write-up (EN/ES) and scientific documents |
| `docs/provenance/` | SHA-256 manifests: original working layout and this public release |
| `archive/` | Scripts not referenced by any document |

## Method

One person directed the work and made every decision; Claude Code (Anthropic) implemented
and ran the code; ChatGPT (OpenAI) acted as external reviewer of designs, code and data
packages. Success and failure criteria ("gates") were written before each run and never
adjusted afterwards. The AIs are tools, not authors.

## Provenance and integrity

- `docs/provenance/RELEASE_MANIFEST.txt` — SHA-256 of every file in this release.
  Check with `grep '^[0-9a-f]' docs/provenance/RELEASE_MANIFEST.txt | shasum -a 256 -c`.
- `docs/provenance/FP_REVIEW_MANIFEST.txt` and `PROVENANCE.txt` — hashes of the original
  working layout (single folder), kept as a historical record. Files that were moved keep
  their content except for a one-line import-path header, so their hashes differ.
- Manifests and logs in `out/` were produced in the original layout and may contain its
  absolute paths.

For this release the benchmark, E1, T2 and the Goldstone test were re-run in the new
layout and compared with the original outputs (see `docs/provenance/RELEASE_VERIFICATION.md`).

## Licence

Code: MIT ([`LICENSE`](LICENSE)). Documentation, figures and data: CC BY 4.0 ([`LICENSE-DOCS.md`](LICENSE-DOCS.md)).
To cite, see [`CITATION.cff`](CITATION.cff).
