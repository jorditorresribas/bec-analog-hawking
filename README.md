# Simulating analogue Hawking radiation with cross-review between two AIs

**An acoustic black hole, a spurious result, and how its origin was identified**

Jordi Torres Ribas · v1.4 · 17 November 2026

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
| Benchmark v1.0 vs Finazzi & Parentani, PRD 83, 084010 (2011) | T_up/(κ̄/2π) ∈ [0.9926, 1.0050] for w ∈ {2,3,4,6,8} (≤ 0.75 %) | `python3 benchmarks/fp_v1/FP_check_reproduction.py` (~30 s) |
| External validation vs Larré, Recati, Carusotto & Pavloff, PRA 85, 013621 (2012) | leading coefficient: relative error 2.4–2.9·10⁻⁵; secondary: < 10⁻³ | `python3 residual_anatomy/gates/gate_E1_lrcp.py` (~1 min) |
| The "collapse" was an artefact (inconsistent background) | with a consistent background T_eff/T_H varies smoothly 0.979 → 0.892 for λ = 0.25 → 4 | `python3 residual_anatomy/diag_T2_qp_consistent.py` (~1 min) |
| Goldstone (symmetry) test identifying the defect | measured / predicted residual = 1.0000 | `python3 residual_anatomy/diag_goldstone.py` (~7 min) |

**Scope.** No new physics is claimed: the robustness of the radiation against short-scale
modifications was already established (Macher, Finazzi, Parentani and others). The
artefact is shown to be one within the tested scope (one flow profile, one range of λ).
The external validation covers the basic solver (λ = 1). See §7 of the write-up.

## Quick start

```bash
# Python 3.9.6 (see .python-version); pinned dependencies:
python3 -m pip install -r requirements.txt

# Always run from the repository root:
python3 benchmarks/fp_v1/FP_check_reproduction.py
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

Code: MIT. Documentation, figures and data: CC BY 4.0. See [`LICENSE`](LICENSE).
To cite, see [`CITATION.cff`](CITATION.cff).
