# Simulating Analogue Hawking Radiation with Cross-Review Between Two AIs

**An acoustic black hole, a spurious result, and how its origin was identified**

**Author:** Jordi Torres Ribas · **Contact:** jorditorresribas+bec@gmail.com
**Version:** 1.5 (English translation of the Spanish version) · **Date:** 18 November 2026

---

The goal was to simulate on a computer a laboratory analogue of a black hole: a gas of ultracold atoms flowing faster than sound which, according to theory, should emit faint thermal radiation, the acoustic counterpart of Hawking radiation. Halfway through, the simulation produced a striking result that could have been presented as a discovery. It wasn't one. This document includes the mistakes made along the way, because they are the most informative part.

## Summary

- I used AI to build a simulator of a laboratory analogue of a black hole, and subjected it to validation rules fixed before every test.
- Halfway through, a striking result appeared: the radiation "collapsed" when one model parameter was changed. **It wasn't physics: it was a flaw in the model.** The symmetry test identified the defect, with agreement to four significant figures between the measured and predicted residuals. Once it was fixed, the "collapse" disappeared.
- The instrument reproduces an **exact** result published by another research group: first with a relative error of 3×10⁻⁵ in the leading coefficient and then, in a second test whose criteria were published before it was run, with an error of 10⁻⁸.
- Once the error was fixed, it turned out that the parameter being varied did not do what it was meant to do: it is exactly equivalent to changing the width of the horizon.
- The audits detected real errors in the checks themselves and in the documentation. They are documented below, together with who detected each one.

## What this project demonstrates

1. **Directing AI on a hard technical problem without errors reaching the conclusions**: one AI builds, another audits, and a person decides.
2. **Turning "it seems to work" into "it has been checked"**: success and failure criteria written down *before* running anything. There are two exceptions, both test tolerances, declared in the code itself with their justification: they were relaxed after a failure, once the numerical noise that explained it had been measured (`tests/validate_lambda.py` and `residual_anatomy/diag_scaling_equivalence.py`).
3. **Spotting a misleading result and identifying where it came from**, instead of publishing it or burying it.
4. **Full traceability**: every result stored with its data, a manifest and cryptographic fingerprints (SHA-256), so that a third party can reproduce it.
5. **Honesty about scope**: keeping what is PROVEN separate from what is a HYPOTHESIS and what is OPEN, and keeping negative results too.

---

## 1. The problem, in plain terms

A black hole traps light because nothing can travel faster than light. Something similar can be built with sound: if a gas flows faster than the speed of sound inside it, sound waves trying to travel upstream are swept back. There is a sonic "horizon".

In 1974, Stephen Hawking predicted that black holes emit faint thermal radiation. Their acoustic analogue should do the same: pairs of quantum vibrations are created, one escapes and the other falls in. In a gas of atoms cooled almost to absolute zero (a *Bose-Einstein condensate*), this effect can be calculated and, in principle, measured.

This project computes that radiation with a numerical solver written from scratch. The key quantity is the **"temperature"** of the emitted radiation, and the question was how it changes when the physics at very short length scales is modified. It is a laboratory version of an open problem in theoretical physics: the "trans-Planckian problem".

## 2. How the work was organised

**Three separate roles:**

- **A person (me)** decides what gets done, authorises every run and settles disagreements.
- **Claude Code (Anthropic)** implements: writes the code, runs it, checks the numbers and drafts the reports.
- **ChatGPT (OpenAI)** acts as **external reviewer**: reviews designs before any code exists, audits the code before it is run and, at the end, audits the data packages.

The work was directed by someone with no academic training in physics. For that reason, the method was designed so that conclusions would not rest on anyone's authority, human or artificial, but on explicit, documented and reproducible checks.

**Rules applied without exception:**

- **Design on paper → review → implementation → audit → run.** No significant computation was launched without an approved design.
- **Criteria fixed in advance** ("gates"). Before running, it is written down which result counts as success, which as failure and which as inconclusive. If the outcome is unfavourable, it is reported as is; the threshold does not move. (There were two declared exceptions in the tests: see "What this project demonstrates", item 2.)
- **Mandatory labels**: PROVEN / HYPOTHESIS / OPEN, and always "validated within the documented scope".
- **Persistence**: every result comes with its CSV files, a manifest recording the software environment, and SHA-256 fingerprints of code and data.
- **Real independence**: two models from the same family do not count as two independent checks, because they can share biases. That is why the reviewer comes from a different company.

In total, the project accumulated about twenty documented gates.

## 3. The result that looked like a discovery

The basic simulator (version 1.0) was validated first and then frozen. It reproduces the theoretical prediction by Finazzi and Parentani (2011): the temperature obtained is between 95.7 % and 98.5 % of the ideal Hawking temperature, approaching 100 % as the horizon becomes smoother. (With the corrected temperature those authors propose, the agreement is better than 0.75 %, but that figure depends on a parameter they fitted to their own profiles.)

A parameter, λ, was then added with the aim of changing only the very-short-scale physics and nothing else. The result was surprising: moving λ away from its original value, in either direction, **made the temperature of the radiation drop to zero**. The expected effect was a few per cent; this was a total "collapse".

A result like that could have been presented as a finding. Instead it was investigated (section 4) and turned out to be a modelling error. The investigation also left a second lesson, less visible but just as important: **once the error was fixed, λ no longer did what it was meant to do.** Changing λ was shown to be exactly equivalent to making the horizon wider or narrower: the results with λ match those of the original simulator on a rescaled horizon with a relative difference of a few parts in a billion (3.5×10⁻⁹). In other words, λ does not isolate short-scale physics; it is another way of moving a parameter that was already known. The original question, as posed, did not have the content it seemed to have.

## 4. How the error was identified

1. **Pause and literature audit.** Before building a more complex solver to study the "collapse" (several more steps of work), we checked whether the question had already been answered. Reading the original papers (Macher and Parentani 2009; Finazzi and Parentani 2012) showed that **very similar calculations, with other flow profiles, had already been done and showed no collapse at all**; our result was in strong tension with them.
2. **A concrete hypothesis.** When λ was changed, the code modified the equation for the vibrations but did not update one term describing the background state of the gas. The simulated gas was no longer in equilibrium.
3. **Symmetry test.** If the background is in equilibrium, there is one specific vibration, the *Goldstone mode*, that must have exactly zero frequency. We measured it:
   - the original simulator passed the test (the residual decreased with grid refinement, as expected from discretisation error);
   - the λ version violated it, and the residual matched the formula derived on paper with a ratio of **1.0000**;
   - correcting the term brought the zero mode back.
4. **The decisive test.** The original calculation was repeated with the corrected background, under the same conditions. **The "collapse" disappeared.** The temperature now varied smoothly, from 0.98 to 0.89 times the reference value, and it matches exactly the original simulator applied to a horizon of rescaled width (section 3).
5. **Additional consequence.** The more complex solver that had been planned (several steps of work) **had the same design flaw**. It was caught before it was built.

## 5. External validation

An instrument agreeing with itself is not enough: it could carry a systematic error. So we looked for an **exact** result from another group, for the same kind of system: the *waterfall* configuration of Larré, Recati, Carusotto and Pavloff (2012): a flow falling down a potential step, whose exact solution is known, with closed-form formulas for the radiation.

With criteria fixed before running, the simulator reproduced:

- the leading coefficient of the radiation with a relative error of **3×10⁻⁵**;
- the subleading coefficient with an error **below 10⁻³**;

in both cases tested, with numerical convergence at the expected order and the sensitivity to the fitting method quantified.

A later audit showed that the 3×10⁻⁵ came almost entirely from how the data were extrapolated to zero frequency, not from the simulator. To check this without adjusting anything after the fact, the criteria for a second test (E1-bis) were written down and stored in a commit **before** it was run; that commit is public (`001a7ac`). In two cases that had not been looked at, the simulator reproduced the leading coefficient with errors of **5×10⁻⁹ and 4×10⁻⁸**, and the subleading one with **8×10⁻⁷ and 3×10⁻⁶**.

The equations taken from the papers were later checked against their literal text. The check confirmed the conventions of the benchmark (Finazzi and Parentani's scale Λ) and found one citation with the wrong equation number, now corrected.

## 6. What went wrong along the way

It is the most revealing part of the method. Seven real errors, chosen from the nineteen on record (the full list is in the appendix):

| What went wrong | Who detected it | What it taught / what changed in the process |
|---|---|---|
| The theoretical reference temperature carried a spurious factor of ½; it looked like a 120 % physical discrepancy | The reviewer insisted on auditing conventions; the implementer found the error by checking against the literature | Before reading a discrepancy as physics, normalisations and conventions are audited |
| A false alarm ("the conditioning blows up") came from a bug in a numerical library (OpenBLAS) when solving the adjoint system | Three-phase adversarial audit (implementer, second model, external reviewer) | One anomalous number is not enough to conclude anything: it is triangulated through an independent route (here, three certified methods) |
| A `max()` call was hiding NaN values in the code of a gate (one of the reasons that gate needed ten versions) | The reviewer, auditing the code before it was run | "Fail-safe" code: mandatory coverage and explicit handling of NaN |
| The "collapse" of the radiation (section 3) | Literature audit requested by the person; the implementer identified the defect | Any modified model must pass a symmetry test before its results are interpreted |
| The first measurement in that test was contaminated by the edges of the domain | The implementer, noticing that the error grew with resolution when it should have shrunk | Measure where it makes sense, and keep the wrong version on record |
| A convergence check returned zero by construction: it checked nothing | The reviewer, auditing the data package | The checks get audited too: convergence orders computed separately and extrapolation declared |
| A channel check did not verify what it claimed to verify | The reviewer, auditing the package | Explicit check, logged frequency by frequency |

None of these errors made it into a final conclusion. The last two affected the checks, not the physics: after fixing them, the recomputed data were byte-for-byte identical.

## 7. What is not being claimed

- **There is no new physics.** The robustness of the radiation against changes at short scales was already established (Macher, Finazzi, Parentani and others). This project reproduces it; it does not discover it.
- That the "collapse" was an artefact is shown **within the scope tested** (one flow profile, one range of λ). With the corrected background the results match the rescaled original simulator up to rounding, so there is no room for another mechanism within the tested scope; convergence with respect to the domain size has still not been studied.
- The parameter λ does not probe short-scale physics at fixed geometry: with the corrected background it is a rescaling of the horizon width. The conclusions about correlations obtained with λ have been withdrawn.
- The external validation covers the basic simulator (λ = 1), in the configuration and the two cases tested.
- Some intermediate results ended up inconclusive or negative, and are recorded as such. For example, a phase-reference gate ended in "class C" (criterion not met) and was archived as it was.

## 8. In one sentence

A person organised two AIs into separate roles with rules fixed in advance. The system produced a misleading result, detected it, identified where it came from, fixed it and validated the instrument against an exact published result, with every step reproducible.

---

*Tools: implementation with Claude Code (Anthropic); external review with ChatGPT (OpenAI). Custom code in Python (NumPy, SciPy).*

---

## Technical appendix (for physicists)

**System.** 1D Bose-Einstein condensate, stationary Bogoliubov-de Gennes (BdG) equations; ℏ = m = g = 1. Density profile n(x) of tanh form (Mach 0.5 upstream, 2 downstream, width w = 2); not to be confused with the Larré et al. *waterfall* configuration used in E1. 3×3 S-matrix by exact modal matching with a second-order finite-difference operator. Observable: n_L(ω), the occupation of the upstream Hawking channel (Finazzi and Parentani 2011, Eq. 42). T_eff = lim ω·n_L.

**v1.0 benchmark (frozen).** T_up/(κ_local/2π) = 0.957 → 0.985 for w ∈ {2, 3, 4, 6, 8}, against Finazzi and Parentani, PRD 83, 084010 (2011). With their broadened-horizon κ̄, T_up/(κ̄/2π) ∈ [0.9926, 1.0050]; this second figure is conditional on their window width d = 1.58 (with d = 1 or 2 the maximum deviation is 2.4 % and 3.4 %). The scale Λ = 2mc²/ℏ of their Eq. 9 has been checked against the paper.

**The artefact.** The operator with dispersion (ω−vk)² = c²k² + λk⁴/4 is equivalent to BdG with mass m′ = 1/√λ and coupling g′n = c²/√λ, but it kept the background quantum pressure Q_p fixed instead of using Q_p,λ = √λ·Q_p. Goldstone residual on Φ_G = (√n, −√n): ½(1 − √λ)(√n)″. Measured / predicted = 1.0000, independent of the grid. With Q_p,λ the residual drops to O(dx²). The planned Padé operator had the residual (K_P − K_1)√n (measured ratio ≈ 1.00).

**T2: n_L with the corrected background** (w = 2, ω/(κ/2π) ∈ [0.05, 0.27], N = 9600 and 14400):

| λ | 0.25 | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|---|
| T_eff/T_H, inconsistent background | −0.006 | −0.027 | 0.957 | −0.053 | −0.017 |
| T_eff/T_H, consistent background | 0.979 | 0.971 | 0.957 | 0.932 | 0.892 |

(T_H = κ_local/2π. With the inconsistent background, n_L ∝ ω at low frequency, so the extrapolated T_eff comes out ≈ 0, slightly negative because of the extrapolation. With the consistent background, the difference between grids is ≤ 6×10⁻⁵.)

**λ ↔ width equivalence.** With the consistent background, x′ = x/√λ, ω′ = √λ·ω map the λ operator on a profile of width w onto the physical operator (λ = 1) on a profile of width w/√λ: M₁(x′) = √λ·M_λ(x), also on the grid (L′ = L/√λ, same N). Checked at the 90 points of T2: relative difference in n_L ≤ 3.2×10⁻⁶ with the standard solve (its rounding noise) and ≤ 3.5×10⁻⁹ with equilibration and iterative refinement. The "consistent background" row is therefore the curve T/T_H against the width w/√λ = 4 → 1.

**E1: validation against Larré et al., PRA 85, 013621 (2012).** a = |f_{u,d2}|² (their Eq. A8) and b = 2Re(f*h) (Eq. A14), extracted from ω·n_L = a + bω + cω² + dω³:

| m_u | a numerical / exact | error | b numerical / exact | error |
|---|---|---|---|---|
| 0.5 | 0.131737 / 0.131734 | 2.4×10⁻⁵ | −0.490187 / −0.489796 | 8.0×10⁻⁴ |
| 0.7 | 0.066910 / 0.066908 | 2.9×10⁻⁵ | −0.498434 / −0.497974 | 9.2×10⁻⁴ |

Γ and T_H (their Eq. 57) reproduced. Integrity: regression of the new driver against v1.0 = 0; Gross-Pitaevskii residual converging at second order; pseudo-unitarity defect with observed order between 1.87 and 1.98; channels verified at 240/240 frequencies.

**E1-bis (criteria pre-registered in commit `001a7ac`, before the run).** The E1 error turned out to be the bias of the ω → 0 extrapolation (window W1 with a cubic fit), not rounding (10⁻⁹–10⁻⁷ on that grid) nor discretisation. Window W2 = [0.05, 0.3]·T_H with a fourth-order fit, dx = 0.01/0.005/0.0025, refined solve and Richardson with p = 2:

| m_u | error in a (A8) | error in b (A14) | observed order |
|---|---|---|---|
| 0.6 | 5.4×10⁻⁹ | 7.9×10⁻⁷ | 2.04 |
| 0.7 | 4.2×10⁻⁸ | 3.1×10⁻⁶ | 2.00 |

Main criterion: error in a ≤ 10⁻⁶ in both, and window check (W2 against W3) ≤ 10⁻⁶: PASS. m_u = 0.5 was repeated as a replica, not counting towards the verdict (1.1×10⁻⁸).

**Note on the two reference temperatures.** In the v1.0 benchmark the ratio is taken against T_κ̄ = κ̄/(2π), the temperature built from the averaged gradient of Finazzi and Parentani (2011); that is why it comes out ≈ 1.00. In T2 the ratio is taken against T_local = κ_local/(2π), built from the gradient at the horizon; that is why λ = 1 gives 0.957. These results are not contradictory: they are ratios with different denominators.

**Full error log (14).** The last four affected the checks and the documentation, not the physics: after fixing them, the recomputed data were byte-for-byte identical.

| What went wrong | Who detected it | What changed in the process |
|---|---|---|
| The theoretical reference temperature carried a spurious factor of ½; it looked like a 120 % physical discrepancy | The reviewer insisted on auditing conventions; the implementer found the error by checking against the literature | Rule: always audit normalisations and conventions before interpreting a discrepancy |
| A fit with a fixed frequency window failed when the scale of the problem changed | The fit residuals (a safeguard requested by the reviewer) | Fit windows scaled with the physics of each case |
| A flux bookkeeping step used the wrong sign (2 % deficit) | The reviewer pointed to the step to check; the implementer fixed it | Always distinguish quantities that look interchangeable but are not |
| Two false alarms: a "blow-up" of the conditioning and an apparent non-convergence | Three-phase adversarial audit (implementer, second model, external reviewer) | Triangulate through an independent route before concluding anything from an anomalous number |
| The underlying defect in one of them: a bug in a numerical library (OpenBLAS) when solving the adjoint system | The same audit | Alternative numerical route certified with three independent methods |
| The "stabilised" solver that had been validated was for a simplified version of the model, not the one that was needed | Adversarial audit | Validate the instrument that will actually be used, not the easiest one |
| Results that were only printed to screen and then lost | Adversarial audit | Every result is stored with CSV, manifest and SHA-256 |
| Robustness bugs in the code of one gate (ten versions before approval), for example a `max()` that hid NaN values | The reviewer, auditing the code before it was run | "Fail-safe" code: mandatory coverage and explicit handling of NaN |
| The "collapse" of the radiation (section 3) | Literature audit requested by the person; the implementer identified the defect | Mandatory symmetry test for any modified model |
| The first measurement in the symmetry test was contaminated by the edges of the domain | The implementer, noticing that the error grew with resolution when it should have shrunk | Measure in the interior and keep the wrong version on record |
| A convergence check returned zero by construction: it checked nothing | The reviewer, auditing the data package | Convergence orders computed separately and extrapolation declared |
| A channel check did not verify what it claimed to verify | The reviewer, auditing the package | Explicit check, logged frequency by frequency |
| The delivered package could not be reproduced on its own | The reviewer | Self-contained packages, with dependencies and checksums |
| A cryptographic fingerprint attributed to the wrong file name | The reviewer | Keep the executed version and the updated version of each document separate |

**Errors found in the November 2026 audit (5).** The second one withdrew conclusions from documents in the repository (not from this text); the others do not change any published number.

| What went wrong | Who detected it | What changed in the process |
|---|---|---|
| The tests printed results but never failed: they had no assertions | The implementer, in the read-only audit requested by the reviewer | Every test checks against a tolerance and returns an error; `tests/run_all.py` runs them all |
| The correlations computed with λ inherited the same background defect as the "collapse" | The implementer, in that same audit, while reviewing which results used the defective operator | Those conclusions were withdrawn, with a public notice |
| The absolute amplitude of the density-density correlation came out at half its value: the paper's prefactor belongs to a different mode normalisation. The equation was also cited with the wrong number | The implementer flagged the possible factor of 4π in the phase 1 audit; the reviewer located the formula (Eq. 25) and the wrong equation number; the implementer settled it with an exact test in a homogeneous condensate | Absolute-normalisation test with a tolerance; equations checked against the literal text of the papers |
| A published figure had been generated with an earlier version of the code than the one described in the text (without the smoothing window) | The implementer, while regenerating it | Figures are regenerated with the repository scripts before publishing |
| An audit note attributed to the E1 calculation a rounding noise measured in a different configuration | The implementer, measuring it in the right configuration | Measure each quantity in the configuration it refers to |

**References.**
- J. Macher, R. Parentani, PRA 80, 043601 (2009)
- A. Recati, N. Pavloff, I. Carusotto, PRA 80, 043603 (2009)
- S. Finazzi, R. Parentani, PRD 83, 084010 (2011)
- S. Finazzi, R. Parentani, PRD 85, 124027 (2012)
- A. Coutant, R. Parentani, S. Finazzi, PRD 85, 024021 (2012)
- P.-É. Larré, A. Recati, I. Carusotto, N. Pavloff, PRA 85, 013621 (2012)

**Where everything is** (in the project repository, <https://github.com/jorditorresribas/bec-analog-hawking>; the detailed project documents are in Spanish):
- `docs/FP_BENCHMARK.md`: v1.0 benchmark.
- `residual_anatomy/NOVELTY_AUDIT.md`: novelty audit, including Appendix A (derivation, Goldstone test, T2).
- `residual_anatomy/E1_RESULT.md` and `residual_anatomy/GATE_EXT_LRCP_DESIGN.md`: external validation; `residual_anatomy/E1BIS_PREREG.md`: E1-bis criteria.
- `docs/CORRELATIONS.md` §9: λ ↔ width equivalence.
- `tests/run_all.py`: all validations.
- `residual_anatomy/P2_RESULT.md`: an archived negative result.
- `out/`: CSV files and manifests.
