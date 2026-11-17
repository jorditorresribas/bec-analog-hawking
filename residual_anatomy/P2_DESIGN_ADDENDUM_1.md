# P2 — ADDENDUM 1 al diseño congelado `P2_DESIGN.md`

**Estado:** APROBADO por el referee externo (ChatGPT) el 2026-11-16, con la precisión sobre la prioridad de
FALLO incorporada en A2. `P2_DESIGN.md` **sigue congelado e intacto**; este addendum lo complementa.

**Motivo:** durante la auditoría del código de `gates/gate_P2_phase_reference.py` (v3) se detectó que el
§3.2(i) solo contempla el caso "Sin I_ij definidas → C2 solo con módulos, y se hace constar", y **no define**
qué hacer cuando las fases plaqueta son evaluables en un tramo y no en el otro. El addendum resuelve esa
ambigüedad **sin cambiar ningún umbral, parámetro ni otra prueba**.

---

## A1. Conjuntos fijos de plaquetas por tramo
- **Q_post** = plaquetas I_ij definidas en X_ref y en **todos** los X0 ∈ {16, 20, 24, 28}.
- **Q_pre** = plaquetas I_ij definidas en X_ref y en **todos** los X0 ∈ {8, 10, 12, 14}.
- E_X^I de cada tramo se calcula **solo** sobre su conjunto, para comparar siempre las mismas cantidades
  (una plaqueta definida solo en algunos X0 de un tramo no entra en la decisión de ese tramo).

## A2. Estado de las fases (mutuamente excluyentes, evaluados en este orden)
1. **NINGUNA:** no hay **ninguna** I_ij definida en **ningún** X0 del barrido (incluido X_ref). Es
   exactamente el caso "Sin I_ij definidas" del §3.2(i).
2. **FALLO (prioritario):** algún tramo evaluable (Q no vacío) incumple su criterio (E_X^I ≤ 1e−6 en X0 ≥ 16;
   decrecimiento **estrictamente** monótono en X0 = 8…14). FALLO prevalece **aunque el otro tramo esté vacío**.
3. **OK:** Q_post y Q_pre no vacíos y ambos tramos cumplen su criterio.
4. **PARCIAL:** cualquier otro caso: existe alguna plaqueta definida en algún X0, pero no se pueden evaluar
   ambos tramos completos (algún Q vacío) y ningún tramo evaluable falla.

## A3. Estado de C2
| Módulos (C2) | Fases | Estado C2 |
|---|---|---|
| incumplen | cualquiera | **FALLA** |
| cumplen | FALLO | **FALLA** |
| cumplen | OK | **PASA** |
| cumplen | NINGUNA | **PASA_SOLO_MODULOS** (anotación obligatoria: "las fases no aportan evidencia favorable") |
| cumplen | PARCIAL | **INCONCLUSO** |

## A4. Clasificación (amplía §4; solo con integridad PASS)
| Estado de fases | Estado C2 | Consecuencia |
|---|---|---|
| OK | PASA | Clasificación ordinaria (B / A / D / D* según C3, C5) |
| NINGUNA | PASA_SOLO_MODULOS | Clasificación ordinaria, **con anotación** |
| FALLO (o módulos incumplen) | FALLA | **Clase C** |
| PARCIAL | INCONCLUSO | **Clase I** — "inconcluso por fases parcialmente evaluables" |

- La **clase I prevalece sobre C3 y C5**: ambos se calculan y persisten, pero **no** modifican la clase.
- La clase I **nunca** puede ser A.

## A5. Alcance
Sin cambios en umbrales, parámetros, mallas ni en ninguna otra prueba. Solo se precisan casos que el texto
congelado dejaba indefinidos. El manifiesto del gate registra el SHA-256 de este addendum junto al de
`P2_DESIGN.md`.

*Addendum creado 2026-11-16 tras aprobación del referee. No autoriza la ejecución del gate.*
