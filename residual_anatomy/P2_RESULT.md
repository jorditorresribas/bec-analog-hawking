# P2 — Resultado archivado (gate de referencia de fase)

**Estado:** CERRADO con **clase C** en ω = 0.1 y ω = 0.05. Revisado por el referee externo (ChatGPT) el
2026-11-16 sobre `P2_para_referee.zip`; el referee recomienda mantener P2 cerrado en C y **no avanzar a P4a**.
Resultado *validado dentro del alcance documentado* (P2_DESIGN.md + ADDENDUM_1, malla y parámetros del manifiesto).

## 1. Ejecución
- Gate `gates/gate_P2_phase_reference.py` v10, sha256 `a0c73727d5a2d7589b948a25ca357d6777058630d49db1cf67079a839f42c51a`
  (aprobado para ejecución).
- 2026-11-16 04:11:41 → 04:35:51 CET, exit 0, pico RSS 6.27 GB, MBP 2015 (i7-4870HQ), Python 3.9.6,
  numpy 2.0.2, scipy 1.13.1 (scipy-openblas 0.3.27).
- **Integridad PASS:** I-C0, I-C1, I-C4 en las 28 ejecuciones; η_back entre 3.8e−17 y 6.2e−17.

## 2. Resultado formal (PROBADO dentro del protocolo congelado)
| | ω = 0.1 | ω = 0.05 | umbral |
|---|---|---|---|
| C2 módulos X0 ≥ 16 | pasa (8.0e−7) | pasa | 1e−6 |
| C2 fases (Q_post) X0 ≥ 16 | **3.14e−6 → FALLO** | **1.36e−6 → FALLO** | 1e−6 |
| C2 fases pre-tramo (máximo de Q_pre, estricto) | OK | OK | — |
| Estado C2 | FALLA | FALLA | |
| C3 decisivo E_X(S_box) | 1.74e−6 (falla) | 8.26e−7 (pasa) | 1e−6 |
| C3 indicativo E_28(cont)/E_28(box) | 2851 | 6424 | ≥ 100 |
| C5 | FALLA: E_ref(12,16) = 8.07e−5 | FALLA: 1.05e−4 | 5e−5 |
| **Clase** | **C** | **C** | |

C2 FALLA → clase C ("diagnosticar causa"); C3 y C5 solo se registran. **H0 (deriva de X0 = re-faseo puro
por canal) NO queda probada.** Sin reclasificación ni cambio de umbrales.

## 3. Diagnóstico post-hoc (HIPÓTESIS — no interviene en la clasificación)
Tabla `out/p2_posthoc_desviacion_vs_X28.csv` (filas `POSTHOC_HIPOTESIS`): desviación de cada fase y módulo frente a X0 = 28.
- **Corrección:** la "tasa 1.09" del manifiesto es un ajuste de los **módulos** en X0 = 8…14 respecto a
  X_ref = 16; **no** es una tasa de las fases y no sirve como evidencia de decaimiento exponencial de las fases.
- Las cuatro fases se aproximan a sus valores en X0 = 28 en ambas ω: frente a 28, ~1e−6–3e−6 en X0 = 16,
  ~5e−8 en 20, ~1e−9 en 24 (≈ tres órdenes entre 16 y 24). La envolvente decrece con tasa ≈ 0.9–1.0 por
  unidad de X0 en 16–24 (referee), **con cambios de signo** y sin caída exponencial limpia fase por fase.
- La monotonía estricta del pre-tramo se cumple para el **máximo** de las fases, **no** fase por fase
  (p. ej. |ΔI_12| sube entre X0 = 12 y 14 con ω = 0.1).
- Hipótesis compatible (no demostrada): X_ref = 16 aún no está suficientemente cerca del régimen asintótico;
  el fallo C2 apuntaría a la elección de X_ref y no necesariamente a un fallo del re-faseo.
- **ABIERTO:** si X0 = 28 representa ya el límite; si existe una cola oscilatoria concreta; el mecanismo de
  los cambios de signo. Pregunta siguiente (referee): ¿persiste la estabilización con una referencia
  independiente y un dominio extendido? — requiere diseño nuevo en papel y revisión; no se ejecuta nada.

## 4. Artefactos (sha256)
| archivo | sha256 |
|---|---|
| out/p2_manifest.json | 97d1cfcb6081176db3912b79a518c7b61c32444e123112ed87b465af4cd24a6f |
| out/p2_X0_sweep.csv | 27b213a7f664f105b08d87efa6adc301c62c1c9ebc0fcc077bd24854988c7e49 |
| out/p2_invariants.csv | bcca1bd43f451fc838e094bfc0804f383f9ca443dd1059f87d2f3ea3b933800e |
| out/p2_dx_sweep.csv | 31ca750bad08131f27185220543518afbaf7c39647a32d11fe68d19cdfcb8612 |
| out/p2_kbox_check.csv | b429bc9cee8071ba3fd84787d92a805e40192f3a545b52daff01d3d304e0a79e |
| out/p2_modal_ref.csv | 5ede578c2dacdef5a2b572ab4aad178e2c7213265c49689cb11bf79efba6a283 |
| out/p2_run.log | 2c51ef83385056ff948491cd8495391bd5e2f347a3e95ccc1cac2c93290e9fb1 |
| out/p2_exit_code.txt | 9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa |
| out/p2_posthoc_desviacion_vs_X28.csv (POST-HOC) | 479255529d9fe6d286b7a5de0c8acccd5af079d13f44debc2eabce2a5197b2ed |
| P2_DESIGN.md | 06b504f9d7d4628824b5c7b2c309b9563c2e6bf16f413ed307e570dbf40064ab |
| P2_DESIGN_ADDENDUM_1.md | 5d663ab1256547856f6d9e3013ec349ea8f1c1b4155da14cc36f20b2a0ed556a |

Los hashes de los tres CSV principales coinciden con los registrados en el manifiesto (verificado por el referee).

## 5. Consecuencia para el roadmap
No se avanza a P4a. Se pausa P4a–P6 para la auditoría de novedad (Macher & Parentani 2009, Finazzi &
Parentani 2012, Coutant–Parentani–Finazzi) y la propuesta de un gate externo (reproducir una figura publicada
con el solver v1.0). `BLOCK2_LOG.md` no se ha modificado; registrar allí solo con autorización.

*Archivado 2026-11-16.*
