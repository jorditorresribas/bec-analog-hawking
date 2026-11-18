# E1 — Gate externo Larré et al. 2012 (*waterfall*): resultado

**v2 (2026-11-16, tras la auditoría del referee del paquete):** tres fallos de los controles, sin afectar a la física.
1. I3: la "extrapolación" de la v1 daba 0 por construcción (p sacado de los mismos dos puntos). Ahora se informan
   los órdenes grueso y fino de W1 y W2, y la extrapolación solo con p = 2 declarado.
2. I4: la v1 contaba modos salientes sin exigir que fueran propagantes. Ahora el guard exige Lout = 1 y Rout = 2, y
   los recuentos por frecuencia van a `e1_lrcp_integrity.csv` (240 filas, todas correctas).
3. El paquete no era autocontenido; ahora incluye las dependencias de v1.0 y `out/waterfall_state.npz`.
Reejecutado: `e1_lrcp_fits.csv` y `e1_lrcp_nL.csv` son **idénticos byte a byte** a los de la v1 (a y b no cambian).
Los resultados de la v1 se conservan como `out/v1_e1_lrcp_*`.

**Clasificación: PASS** (integridad PASS; a y b PASAN en m_u = 0.5 y m_u = 0.7).
*Validado dentro del alcance documentado:* λ = 1 (operador v1.0), configuración *waterfall* de LRCP, m_u ∈ {0.5, 0.7},
coeficientes de baja frecuencia de n_L = |S_{u,d2}|², ventanas ω/T_H ∈ [0.05, 0.6].
Ejecutado el 2026-11-16 (48 s). Diseño: `GATE_EXT_LRCP_DESIGN.md` (I3 del referee + §7, fijado antes de implementar).
El referee **no** ha revisado el §7 ni el código antes de ejecutar (decisión del usuario: una sola sesión); se le
entrega después.

**Estado tras la revisión provisional del referee (2026-11-16):** "PASS comunicado, auditoría externa pendiente".
El referee revisó la transcripción, **no** el paquete. Formulación que se adopta hasta su auditoría:
> El solver v1.0 reproduce los coeficientes analíticos de baja frecuencia de Larré et al. para la configuración
> *waterfall* y los dos números de Mach ensayados. Los controles muestran convergencia espacial de segundo orden y
> convergencia de la pseudo-unitariedad. Pendiente de revisión independiente del diseño detallado (§7), la
> implementación y los datos.

**Sensibilidad a la ventana y al orden del ajuste** (pedida por el referee; malla fina dx = 0.005, L = 30;
datos de `out/e1_lrcp_fits.csv`, error relativo frente a la teoría):

| m_u | W1, orden 2 | W1, orden 3 (referencia) | W2, orden 2 | W2, orden 3 |
|---|---|---|---|---|
| 0.5 — a | −1.3e−4 | −1.1e−7 | −3.1e−5 | −2.3e−5 |
| 0.5 — b | +2.3e−3 | −7.8e−4 | +3.3e−4 | −8.0e−5 |
| 0.7 — a | −1.0e−4 | +3.6e−5 | +4.3e−7 | +9.7e−6 |
| 0.7 — b | +2.5e−3 | −9.3e−4 | +3.3e−4 | −1.2e−4 |

En las cuatro combinaciones: a dentro de 1.3e−4 y b dentro de 2.5e−3 (piso de b: 5e−3). La combinación más alejada
es W1 de orden 2 (ventana más alta con ajuste de menor orden: truncación esperable del desarrollo en ω).

## Referencia
P.-É. Larré, A. Recati, I. Carusotto, N. Pavloff, PRA 85, 013621 (2012), arXiv:1110.4464 — Ecs. 17, 44, 56, 57, A8, A14.
Resultado **analítico exacto** del modelo BdG del *waterfall* (grupo independiente).

## Integridad
| Prueba | m_u = 0.5 | m_u = 0.7 |
|---|---|---|
| I1 regresión del driver frente a v1.0 `solve_inject` (ω = 0.1) | **0.0** relativo (idéntico) | — |
| I2 residuo GP del fondo (máx. interior, dx = 0.02/0.01/0.005) | 2.4e−4 → 6.0e−5 → 1.5e−5 (orden 1.99, 2.00) | 3.5e−5 → 8.8e−6 → 2.2e−6 (orden 2.00) |
| I3 defecto de pseudo-unitariedad, W1 (máx.) | 6.4e−3 → 1.7e−3 → 4.5e−4; p grueso/fino 1.87/1.96; extrap. (p = 2 declarado) 1.7e−5 | 4.9e−3 → 1.3e−3 → 3.3e−4; p 1.91/1.98; extrap. 7.1e−6 |
| I3 defecto de pseudo-unitariedad, W2 (máx.) | 1.6e−2 → 4.3e−3 → 1.1e−3; p 1.93/1.98; extrap. 2.0e−5 | 1.4e−2 → 3.7e−3 → 9.3e−4; p 1.95/1.97; extrap. 2.3e−5 |
| I4 canales (Lin = 1, Rin = 2 con 1 de norma −, salientes **propagantes** Lout = 1 y Rout = 2) | 120/120 ω·malla | 120/120 ω·malla |

La discontinuidad del potencial **no** degradó el orden: es 2 en todas las pruebas.

## Ciencia
| | m_u = 0.5 | m_u = 0.7 |
|---|---|---|
| a = lím ω·n_L, numérico / teoría (A8) | 0.131737 / 0.131734 | 0.066910 / 0.066908 |
| error relativo de a / presupuesto δ | **2.4e−5** / 1.3e−4 | **2.9e−5** / 1.4e−4 |
| b, numérico / teoría (A14) | −0.490187 / −0.489796 | −0.498434 / −0.497974 |
| error relativo de b / presupuesto δ | **8.0e−4** / 3.1e−3 | **9.2e−4** / 3.4e−3 |
| Γ = −2b, numérico / teoría | 0.98037 / 0.97959 | 0.99687 / 0.99595 |
| T_H = a/Γ, numérico / Ec. 57 | 0.13437 / 0.13448 | 0.06712 / 0.06718 |

- a coincide con la fórmula exacta en **3 × 10⁻⁵ relativo**, por debajo de su propio presupuesto numérico.
- b (subdominante) coincide en **< 10⁻³**. Su presupuesto está dominado por la elección de ventana y orden del ajuste
  (δ_est ≈ 1.5e−3), no por la malla ni el dominio.
- Convergencia: Richardson en dx con orden observado 2.0 (a) y 1.8–2.1 (b); dependencia del dominio ≤ 1.4e−6.

## Qué valida y qué no
- **Valida:** el instrumento de v1.0 (operador, *matching* modal, normalización por corriente, identificación de
  canales y extracción de n_L) reproduce un resultado exacto publicado por un grupo independiente, para un perfil
  distinto al nuestro (discontinuo y analítico) y en dos números de Mach.
- **No valida:** el operador con λ ≠ 1 (T2 solo queda cubierto de forma indirecta, por compartir normalización y
  extracción); la curva completa n_L(ω) hasta Ω (LRCP solo da los coeficientes de baja frecuencia); el régimen
  m_u = 0.3 (no ejecutado: ampliación opcional del diseño).

## Artefactos
`out/e1_lrcp_nL.csv`, `out/e1_lrcp_fits.csv`, `out/e1_lrcp_integrity.csv`, `out/e1_lrcp_manifest.json` (sha256 de
los tres CSV, del driver, de `fd_bdg_operator.py`, de `mode_matching.py` y del diseño; verificados tras la ejecución),
`out/e1_lrcp_run.log`. Un primer intento abortó **al importar**, sin calcular nada ni escribir resultados (había que
ejecutarlo desde la raíz del proyecto por una ruta relativa de v1.0); su log está en `out/e1_lrcp_run_intento1_import.log`.

## Nota de auditoría (fase 1, 2026-11-17, corregida la misma tarde): de dónde sale el error en a

**Corrección.** Una primera versión de esta nota decía que el error en a (2,4·10⁻⁵) estaba solo 3–4 veces por
encima de un «suelo de redondeo» de ~7·10⁻⁶. **Era incorrecto:** ese 7·10⁻⁶ se midió en el perfil tanh (w = 1,
N = 14400; `diag_scaling_equivalence.py`), no en la malla de E1. Medido en la configuración de E1:

- **Redondeo en la malla de E1** (`residual_anatomy/diag_conditioning.py`, `out/diag_conditioning.csv`,
  m_u = 0,5, L = 30, ventana W1): error del solve `spsolve` frente a una solución con refinamiento iterativo
  (residuo en long double) de **1,1·10⁻⁹ (dx = 0,02), 1,4·10⁻⁸ (0,01) y 1,1·10⁻⁷ (0,005)**. Mover L en 1 ulp
  cambia n_L ~10⁻¹² con la solución refinada: el problema está bien condicionado y el ruido es del LU sin
  equilibrar. En a, el efecto es ~10⁻⁷; el Richardson de E1 sale igual con cualquier solver
  (2,41·10⁻⁵ con `spsolve`, 2,39·10⁻⁵ con equilibrado + refinamiento).
- **El error de 2,4·10⁻⁵ es casi entero sesgo de la extrapolación ω → 0** (ventana W1 con ajuste cúbico,
  fijados antes de ejecutar E1), no discretización ni redondeo (`residual_anatomy/diag_convergence.py`,
  `out/diag_convergence_e1.csv`). La discretización es de orden 2 limpio y el Richardson la elimina.
- Con ventanas más bajas u orden 4, a coincide con A8 a 1–2·10⁻⁸ (m_u = 0,5). **Esto se exploró después de
  conocer el resultado de E1**: es un diagnóstico, no forma parte del veredicto de E1. Para afirmarlo hace
  falta el gate E1-bis, con criterios fijados y en commit antes de ejecutar (`residual_anatomy/E1BIS_PREREG.md`).
- El error en b (8·10⁻⁴) no se ha reanalizado; está muy por encima del redondeo.
- El veredicto de E1 no cambia.

## Gate E1-bis (pre-registrado, 2026-11-17): **PASS**
Criterios en `residual_anatomy/E1BIS_PREREG.md`, en commit (`001a7ac`) **antes** de ejecutar; el script lo comprobó (I0) y
guardó el hash en `out/e1bis_manifest.json`. Ventana W2 con ajuste de orden 4, dx = 0,01/0,005/0,0025, solve refinado,
Richardson con p = 2 (`residual_anatomy/gates/gate_E1bis_lrcp.py`, 110 s).

| m_u | error relativo en a (A8) | \|W2 − W3\| | error relativo en b (A14) | p observado |
|---|---|---|---|---|
| 0,6 (confirmatoria) | 5,4·10⁻⁹ | 8,5·10⁻⁹ | 7,9·10⁻⁷ | 2,04 |
| 0,7 (confirmatoria) | 4,2·10⁻⁸ | 3,7·10⁻⁸ | 3,1·10⁻⁶ | 2,00 |
| 0,5 (réplica, no cuenta) | 1,1·10⁻⁸ | 2,7·10⁻⁹ | 8,9·10⁻⁷ | 2,01 |

C1 (≤ 10⁻⁶), C2 (control de ventana), C3 (nivel fuerte ≤ 10⁻⁷, solo informativo) y C4 (b ≤ 10⁻³) se cumplen; integridad
I0–I3 correcta. Con la extrapolación en ω controlada, A8 y A14 quedan confirmadas a 10⁻⁸ y 10⁻⁶ en dos números de Mach que
no se habían explorado. Salidas: `out/e1bis_nL.csv`, `out/e1bis_fits.csv`, `out/e1bis_result.csv`, `out/e1bis_manifest.json`,
`out/e1bis_run.log`.
