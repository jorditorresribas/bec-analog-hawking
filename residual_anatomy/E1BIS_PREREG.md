# Gate E1-bis — criterios pre-registrados (Larré et al. 2012, coeficiente a con extrapolación controlada)

**Estado:** PRE-REGISTRO. **No ejecutado.** Este archivo, `residual_anatomy/gates/gate_E1bis_lrcp.py` y
`residual_anatomy/diag_conditioning.py` deben estar en un commit **antes** de la ejecución. El script lo comprueba:
si alguno de los tres no está en git o tiene cambios sin commit, se niega a ejecutar (exit 2) y no calcula nada.
El manifiesto de salida guarda el hash del commit, que demuestra que los criterios son anteriores al resultado.

Redactado: 2026-11-17. Autor de los criterios: el usuario (a propuesta del revisor); implementación: Claude Code.

## 1. Pregunta
¿Reproduce el solver v1.0 el coeficiente a = lim_{ω→0} ω·n_L de la Ec. A8 de Larré, Recati, Carusotto y Pavloff
(PRA 85, 013621, 2012) a mejor que 10⁻⁶ relativo, una vez controlados el redondeo (solve refinado), la
discretización (Richardson) y el sesgo de la extrapolación en ω?

E1 dio 2,4·10⁻⁵ (m_u = 0,5) y 2,9·10⁻⁵ (m_u = 0,7). El diagnóstico `diag_convergence.py` atribuyó ese error al
ajuste en ω (ventana W1 con orden 3). E1-bis comprueba esa atribución en configuraciones que no se han mirado.

## 2. Qué se ha visto ya (declaración honesta)
- **m_u = 0,5:** en `diag_convergence.py` (exploración a posteriori) se vio a con W2/W3 y orden 4 a 1,1·10⁻⁸ y 1,3·10⁻⁸
  de A8. La elección de W2 y del orden 4 está **informada por esa exploración**. Por eso m_u = 0,5 se ejecuta como
  **réplica, no confirmatoria**: su resultado se publica pero **no cuenta para el veredicto**.
- **m_u = 0,7:** E1 calculó n_L en W2 con dx = 0,02/0,01/0,005 (`out/e1_lrcp_nL.csv`) y lo ajustó solo con órdenes
  2 y 3. Nunca se ha ajustado con orden 4, ni con dx = 0,0025, ni con solve refinado. **Confirmatoria.**
- **m_u = 0,6:** nunca se ha calculado nada. **Confirmatoria.**
- Para verificar el script antes del commit se ejecutó su ruta principal **solo con m_u = 0,5** (ventanas W1 y W2 de E1,
  dx = 0,02/0,01/0,005, orden 4, salida al directorio temporal, no al repositorio). Ahí se vio por primera vez b con
  orden 4 en m_u = 0,5 (3,5·10⁻⁵). Ninguna configuración de m_u = 0,6 o 0,7 se ha calculado.
  Después de esa prueba solo se añadió esta declaración y esta referencia: no cambió ningún umbral ni criterio. La versión
  anterior a la prueba se conserva literal en `residual_anatomy/E1BIS_PREREG_pre-prueba.md`; un `diff` entre los dos
  archivos muestra únicamente estas líneas añadidas.
- Antes de este pre-registro se han ajustado dos tolerancias de tests después de un fallo (declarado en
  `tests/validate_lambda.py` y `residual_anatomy/diag_scaling_equivalence.py`). Este gate no tiene esa opción: si
  falla, falla.

## 3. Configuración fijada
| Parámetro | Valor |
|---|---|
| Fondo | waterfall LRCP de `gate_E1_lrcp.bg_lrcp` (sin cambios), L = 30 |
| m_u confirmatorias | **0,6 y 0,7** |
| m_u de réplica (no cuenta) | 0,5 |
| Mallas | dx = 0,01; 0,005; 0,0025 (`gate_E1_lrcp.grid`) |
| Ventana principal | **W2 = [0,05; 0,3]·T_H**, 12 ω en progresión geométrica (T_H de la Ec. 57) |
| Ventana de control | W3 = [0,025; 0,15]·T_H, 12 ω |
| Ajuste | polinomio de **orden 4** de ω·n_L en ω; a = término constante, b = término lineal (`gate_E1_lrcp.fit_ab`) |
| Solve | `diag_conditioning.solve(..., 'eqref')`: equilibrado de Ruiz + 4 refinamientos con residuo en long double, sobre el sistema ensamblado por `diag_conditioning.assemble` (réplica bit a bit de `solve_d2`). **El operador y la discretización son los de v1.0, sin cambios.** |
| Richardson | p = 2 **declarado**: a_R = a(0,0025) + [a(0,0025) − a(0,005)]/3 |
| Teoría | `gate_E1_lrcp.lrcp_theory` (A8, A14, Ec. 57), sin cambios |

## 4. Criterios de integridad (si falla alguno: STOP, sin veredicto, exit 2)
- **I0** — los tres archivos citados arriba están en git y sin cambios respecto a HEAD.
- **I1** — en m_u = 0,6, dx = 0,01 y la primera ω de W2, `eqref` y `spsolve` dan el mismo n_L a ≤ 1·10⁻⁷ relativo
  (comprueba que se resuelve el mismo sistema).
- **I2** — en cada (m_u, dx, ω): canales Lin = 1, Rin = 2, Rin_neg = 1, Lout = 1, Rout = 2 (como E1 v2).
- **I3** — orden observado p_obs = log₂(|a(0,01) − a(0,005)| / |a(0,005) − a(0,0025)|) en [1,5; 2,5] para la ventana
  principal en cada m_u (si no, el p = 2 declarado no está justificado).

## 5. Criterios de veredicto (solo m_u = 0,6 y 0,7)
Sea ε_a = |a_R − a_A8| / a_A8 con la ventana principal W2.

| Criterio | Condición | Papel |
|---|---|---|
| **C1 (principal)** | ε_a ≤ **1·10⁻⁶** en m_u = 0,6 **y** en m_u = 0,7 | decide PASS / FAIL |
| **C2 (control de ventana)** | \|a_R(W2) − a_R(W3)\| / a_A8 ≤ 1·10⁻⁶ en las dos m_u | decide PASS / FAIL junto con C1 |
| C3 (nivel fuerte) | ε_a ≤ 1·10⁻⁷ en las dos m_u | solo se informa |
| C4 (secundario, b) | \|b_R − b_A14\| / \|b_A14\| ≤ 1·10⁻³ en las dos m_u (b_R con el mismo Richardson) | veredicto propio, no cambia C1/C2 |

**PASS** = C1 y C2. **FAIL** = cualquiera de los dos no se cumple. No hay segunda ejecución con otros parámetros:
cualquier cambio exige un E1-ter con su propio pre-registro en commit.

**Justificación de la tolerancia 1·10⁻⁶:** 24 veces por debajo del error de E1 (2,4·10⁻⁵), así que un PASS es una
mejora real de la validación; y 100 veces por encima de lo visto en m_u = 0,5 (1,1·10⁻⁸), para dejar margen a la
dependencia en m_u, que no se conoce. Está informada por la exploración de §2, y se declara.

## 6. Qué se afirmará según el resultado
- **PASS:** «Con redondeo, discretización y extrapolación en ω controlados, el coeficiente a coincide con la Ec. A8
  de Larré et al. a ≤ 10⁻⁶ en dos números de Mach no explorados (criterios en commit antes de ejecutar).» Más C3 y C4
  tal como salgan.
- **FAIL:** se publica igual, con la tabla. Significa que la atribución del error de E1 al ajuste en ω no se
  generaliza (o que A8 está mal transcrita para m_u ≠ 0,5) y hay que investigarlo. El veredicto de E1 no cambia.

## 7. Salidas (el script se niega a sobrescribirlas)
`out/e1bis_nL.csv` (n_L por punto), `out/e1bis_fits.csv` (a, b por m_u, ventana y dx), `out/e1bis_result.csv`
(criterios), `out/e1bis_manifest.json` (commit, SHA-256 de este archivo, del script, de `diag_conditioning.py`,
`gate_E1_lrcp.py`, `fd_bdg_operator.py`, versiones, tiempo). Ejecución, desde la raíz del repositorio:
`python3 residual_anatomy/gates/gate_E1bis_lrcp.py 2>&1 | tee out/e1bis_run.log`.

Antes del commit se puede comprobar el script **sin tocar ninguna configuración del gate** con
`python3 residual_anatomy/gates/gate_E1bis_lrcp.py --selftest`: reproduce el ajuste ya publicado de E1
(m_u = 0,5, W1, orden 3, dx = 0,02) y no escribe nada.
