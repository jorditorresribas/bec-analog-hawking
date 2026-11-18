# Fase 2 — Resultado de los gates técnicos F0–F3 (2026-11-18)

Criterios: `docs/FASE2_DISENO.md` §5, en commit (`55cf395`) antes de implementar. Código: `residual_anatomy/fase2/sextic.py`
y `residual_anatomy/fase2/gates_F.py`. Salida completa: `out/fase2_gates_run.log`; datos: `out/fase2_F0.csv`,
`out/fase2_F1.csv`, `out/fase2_F2.csv`.

## Veredicto
| Gate | Veredicto | En una línea |
|---|---|---|
| **F0** fondo y Goldstone | **FALLA** | el fondo es consistente (control negativo: ×10⁴–10⁵), pero en mallas finas el residuo lo domina un suelo de redondeo ∝ N⁴ |
| **F1** homogéneo | **FALLA** | raíces, normas y recuento correctos; transmisión correcta con γ = 1. Los residuos de diferencias finitas y la pseudo-unitariedad en el horizonte, contaminados por el mismo suelo (con γ = 10 y 2,25·N₀, defecto 0,16) |
| **F2** límite γ → 0 | **PASA (a)** | con γ = 0, el driver nuevo reproduce v1.0 a **3,3·10⁻¹²** (tol. 10⁻⁸); la parte (b), informativa, no se ejecuta porque usa el horizonte con γ > 0 |
| **F3** no degeneración | **NO EJECUTADO** | depende de F0 |

**La medición física no se puede ejecutar con esta formulación.** Siguiendo la regla del encargo, ese camino se detiene
aquí.

## F0 en detalle
| N | residuo v1.0 | γ = 10, fondo consistente | γ = 10, fondo NO consistente (V₀) |
|---|---|---|---|
| 2000 | 1,505·10⁻⁶ | 1,505·10⁻⁶ | 0,142 |
| 4000 | 3,761·10⁻⁷ | 3,773·10⁻⁷ | 0,142 |
| 8000 | 9,40·10⁻⁸ | **5,20·10⁻⁷** | 0,142 |

- (c) **pasa**: con el fondo no consistente, el residuo es ~10⁵ veces mayor y no baja con dx. El fondo V_γ = V₀ − γgn⁗ es
  el correcto, y el gate detecta un fondo erróneo (la lección de la fase 1).
- (a) y (b) **fallan** en N = 8000 (y en N = 16 000 el residuo sube a 7,9·10⁻⁶).

## Diagnóstico (a posteriori; no cambia ningún veredicto)
`residual_anatomy/fase2/diag_F0_redondeo.py`, `out/fase2_diag_F0_redondeo.log`.

- Algebraicamente, M_γΦ_G = M_v1.0Φ_G, porque el bloque nuevo actúa sobre u + v = 0. Numéricamente, el bloque tiene entradas
  de orden γg·6/dx⁴ ≈ 10⁹ (γ = 10, dx = 0,01). Al ensamblar M_v1.0 + γg[[P, P], [−P, −P]] en doble precisión, las entradas
  de P se suman a entradas distintas de v1.0 en los bloques uu y uv y se redondean de forma distinta, así que la cancelación
  P√n − P√n deja un error ~ε·γg‖D₄‖ ∝ dx⁻⁴ ∝ N⁴.
- Pruebas: (1) el residuo crece como N⁴ (×15 de N = 8000 a 16 000; se esperaba ×16); (2) hacer el producto
  matriz-vector en long double **no** lo elimina (3,0·10⁻⁷ y 4,2·10⁻⁶), así que el error ya está en las entradas
  ensambladas, no en el producto; (3) en N ≤ 4000, donde el suelo está por debajo del error de discretización, el residuo
  coincide con el de v1.0 a 3·10⁻³.
- **El mismo suelo explica F1:** un ∂⁴ en diferencias finitas aplicado a ondas de k pequeño resta números de orden 1/dx⁴
  para obtener k⁴. En los residuos de las ondas planas aparece como un error que crece con N, y en el scattering del
  horizonte, el condicionamiento ∝ dx⁻⁴ lo amplifica (defecto de pseudo-unitariedad 1,0·10⁻⁴ → 2,1·10⁻⁴ → 0,16 con
  N₀ → 1,5N₀ → 2,25N₀, γ = 10).
- **No es un error de la física ni del fondo**: con γ = 1 y mallas moderadas, el modelo se comporta (transmisión 10⁻⁸,
  defecto en el horizonte 2,5–5,5·10⁻⁴). Es un límite de la discretización elegida.

## Opciones para seguir (decisión del usuario y del revisor; ninguna aplicada)
1. **Variables (s, d) = (u + v, u − v)** en el ensamblado: P aparece una sola vez, actuando sobre s, así que la cancelación
   del Goldstone es exacta. **Y mallas más gruesas con esténciles de cuarto orden** (dx ≈ 0,02–0,03: el suelo baja como dx⁴,
   y el error de discretización de las partes de segundo orden se compensa con el esténcil de cuarto orden). Contra: el
   operador de γ = 0 ya no es bit a bit el de v1.0. F2 pasaría a comparar en el límite del continuo (Richardson), no a 10⁻⁸.
   **Es la que recomiendo.**
2. **γ menor** (γ = 1: suelo 10 veces menor). Contra: a γ = 1, k⁶ no domina sobre k⁴ en las anchuras ensayadas (régimen de
   cruce), así que la predicción de exponente 1/5 no se pondría a prueba de forma limpia.
3. **Representación espectral (Fourier) de ∂⁴** sobre la parte suave. Contra: no local, matrices densas, y en la fase 1 los
   operadores no locales tuvieron artefactos de caja finita.
4. **Precisión extendida en el ensamblado y en el solve.** Contra: SciPy no tiene solve disperso en long double; habría que
   salir del stack actual.

Cualquiera de ellas exige reescribir los criterios de F0–F3 (F0′–F3′), ponerlos en commit antes de ejecutar y repetir los
gates. La medición física sigue bloqueada hasta que pasen.

## Lo que se ha calculado de γ > 0 (para que conste antes del pre-registro)
Durante la depuración del *matching* y en estos gates se calcularon: n_L en el horizonte w = 2 para γ = 1 y 10 en
ω ≈ 0,01, 0,05·T₀, 0,15·T₀ y 0,27·T₀ (para el defecto de pseudo-unitariedad), y casos homogéneos. **No se ha calculado
ninguna temperatura T_up, ningún κ̄ ni ninguna anchura d\*, ni se ha ajustado ningún exponente.**

---

## Gates revisados F0′–F3′ (formulación de `docs/FASE2_DISENO.md` §6; criterios en commit `39f7607`)
Código: `residual_anatomy/fase2/sextic4.py` y `gates_F4.py` (en commit `80014e4` antes de ejecutar).

### Ejecución 1: F1′ FALLA (se paró ahí, como exigen los criterios)
`out/fase2p_gates_run1_FALLA.log`. El criterio (a) falló: las raíces de la red se alejaban del continuo como h⁻⁴ (error de
2·10⁻⁹ a 4·10⁻⁸ con γ = 10). **Causa: un error mío de implementación.** Los coeficientes de los esténciles en alta precisión
se creaban al importar el módulo con los 15 dígitos por defecto de mpmath, y su redondeo, multiplicado por γg/h⁴ ≈ 10⁷,
desplazaba las raíces. Se comprobó que la referencia del continuo era exacta (mpmath frente a `np.roots`: 10⁻¹⁶) y que
`polyroots` convergía (error estimado 10⁻⁶¹). Corrección en `f64c58e`: coeficientes como racionales exactos convertidos
dentro de `workdps`. **Criterios sin cambios.**

### Ejecución 2: TODOS PASAN (`out/fase2p_gates_run.log`, `out/fase2p_F*.csv`)
| Gate | Veredicto | Datos |
|---|---|---|
| F1′ homogéneo | **PASA** | raíces de la red → continuo con orden 4,00 (γ = 1 y 10, sub/supersónico); normas exactas; 10 raíces y 5 «away» por lado; \|T − 1\| ≤ 1,6·10⁻¹¹ y resto ≤ 2·10⁻¹⁵; suelo 0 |
| F0′ Goldstone | **PASA** | residuo idéntico para γ = 0, 1, 10 (dif. 0); 1,03·10⁻⁹ → 3,27·10⁻¹⁰ → 6,47·10⁻¹¹, orden 4,00; fondo NO consistente: 1,4·10⁻² (γ = 1) y 0,14 (γ = 10), 10⁷–10⁸ veces mayor |
| F1′(e) horizonte | **PASA** | defecto de pseudo-unitariedad 4,8·10⁻⁷ → 3,0·10⁻⁸ (γ = 1) y 2,1·10⁻⁷ → 1,3·10⁻⁸ (γ = 10), orden 4,0; suelo de n_L ≤ 3,7·10⁻⁹ |
| F2′ γ → 0 | **PASA** | γ = 0: orden 4,00; límite del continuo frente al de v1.0: 5,9·10⁻⁹ (tol. 10⁻⁶). Informativo: n_L(γ) − n_L(0) = −1,0·10⁻⁴ / −2,7·10⁻⁴ / −6,7·10⁻⁴ para γ = 0,01 / 0,03 / 0,1 (orden en γ 0,8–0,9) |
| F3′ no degeneración | **PASA** | γ = 10: D_min = 3,1·10⁻⁴ en s = 1,303, ruido 5,9·10⁻⁸ → 5340 × ruido; control λ = 4: D_min = 6·10⁻⁷ en s = 1,9999 (degenerado, como debe) |

**Advertencia para la medición física (ver `FASE2_PREREG_ADDENDUM.md`):** la no degeneración es clara numéricamente
(5340 veces el ruido), pero **pequeña en magnitud**: con γ = 10 y w = 2, el espectro n_L(ω) coincide a 3·10⁻⁴ con un
reescalado del modelo γ = 0 (s ≈ 1,30). El efecto de la forma k⁶ sobre el espectro, en esta anchura, es casi un cambio de
escala efectivo.
