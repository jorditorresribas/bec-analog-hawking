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
