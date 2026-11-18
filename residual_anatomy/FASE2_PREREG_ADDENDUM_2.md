# Fase 2 — Segundo addendum al pre-registro (2026-11-18): en qué régimen está el diseño

**No cambia ningún criterio** de `FASE2_PREREG.md` (K, P1, P0, tolerancias, γ y anchuras). Declara, antes de cualquier
medición, una estimación **analítica** del régimen (sin resolver BdG): `residual_anatomy/fase2/regimen_analitico.py`,
salida `out/fase2_regimen_analitico.log`. Se hace a petición del asesor.

## Método
Con la dispersión mixta c²k² + k⁴/4 + γc²k⁶, Ω ≈ cp + p³/(8c) + γcp⁵/2, y la fase dispersiva del argumento de Coutant y
Parentani (diseño §3) es p³/(24cκ) + γcp⁵/(10κ). La escala del horizonte ensanchado p\* es aquella en la que la fase llega a
un umbral de orden 1 (se prueban 0,5, 1 y 2). El peso relativo de los dos términos en p\* es R = 2,4·γc_H²·p\*². Criterio: k⁴
si R < 0,1, k⁶ si R > 10, transición en medio. Con c = c_H = 0,794 y κ = κ_local de las cinco anchuras.

## Resultado
| γ | R en p\* (w = 2 → 8) | Régimen | Exponente efectivo de la estimación |
|---|---|---|---|
| 0 | 0 | k⁴ en las cinco anchuras | 0,333 (exacto, por construcción) |
| 1 | 0,8–2,9 | **transición** en las cinco | 0,184–0,187 |
| 10 | 3,9–12,6 | **transición** (con umbral 2: k⁶ en w = 2 y 3, transición en el resto) | 0,189–0,193 |

- **El diseño no está en el régimen k⁶ puro:** en γ = 10, el término k⁶ pesa entre 4 y 13 veces más que el k⁴ en la escala
  relevante, no ≫ 10 en todas las anchuras. Y R decrece con w: las anchuras grandes están más cerca del régimen cuártico.
- **Aun así, la propia estimación da un exponente efectivo ≈ 0,19** con γ = 10 (y también ≈ 0,18 con γ = 1), casi
  insensible al umbral y dentro de la ventana de P1 (1/5 ± 0,05). En la transición, la estimación no da un valor intermedio
  entre 1/3 y 1/5, sino ≈ 0,19.
- **Cautelas:**
  - es la misma heurística que la predicción del diseño (§3), no una comprobación independiente;
  - los prefactores de orden 1 no están fijados;
  - en la transición, la medición real puede depender de detalles que la estimación no recoge;
  - F3′ ya mostró que en w = 2 el efecto del k⁶ es casi un reescalado (addendum 1).
  - Por todo ello, un resultado INCONCLUSO es plausible y estaría bien explicado.
- Que γ = 1 dé, según la estimación, un exponente parecido al de γ = 10 es una predicción de la heurística: el pre-registro
  lo trata como informativo, y así sigue.

## Qué no cambia
Los criterios K, P1 y P0, γ ∈ {0, 1, 10}, w ∈ {2, 3, 4, 6, 8} y la malla. No se ha calculado ninguna temperatura, κ̄, d\* ni
exponente medido. El script de la medición (`measure_F2.py`) comprueba en I0 el pre-registro y el addendum 1, pero no este
addendum 2, que solo declara.
