# Fase 2 — Tercer addendum al pre-registro (2026-11-18): correcciones tras la revisión externa

**No cambia ningún criterio** de `FASE2_PREREG.md` ni edita los addenda 1 y 2. Se escribe tras la revisión del revisor
externo (ChatGPT), que concluye que está justificado ejecutar la medición sin cambiar criterios, y **antes** de ejecutarla.

## 1. Corrección de una frase
En el resumen del addendum 2 y en el mensaje que lo acompañó se dijo que «los diez casos con γ > 0 están en transición».
**No es exacto:** con umbral de fase 2, γ = 10 da régimen k⁶ (R > 10) en w = 2 y w = 3, y transición en w = 4, 6 y 8. La
tabla del addendum 2 ya lo recogía; la frase resumen, no.

## 2. Término que faltaba en la expansión (señalado por el revisor)
Desarrollando Ω = cp·√(1 + p²/(4c²) + γp⁴), el segundo orden de la raíz, −⅛·(p²/(4c²))², da en Ω un término −p⁵/(128c³)
que se había omitido. El coeficiente de p⁵ pasa de γc/2 a a = γc/2 − 1/(128c³).
Análisis de sensibilidad: `residual_anatomy/fase2/regimen_analitico_v2.py`, salida `out/fase2_regimen_analitico_v2.log`.

| γ | Sin el término (addendum 2) | Con el término | Cambio |
|---|---|---|---|
| 10 | R = 3,9–12,6; p_eff = 0,189–0,193 | R = 3,9–12,5; p_eff = 0,189–0,193 | R < 0,5 %; p_eff igual a tres decimales |
| 1 | R = 0,8–2,9; p_eff = 0,184–0,187 | R = 0,75–2,8; p_eff = 0,184–0,187 | R ≈ 2–3 %; p_eff igual a tres decimales |
| 0 | p_eff = 0,333 (exacto) | p_eff = 0,40 / 0,47 / 0,74 (umbral 0,5 / 1 / 2) | **grande** |

- **γ = 1 y 10: el término no cambia nada.** El revisor tenía razón en que faltaba, y el efecto es despreciable, como
  estimaba el asesor (cociente con el término séxtico ~1/(64γc⁴)).
- **γ = 0: la heurística no es fiable para el control.** Con γ = 0, la escala de la heurística es p\* ≈ 1,1–2,2, es decir
  p/(2c) ≈ 0,7–1,4: la expansión de √(1 + p²/(4c²)) no converge bien ahí, y truncarla en un orden u otro cambia mucho el
  exponente. **Esto no afecta al criterio K**, que no se basa en esta heurística sino en el benchmark v1.0 y en el ajuste de
  Finazzi y Parentani (p(0) ∈ [0,25; 0,42]). Pero significa que la heurística **no** predice el control; si K fallase, este
  sería un motivo plausible a examinar, no una excusa para cambiarlo.
- Para γ = 10 la escala relevante es menor (p\* ≈ 0,5–0,9, p/(2c) ≈ 0,3–0,6), donde la expansión es más fiable; por eso la
  estimación ≈ 0,19 sí es estable.

## 3. Qué no cambia
Criterios K, P1 y P0; γ ∈ {0, 1, 10}; w ∈ {2, 3, 4, 6, 8}; malla y estimadores. La medición se ejecuta a continuación con
`measure_F2.py`, sin modificar (su I0 comprueba el pre-registro, el addendum 1, el propio script y `sextic4.py`).
