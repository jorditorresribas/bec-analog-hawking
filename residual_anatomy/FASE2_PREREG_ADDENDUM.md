# Fase 2 — Addendum al pre-registro (2026-11-18): lo que han revelado los gates F0′–F3′

`FASE2_PREREG.md` está en commit (`870cb60`) y **no se edita**. Este addendum declara información obtenida **después** de
ese commit y **antes** de cualquier medición física, para que conste quién sabía qué. No cambia la predicción ni los
criterios K, P1 y P0.

1. **Los gates revisados F0′–F3′ pasan** (`residual_anatomy/FASE2_GATES_RESULT.md`), en la formulación de
   `docs/FASE2_DISENO.md` §6: (s, d), cuarto orden, modos exactos de la red, h ∈ {0,04; 0,03; 0,02}. Se cumple, por tanto,
   la condición 2 del pre-registro. La condición 1 (revisión del revisor externo) sigue pendiente; la 3 (script de medición
   en commit antes de ejecutar), también.
2. **Visto en F3′ (w = 2):** con γ = 10, n_L(ω) en la banda del benchmark coincide a **3·10⁻⁴** relativo con
   n_L^{(γ=0, w/s)}(s·ω) para **s ≈ 1,30**. Es decir: en w = 2, el término k⁶ actúa casi como un reescalado efectivo de
   la anchura, w → w/1,30. Implica, de forma aproximada, que T_up(γ = 10, w = 2) es la del modelo γ = 0 con w ≈ 1,54, en
   unidades reescaladas. **No se ha calculado ninguna T_up, κ̄, d\* ni exponente**; esto es una inferencia a partir de
   un número del gate.
3. **Lo que eso anticipa:** la predicción P1 (p → 1/5) solo puede cumplirse si ese factor de escala efectivo s **depende de
   κ** (de w) del modo que predice el argumento de la fase. Si s fuese el mismo para todas las anchuras, el exponente no
   cambiaría (P0). Ninguno de los dos casos está decidido por lo visto, que es una sola anchura.
4. Visto en F2′ (b): n_L(γ) − n_L(0) en w = 2, ω = 0,1·T₀, para γ = 0,01; 0,03; 0,1 (−1,0·10⁻⁴; −2,7·10⁻⁴; −6,7·10⁻⁴).
5. Malla para la medición: la de los gates revisados (h ∈ {0,04; 0,03; 0,02}, Richardson con p = 4). Sustituye a la frase
   del §3 del pre-registro «los que fijen los gates F0′–F3′ revisados». Convergencia exigida para T_up: ≤ 10⁻⁴ (sin cambios).
