# Auditoría de novedad — ¿está ya respondida nuestra pregunta?

**Fecha:** 2026-11-16 · **Estado:** BORRADOR para el referee externo. No autoriza ejecutar nada.
**Contexto:** P2 cerrado en clase C (`P2_RESULT.md`). El roadmap P4a–P6 queda en pausa hasta cerrar esta auditoría.
**Pregunta auditada (formulación vigente):** ¿persiste el colapso de n_L(ω) cuando se varía la dispersión UV
(coeficiente cuártico λ) manteniendo fijos el fondo n(x), v(x), c(x) y el acoplo de Bogoliubov G = gn?

Etiquetas: **PROBADO** (verificado aquí, en el texto original o algebraicamente), **HIPÓTESIS**, **ABIERTO**.
Fuentes leídas en su versión original de arXiv (texto completo): MP2009 (0905.3634v4), FP2012 (1202.6015v2),
CPF2012 (1108.1821v3), LRCP2012 (1110.4464v2). Las demás solo por su resumen (se indica).

---

## 1. Qué establecen los tres artículos pedidos (PROBADO sobre el texto original)

### 1.1 Macher & Parentani 2009 — PRA 80, 043601 (arXiv:0905.3634)
- BdG 1D estacionario **sin aproximaciones adicionales** (Ecs. 12–16), con transformación de Bogoliubov 3×3.
- **La ecuación de las fluctuaciones relativas φ (Ψ = Ψ0(1+φ)) solo depende de c(x) y v(x)**: el potencial
  externo y la presión cuántica del fondo **desaparecen** (Ec. 12, con T_ρ = −(ℏ²/2mρ0)∂ρ0∂). Es relevante para §3.
- Perfil: c+v = c0·D·sign(x)·tanh^{1/n}[(|κx|/Dc0)^n] con n = 2 (Ec. 21); reparto c/v con q (Ec. 22);
  κ = d(c+v)/dx en el horizonte (Ec. 20, flujo hacia la izquierda, v < 0); Λ = √2c0/ξ0, λ = Λ/κ.
- **Fig. 7** (confirmado): f_ω y T_ω/T_H frente a ω/κ para (D, λ) = (0.1, 70), (0.4, 18), (0.7, 10), q = 0.3;
  ω_max/κ = 1.17, 2.2, 2.54; T0/T_H = 1.08, 0.996, 1.0004; |ΔT0| ≃ 0.1 %.
- **§V.C.2 (clave para la novedad):** "relajamos la ligadura (73) y permitimos **valores arbitrarios de λ para un D
  dado**". Es decir, **variar la escala de dispersión con el perfil c(x), v(x) fijo**. Conclusión: la desviación
  respecto a Hawking la controla **ω_max/κ**, no λ (Fig. 8).
- **Fig. 10** (confirmado): T_ω/T_H para D = 0.4 y (q, λ) = (0, 19.54), (0.5, 17.15), (1, 15.47), con ω_max/κ = 2.2
  fijo: el reparto c/v (q) cambia apreciablemente el espectro; q = 0.5 minimiza la mezcla u–v.

### 1.2 Finazzi & Parentani 2012 — PRD 85, 124027 (arXiv:1202.6015)
- **Mismo BdG** (Ec. 6: (ω−vk)² = c²k² + c0⁴k⁴/Λ²), perfiles c/v con q (Ec. 8).
- **Barrido en el plano (Λ/κ, ω_max/κ)** (Figs. 1–2): "ω_max/κ es el parámetro espectral más relevante **en todo
  el plano**". Hawking (T_eff/T_H dentro del 10 %) para ω_max/κ ≳ 0.4; régimen escalonado para ω_max/κ ≲ 0.1.
- En **todos los casos estudiados** n_ω diverge como 1/ω (T_eff = lím ω·n_ω finito), incluso para Λ/κ < 1.
  **Corrección (referee, 2026-11-16):** que diverja como 1/ω **no fija el coeficiente** A(λ) en n_L ~ A(λ)/ω,
  ni descarta supresiones a frecuencia finita. Lo que sí aporta FP2012 es el **mapa medido** de T_eff/T_H
  en el plano (Λ/κ, ω_max/κ) para sus perfiles (Fig. 1).
- Régimen escalonado: T_step = (3√3/8)·ω_max·[1 − (q − 1/6)D + …] (Ec. 15). **Finito, no nulo.**
- Corrección menor al resumen del referee: la longitud d_ξ = d·c0(κΛ²)^{−1/3} con d ≃ 2^{2/3} es de
  **Finazzi–Parentani 2011** (PRD 83, 084010, arXiv:1012.1556, "broadening the horizon"). FP2012 la cita y
  da su escalado ξ^{2/3}κ^{−1/3}. Es la referencia que reproduce nuestro benchmark v1.0.

### 1.3 Coutant, Parentani & Finazzi 2012 — PRD 85, 024021 (arXiv:1108.1821)
- Matriz S analítica por fórmulas de conexión. La robustez del espectro y de la función de dos puntos se explica por
  la extensión de la región cercana al horizonte. Las primeras desviaciones decrecen como la inversa de esa extensión.
- **Fases:** las de los coeficientes de Bogoliubov importan para la función de dos puntos (agujeros negros y blancos,
  láser). Además, **arg(β/α) por sí sola no es observable**: depende de la elección de fase de los modos. Solo son
  físicas ciertas combinaciones. Esto respalda el uso de invariantes plaqueta en P2.
- Estudiar las fases de S **no es novedoso** por sí solo (coincido con el referee).

**Veredicto sobre las afirmaciones del referee:** todas confirmadas en el original, salvo la atribución de d_ξ (§1.2).

---

## 2. Comparación con nuestra pregunta

| Elemento de nuestra pregunta | ¿Cubierto? | Dónde |
|---|---|---|
| Robustez de n_ω frente a dispersión (general) | **Sí** | Unruh 1995, Corley–Jacobson 1996, MP2009, FP2011/2012, CPF2012 |
| Variar la escala de dispersión con c(x), v(x) fijos, en BdG | **Sí** | MP2009 §V.C.2 (Fig. 8); FP2012 Figs. 1–2 |
| Qué parámetro gobierna la respuesta | **Sí**: ω_max/κ (y d_ξ / κ̄) | FP2012; FP2011; Coutant–Parentani 2014 (PRD 90, 121501, solo resumen) |
| Dispersiones superlumínicas genéricas (convexas) | **Parcialmente** | FP2012 (correspondencia entre dispersiones); Del Porro–Liberati–Schneider 2026 (arXiv:2604.24861, solo resumen) |
| Fases de S y correlaciones | **Sí** | CPF2012; MP2009 nota 3 |
| Factores grises y mezcla u–v a baja ω en BdG completo | **Sí** | Coutant–Weinfurtner 2018 (PRD 97, 025006, solo resumen); MP2009 Figs. 10–11 |
| Configuración tipo *waterfall* (solución GP exacta) | **Sí**, analítico | Larré–Recati–Carusotto–Pavloff 2012 (PRA 85, 013621) |
| **"Colapso" de n_L al variar λ** | **Contradicho** por la literatura | ver §3 |
| Dispersión variable con **G fijo** (operador no local, fuera de BdG) | **No encontrado** | ABIERTO (§4) |

---

## 3. Hallazgo principal: el "colapso" de n_L choca con la literatura, y hay una causa candidata

### 3.1 El knob antiguo es el barrido de MP/FP (PROBADO, algebraico)
`build_op_lambda` (`correlations.py:89`) usa t = √λ/2 y G = gn/√λ, con c² = gn fijo y convección de coeficiente 1.
**Precisión (referee):** la equivalencia exige cambiar **a la vez** la masa y el acoplo: m′ = 1/√λ **y**
g′n = m′c² = c²/√λ (que es justo G = gn/√λ). Con eso la dispersión homogénea es (ω−vk)² = c²k² + λk⁴/4, o sea
Λ_eff = Λ/√λ. **Solo la parte homogénea** queda probada así. Que el operador inhomogéneo corresponda a esa teoría,
y que el fondo sea estacionario para ella, es lo que examinan §3.3 y el Apéndice A. Por tanto, el knob antiguo es el barrido "λ arbitrario a D fijo"
de MP2009 §V.C.2 y FP2012 Fig. 1. La conclusión del 2026-11-15 ("el knob mezcla dispersión con reescalado del acoplo")
es correcta en la forma, pero ese "reescalado" es justo el que hace la literatura al cambiar Λ.

### 3.2 La literatura predice que no hay colapso en nuestro régimen (PROBADO para su modelo; inferencia para el nuestro)
- Nuestros números (w = 2): c_H² = n_H = 0.630 → Λ = 2c_H² = 1.26; κ_local = 0.270 → Λ/κ ≈ 4.7; ω_max = 0.468 →
  **ω_max/κ ≈ 1.73**: régimen de Hawking según FP2012. Coherente con v1.0 (T_up/κ_local = 0.957).
- Con λ ∈ [0.25, 4], ω_max/κ ∈ ~[0.9, 3.5] (escala ∝ Λ_eff). En el mapa de FP2012 Fig. 1, esa zona tiene
  T_eff/T_H de ~0.9 a 1.0, que varía suavemente, **para sus perfiles**. Es una tensión fuerte con un colapso
  de T_eff a 0, pero **no una imposibilidad demostrada** para nuestro observable y perfil (corrección del referee:
  la literatura no excluye por sí sola un colapso de A(λ)).

### 3.3 Causa candidata: el modo de Goldstone queda roto para λ ≠ 1
**PROBADO (álgebra sobre el código):** en `build_op_lambda` la presión cuántica del fondo se mantiene
independiente de λ (Qp = V − μ + gn + v²/2). Aplicado al modo de fase del condensado (u, v) = (√n, −√n) a ω = 0:

    fila u:  −t(√n)'' + Qp·√n − i(v(√n)' + v'√n/2)  =  (½ − t)(√n)''  =  ½(1 − √λ)(√n)''
    (el término convectivo se anula porque n·v = j es constante; Qp·√n = ½(√n)'' por la GP con m = 1)

Ese residuo es **distinto de cero para λ ≠ 1** allí donde varía la densidad (alrededor del horizonte). Es decir:
el fondo **no es un estado estacionario** de la dinámica con masa m′, y se linealiza alrededor de algo que no es
solución. El BdG consistente con masa m′ lleva Qp_λ = √λ·Qp, que coincide con el de MP2009, donde Qp ni siquiera aparece.

**PROBADO (álgebra sobre la especificación):** `NONLOCAL_SPEC.md` §2 y §3 (y por tanto la realización Padé
K_P) mantiene la misma decisión ("Qp λ-indep, fondo impuesto"). Su residuo es (ε_λ(−i∂) − ε_1(−i∂))√n = (K_P − D/2)√n,
también ≠ 0 para λ ≠ 1. La especificación es además internamente inconsistente: dice que V_ext "se ajusta para
sostener el fondo", y ajustarlo es precisamente lo que cambia Qp.

**HIPÓTESIS (mecanismo, sin verificar):** el residuo actúa como un "término de masa" localizado para el campo de
fase. A baja ω ese término no se apaga mientras la parte cinética sí lo hace (∝ ω²), así que dispersa con fuerza los
fonones de baja frecuencia y ω·n_L → 0, es decir, **T_eff → 0**. El mecanismo encaja con las **tres firmas de
artefacto** anotadas el 2026-11-15:
1. colapso de T a 0 frente a un efecto físico esperado del orden del %;
2. colapso para los dos signos de λ − 1 (cualquier residuo ≠ 0 lo produce);
3. λ = 1 es el único punto sano (único λ con residuo nulo).

También explica por qué los tests de modos homogéneos no lo detectaron: en las zonas asintóticas (√n)'' = 0
(lección 6 del proyecto: "validar donde importa").

**Consecuencias si se confirma (HIPÓTESIS):**
- El colapso de n_L no sería física de "esta familia", sino un artefacto de linealización.
- La pregunta "¿persiste el colapso con G fijo?" **perdería su motivación**.
- El instrumento Padé en construcción (P4a–P6) **heredaría el mismo defecto** tal como está especificado.
- Los resultados de robustez con λ de CORRELATIONS §7–8 usaron el mismo operador, así que habría que revisarlos.
- El benchmark v1.0 (λ = 1) **no se ve afectado** (residuo nulo).

### 3.4 Prueba discriminante propuesta (barata; requiere diseño, referee y autorización)
- **T1 (Goldstone):** ‖M_λ·(√n, −√n)‖ / ‖M_λ‖‖√n‖ para λ ∈ {0.25, 0.5, 1, 2, 4}, en las dos versiones de Qp
  (independiente de λ, y √λ·Qp). Predicción: ∝ |1 − √λ| frente a ~dx² (solo discretización).
- **T2 (refutación):** n_L(ω; λ) con Qp_λ = √λ·Qp (BdG consistente = barrido MP/FP), mismas mallas que el Paso 0.
  - Si T_eff(λ) varía de forma suave y compatible con FP2012 Fig. 1 (gobernada por ω_max/κ), el colapso era artefacto.
  - Si el colapso persiste con Qp consistente, la hipótesis queda **refutada** y hay que buscar otra causa.
- Coste estimado (no medido): minutos con el solver por *matching* modal de v1.0, sin tocar v1.0
  (script nuevo que importa `build_MBdG`/funciones existentes o reimplementa el operador).

---

## 4. ¿Qué novedad sigue siendo defendible? (ABIERTO)
1. **No defendible:** "la radiación de Hawking es robusta frente a la dispersión"; "estudiamos las fases de S";
   "variamos la escala de dispersión a métrica fija" (MP2009, FP2012).
2. **Defendible solo como validación numérica:** nuestro benchmark v1.0 (T0 = κ̄/2π ≤ 0.75 %, FP2011) y la
   metodología de des-rotación e invariantes plaqueta (P2). No es física nueva.
3. **Posible novedad acotada (sin establecer):** un modelo **con dispersión deformada y acoplo G fijo**, fuera de BdG,
   **formulado de forma que conserve el modo de Goldstone** (por ejemplo, en fluctuaciones relativas, con un operador
   cinético no local que anule las constantes, como T_ρ de MP2009). Pregunta: a igual ω_max/κ, ¿cambia n_L la
   estructura interna u–v (|V/U| a k finito, factores grises) respecto a BdG? La literatura (FP2012, Del Porro et al.
   2026) sugiere que dominaría ω_max/κ, así que el efecto esperado es pequeño. Antes de afirmar novedad hay que revisar
   el texto completo de Coutant–Parentani 2014, Coutant–Weinfurtner 2018, Del Porro et al. 2026 y Michel–Parentani–Zegers
   2016 ("no-hair", PRD 93, 065039; no leído).
4. **Recomendación:** no invertir en el solver Padé (P4a–P6) hasta resolver §3 con T1/T2.

---

## 5. Gate externo con v1.0 — evaluación de candidatos

| | **LRCP 2012 *waterfall*** (propuesta nueva) | **MP2009 Fig. 7** (propuesta del referee) |
|---|---|---|
| Sistema | g constante, escalón U(x) = −U0Θ(x), solución GP **exacta** (medio solitón oscuro aguas arriba) | c(x), v(x) prescritos con q = 0.3: requiere **g(x) variable** (Feshbach) |
| Parecido a nuestro modelo | Muy alto (misma clase: g fija, potencial externo) | Medio |
| Observable | \|S_{u,d2}\|² = **nuestro n_L** (canal upstream desde la entrada de norma −) | n_ω, f_ω, T_ω/T_H |
| Referencia cuantitativa | **Analítica:** k_BT_H/(mc_u²) = ½(1−m_u⁴)^{3/2}/[(2+m_u²)(1+2m_u²)] (Ec. 57); con m_u = 0.5: m_d = 4, Γ ≃ 0.980, T_H/(mc_u²) ≈ 0.1345 (cálculo mío con la Ec. 57) | Numérica: T0/T_H = 1.08, 0.996, 1.0004; curvas leídas de la figura |
| Ventajas | Fórmula cerrada en función de m_u (curva Fig. 5) → compara una familia, no un punto; grupo independiente (Pavloff) | Tres perfiles y varias escalas dispersivas |
| Dificultades | Potencial discontinuo (ψ'' salta en x = 0: puede degradar el orden de convergencia); hay que separar Γ y T_H con el desarrollo a baja ω (Γ·[T/ω − ½]) | Perfil c/v con g(x), flujo con signo contrario (espejo x → −x), n = 2 |
| Viabilidad con v1.0 (sin tocarlo) | `build_MBdG` sirve tal cual (g = 1, V escalón); hace falta un **driver nuevo**: `mode_matching.solve` tiene fijos `bg`, ω, μ y las asíntotas | `build_MBdG` solo depende de g·n y v: valdría pasando n_arg = c²/g y V coherente; también requiere un **driver nuevo** |

**Recomendación (HIPÓTESIS de diseño, para el referee):** usar **LRCP 2012 como gate externo principal** (T_H y Γ
frente a m_u, contra la Ec. 57 y la Fig. 5) y **MP2009 Fig. 7 como secundario**. En ambos casos "con el solver v1.0"
significa reutilizar su operador sin modificarlo, con un driver nuevo auditado. No es reejecutar v1.0 tal cual.
Queda pendiente leer el Apéndice A de LRCP (f_{u,d2}, h_{u,d2} analíticos) antes de fijar las tolerancias.

---

## 6. Orden de trabajo propuesto (cada paso: diseño en papel → referee → autorización)
1. **T1/T2 (§3.4):** decide si el colapso es artefacto. Es lo de mayor impacto y lo más barato.
2. **Gate externo LRCP (§5):** diseño con entradas, observables (T_H, Γ, \|S_{u,d2}\|²(ω)), malla, tolerancias y criterios.
3. Según T2: (a) si el colapso era artefacto, reformular la pregunta física (§4.3) o cerrar la línea; (b) si persiste,
   reanudar el diagnóstico antes de cualquier P4a.

P2 sigue cerrado en clase C. P4a sigue bloqueado. No se ha ejecutado nada en esta auditoría.

### Referencias
- J. Macher, R. Parentani, PRA 80, 043601 (2009), arXiv:0905.3634
- S. Finazzi, R. Parentani, PRD 85, 124027 (2012), arXiv:1202.6015
- A. Coutant, R. Parentani, S. Finazzi, PRD 85, 024021 (2012), arXiv:1108.1821
- S. Finazzi, R. Parentani, PRD 83, 084010 (2011), arXiv:1012.1556
- P.-É. Larré, A. Recati, I. Carusotto, N. Pavloff, PRA 85, 013621 (2012), arXiv:1110.4464
- A. Coutant, R. Parentani, PRD 90, 121501 (2014), arXiv:1402.2514 (solo resumen)
- A. Coutant, S. Weinfurtner, PRD 97, 025006 (2018), arXiv:1707.09664 (solo resumen)
- F. Del Porro, S. Liberati, M. Schneider, arXiv:2604.24861 (2026, preprint; solo resumen)

---

## Apéndice A — Derivación del residuo de Goldstone (pedido por el referee, 2026-11-16)

Convenciones del código: ℏ = 1, campo con la fase del condensado extraída, Φ_G = (ψ0, −ψ0*) → (√n, −√n) en esta
convención (comprobado: es el núcleo de `build_MBdG` con λ = 1, ver A.3).

**A.1 GPE con masa m′ y acoplo g′ (teoría consistente con el knob λ).** iΨ_t = [−(1/2m′)∂² + V_λ + g′|Ψ|²]Ψ.
Estacionario Ψ0 = √n·e^{iθ(x)}·e^{−iμ′t}, θ′ = m′v:
- parte imaginaria: (n v)′ = 0 → n v = j (la continuidad no depende de m′);
- parte real: −(1/2m′)(√n)″ + (m′v²/2 + V_λ + g′n − μ′)√n = 0 → **V_λ + g′n + m′v²/2 − μ′ = Qp_λ**, con
  **Qp_λ ≡ (1/2m′)(√n)″/√n**.
Con m′ = 1/√λ y g′n = c²/√λ (para mantener c² = g′n/m′ = gn): **Qp_λ = √λ·Qp_1**, Qp_1 = ½(√n)″/√n.
Para sostener el **mismo** n(x) hay que cambiar el potencial: V_λ = μ′ − g′n − m′v²/2 + √λ·Qp_1.

**A.2 Linealización (BdG) alrededor de ese fondo.** Con t = 1/(2m′) = √λ/2 y 2tθ′ = v:
L_± = −t(∂ ± iθ′)² + 2g′n + V_λ − μ′ = −t∂² ∓ i(v∂ + v′/2) + [g′n] + [g′n + tθ′² + V_λ − μ′]
    = −t∂² ∓ i(v∂ + v′/2) + G_λ + Qp_λ,   G_λ = g′n = gn/√λ,   M = [[L_+, G_λ], [−G_λ, −L_−]].
Es la estructura de `build_op_lambda` (t, G, convección), **salvo** que el código usa **Qp_1** en lugar de Qp_λ.
Los términos de flujo (m′v²/2) y de potencial (V_λ, μ′) quedan absorbidos exactamente en Qp_λ; no hay otros.

**A.3 Aplicación a Φ_G = (√n, −√n).** Fila u: L_+√n − G√n = −t(√n)″ − i(v(√n)′ + v′√n/2) + Qp_code·√n.
- Convección: v(√n)′ + v′√n/2 = (nv)′/(2√n) = 0 (exacto en el continuo).
- Qp_code = V − μ + gn + v²/2 = Qp_1 **por construcción del fondo** (`bg`: V = μ − gn − v²/2 + Qp analítico).
- Luego fila u = −t(√n)″ + ½(√n)″ = **½(1 − √λ)(√n)″**. Fila v: −G√n + L_−√n = el mismo valor (con +i·0).
- Con Qp_λ (teoría consistente A.1–A.2): fila u = −t(√n)″ + t(√n)″ = 0. **El residuo completo, flujo y potencial
  incluidos, es ½(1−√λ)(√n)″ en ambas componentes, y se anula para λ = 1 o con Qp_λ.**

**A.4 Operador no local (Weyl) y Padé con G fijo.** Aquí no hay una GPE subyacente definida (el símbolo ε_λ(x,k)
depende de G(x) = gn(x)), y la especificación solo define el operador de fluctuación. La condición mínima para que
Φ_G sea modo cero, con la misma convección y la misma diagonal G + Qp_1, es K_λ√n + Qp_1√n = 0. Como
K_1√n + Qp_1√n = 0 (K_1 = −½∂², salvo discretización), **el residuo es (K_λ − K_1)√n**, distinto de cero para
λ ≠ 1 donde n varía. Para Padé: K_P − D/2 = (λ−1)/4·D(2G+D)⁻¹D − (λ−1)²/16·D²(2G+D)⁻³D², con D = −∂².
Ver la comprobación numérica en §A.5 (diagnóstico `diag_goldstone.py`).

**A.5 Comprobación numérica (diagnóstico, no gate)** — `residual_anatomy/diag_goldstone.py` v2;
`out/diag_goldstone.csv` + `out/diag_goldstone_manifest.json` (+ log `out/diag_goldstone_run.log`).
Residuo ‖M·Φ_G‖/‖Φ_G‖ medido **solo en el interior** |x| ≤ L − 10. La v1, conservada como
`out/diag_goldstone_v1_bordes.*`, estaba dominada por las filas de borde de la FD (crecía como N^1.5 incluso en v1.0),
que el *matching* modal sustituye; error de medida mío, corregido antes de interpretar.

| Operador | λ = 1 | λ ≠ 1 | Escalado en N | Frente a la predicción |
|---|---|---|---|---|
| A v1.0 `build_MBdG` | 1.5e−6 → 3.8e−7 → 9.4e−8 (N = 2000/4000/8000) | — | ×4 por duplicación: **O(dx²)** | — |
| B `build_op_lambda` (Qp fijo) | igual que A | 2.56e−3 (0.25), 1.50e−3 (0.5), 2.12e−3 (2), 5.12e−3 (4) | **independiente de N** | ½\|1−√λ\|‖(√n)″‖/‖√n‖: cociente **1.0000** |
| B′ ídem con Qp_λ = √λ·Qp | igual que A | ~1e−7 a N = 8000 | O(dx²) | se anula como predice A.1–A.3 |
| D Padé 2º orden (Qp fijo) | 1.4e−5 → 3.4e−6 (N = 400/800): O(dx²) | 8.5e−4, 5.6e−4, 1.0e−3, 2.6e−3 | estable | (K_P − K_1)√n: cociente 0.997–1.005 a N = 800 |
| C Weyl no local (Qp fijo) | **1.8e−3, independiente de N** | 1.1e−3 a 4.3e−3 | — | cociente 1.0–2.4 |
| Convección sobre √n | 1.4e−6 → 8.5e−8 | — | O(dx²) | 0 en el continuo |

**Lectura:**
- **PROBADO numéricamente (dentro del alcance documentado):** `build_op_lambda` rompe el modo de Goldstone para
  λ ≠ 1 con exactamente el residuo de A.3. Con Qp_λ = √λ·Qp se recupera el modo cero hasta la discretización.
- **PROBADO numéricamente:** el operador **Padé** (Kpade2) rompe el modo cero con el residuo (K_P − K_1)√n.
- **ABIERTO:** el operador Weyl de referencia (gate0) **no** satisface Goldstone ni siquiera con λ = 1 (residuo
  1.8e−3, que no baja con N). Causa no diagnosticada; candidata: truncación del núcleo o efecto del dominio finito
  (núcleo con colas). No se puede usar como referencia limpia en este test.
- **Sigue siendo HIPÓTESIS:** que esta ruptura **cause** el colapso de n_L. Lo decide T2 (n_L(λ) con Qp_λ).

**A.6 Diagnóstico T2 — n_L(λ) con fondo consistente** (`residual_anatomy/diag_T2_qp_consistent.py`;
`out/diag_T2.csv`, `out/diag_T2_nL.csv`, `out/diag_T2_manifest.json`, log `out/diag_T2_run.log`).
Mismas condiciones que el Paso 0 (w = 2, L = 40, N = 9600 y 14400, 9 frecuencias ω/(κ/2π) ∈ [0.05, 0.27], estimador
de T0 por extrapolación cuadrática, ruta `scatter_field`). Dos variantes **emparejadas** del mismo operador:
Qp fijo (el del Paso 0) y Qp_λ = √λ·Qp. Integridad: el envoltorio con λ = 1 es idéntico al original (0.0) y la
regresión al benchmark congelado da 4.6e−11 → PASS.

| λ | T_eff/T_H, Qp fijo (Paso 0) | T_eff/T_H, **Qp_λ (consistente)** | Δ entre N | *spread* entre órdenes (Qp_λ) |
|---|---|---|---|---|
| 0.25 | −0.006 | **0.97900** | 5.7e−5 | 3e−5 |
| 0.5 | −0.027 | **0.97139** | 2.1e−5 | 2e−5 |
| 1 | 0.95724 | **0.95724** | 8.9e−6 | 3e−5 |
| 2 | −0.053 | **0.93241** | 6.4e−6 | 3e−5 |
| 4 | −0.017 | **0.89175** | 2.2e−6 | 2e−5 |

(Columna «Δ entre N» = |T_eff/T_H(N = 9600) − T_eff/T_H(N = 14400)| con Qp_λ, de `out/diag_T2.csv`; *corregida en la auditoría de números del 2026-11-18*: antes decía 5e−5, 3e−5, 1e−5, <1e−5 y 1e−5. Ninguna conclusión cambia.)

(T_H ≡ κ_local/2π. Pseudo-unitaridad de la columna d2: ≤ 6.7e−3 a N = 9600 y ≤ 3.0e−3 a N = 14400, mejora con N.)

**Lectura:**
- **PROBADO (dentro del alcance documentado: w = 2, este perfil, esta banda, λ ∈ [0.25, 4]):**
  - con Qp fijo se reproduce el colapso del Paso 0: n_L **se anula** al bajar ω (n_L ∝ ω; p. ej. λ = 4: 0.16 en ω/T_H = 0.05,
    frente a 17 con Qp_λ), así que T_eff → 0;
  - con el fondo consistente, el colapso **desaparece**: n_L recupera la divergencia 1/ω y T_eff/T_H varía
    **suave y monótonamente** de 0.979 a 0.892, convergido en N (≤ 6e−5) y estable frente al orden de extrapolación.
- **Conclusión:** el "colapso de n_L" fue un **artefacto** del fondo no estacionario (modo de Goldstone roto, A.3–A.5),
  no física. Los resultados de robustez en λ de CORRELATIONS §7–8 usaron el mismo operador y deben marcarse como
  afectados (no reevaluados aquí).
- **Compatible con la literatura (inferencia, no comparación cuantitativa):** aumentar λ reduce Λ_eff y ω_max/κ, y
  T_eff/T_H baja suavemente, en el rango ~0.9–1.0 del mapa de FP2012 Fig. 1. Una comparación cuantitativa exigiría
  sus perfiles.
- **Consecuencia para la novedad (§4):** la pregunta "¿persiste el colapso con G fijo?" **pierde su objeto**. Lo que
  queda es la robustez conocida (MP2009, FP2012), ahora reproducida por nuestro instrumento con el fondo correcto.

**A.7 Revisión provisional del referee (2026-11-16) — alcance que se adopta**
(Revisión hecha sobre la transcripción, **no** sobre scripts/CSV/manifiestos; pendiente revisión del paquete.)
- Goldstone: probado **en el interior** examinado (|x| ≤ L − 10), no en todo el dominio ni en sus condiciones de contorno.
- T2: convergencia probada en **N** (9600 → 14400), **no en L** (dominio fijo L = 40). La estabilidad frente al orden de
  extrapolación no prueba que la extrapolación infrarroja sea correcta.
- Causalidad: se acepta "el colapso es un artefacto del operador antiguo **dentro del alcance ensayado**"; **no** está
  demostrado que la ruptura de Goldstone sea la **única** causa posible.
- Padé (formulación anterior): ruptura de Goldstone con evidencia numérica dentro del test descrito; P4a–P6 siguen
  bloqueados en esa formulación. Weyl: ABIERTO e independiente de T2.
- Recomendación del referee: **no abrir más diagnósticos sobre el colapso**; revisar el paquete T2 y el diseño LRCP.
