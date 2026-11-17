# Gate externo E1 — reproducir Larré–Recati–Carusotto–Pavloff 2012 (*waterfall*) con el operador v1.0

**Estado:** BORRADOR DE DISEÑO EN PAPEL para el referee. **No implementado ni ejecutado.**
**Referencia:** P.-É. Larré, A. Recati, I. Carusotto, N. Pavloff, PRA 85, 013621 (2012), arXiv:1110.4464v2
(texto original leído: §II.C, Ecs. 16–17, 44, 56–57, Apéndice A.2).
**Motivo:** validar el instrumento contra un resultado **analítico exacto de un grupo independiente**, en la misma
clase de sistema que el nuestro (g constante, potencial externo, solución GP exacta). λ = 1: el posible artefacto de
`NOVELTY_AUDIT.md` §3 no le afecta.

## 1. Sistema (LRCP §II.C), unidades ℏ = m = 1, g = 1, n_u = 1 → c_u = 1, ξ_u = 1
- Potencial: U(x) = −U0·Θ(x), con U0/(g n_u) = m_u²/2 + 1/(2m_u²) − 1 (Ec. 17).
- Aguas arriba (x < 0), medio solitón oscuro (Ec. 11 con x0 = 0, sin θ = m_u):
  φ_u(x) = cos θ·tanh(x cos θ/ξ_u) − i sin θ, de modo que n(x) = n_u[cos²θ·tanh²(x cos θ) + sin²θ].
- Aguas abajo (x > 0): uniforme, n_d = n_u·m_u², V_d = V_u/m_u², m_d = 1/m_u² (Ec. 16).
- Continuidad n·v = n_u·V_u = j. Flujo hacia +x, aguas arriba a la izquierda (**mismo sentido que v1.0**).
- **Comprobación previa (integridad):** el residuo de la GP estacionaria con este n, v, U debe ser O(dx²) en el
  interior. n es C¹ en x = 0 y n'' salta (U es discontinuo): hay que **medir** el orden de convergencia, no suponerlo.

## 2. Observables y valores de referencia (analíticos)
LRCP Ec. 44: a baja ω, S_{u,d2}(ω) = f/√ω + h√ω + …, luego |S_{u,d2}|² = |f|²/ω + 2Re(f*h) + O(ω).
Γ = −4Re(f*h) y T_H = |f|²/Γ (Ec. 56). |S_{u,d2}|² **es nuestro n_L** (canal u saliente aguas arriba al inyectar el
entrante de norma negativa d2), con modos normalizados a corriente unidad.

| m_u | m_d | **a ≡ |f|²** (A8) | **b ≡ 2Re(f*h)** (A14) | Γ | T_H/(mc_u²) (Ec. 57) |
|---|---|---|---|---|---|
| 0.5 | 4.000 | 0.13173 | −0.48980 | 0.97959 | 0.13448 |
| 0.7 | 2.041 | 0.06691 | −0.49797 | 0.99595 | 0.06718 |
| 0.3 | 11.11 | 0.18153 | −0.45317 | 0.90634 | 0.20028 |

(Calculado por mí con A8/A14. Γ(0.5) coincide con el 0.980 del texto y T_H con la Ec. 57 a 5 cifras.)
Observables **primarios**: a y b, extraídos de ω·n_L(ω) = a + b·ω + c·ω² (+ d·ω³) en una ventana baja.
Se comparan **a y b directamente**: son menos dependientes del estimador que T_H y Γ. T_H y Γ se derivan de ellos.

## 3. Instrumento
- **Operador:** `fd_bdg_operator.build_MBdG` **sin modificar** (g = 1, V(x) escalón, μ = g n_u + V_u²/2).
- **Driver nuevo** (`gates/gate_E1_lrcp.py`): `mode_matching.solve_inject` tiene fijos `bg`, ω, μ y las asíntotas,
  así que se reimplementa el *matching* modal con fondo genérico. El driver usa `plane_wave_mode` y `quartic_roots`
  de v1.0 y conserva la normalización por corriente.
- **Control de regresión obligatorio:** el driver nuevo, alimentado con el fondo tanh de v1.0, debe reproducir
  n_L del benchmark congelado (`out/fp_benchmark.csv` / `solve_inject`) con error ≤ 1e−8 relativo en ω = 0.1.
  Si no, STOP: el driver no es v1.0.
- Rejilla: x = 0 **en un nodo**. Dominio [−L, L].

## 4. Pruebas (fijadas antes de ejecutar)
**Integridad (si falla cualquiera: STOP sin veredicto)**
- I1: regresión del driver (≤ 1e−8).
- I2: residuo GP del fondo LRCP en el interior; se informa el orden observado.
- I3 — Convergencia de pseudo-unitariedad (redacción del referee, 2026-11-16; sustituye al umbral 1e−6 anterior).
  Se evaluará el defecto de pseudo-unitariedad de la columna d2, |Σ_out η|a|² − η_in|, en al menos tres resoluciones,
  con razón de refinamiento 2. Se exigirá que tienda a cero bajo refinamiento y que el orden observado sea compatible
  con el orden del esquema. Se informarán el defecto en cada resolución, el orden observado y la extrapolación.
  No se impondrá un umbral absoluto 1e−6 sobre la malla más fina sin justificar que está por debajo del error de
  discretización alcanzable. (Definición del defecto: idéntica a `flux_out − eta_in` de `correlations.scatter_field`
  y a `pu` de los scripts FP de v1.0; ver la nota de equivalencia del 2026-11-16 en la memoria del proyecto.)
- I4: número de canales = 3 entrantes / 3 salientes (ω < Ω) en todas las ω usadas.

**Ciencia (para cada m_u ∈ {0.5, 0.7}; 0.3 solo como ampliación si las dos primeras pasan)**
- E1-a: |a_num − a_th| / a_th ≤ **tol_a**.
- E1-b: |b_num − b_th| / |b_th| ≤ **tol_b**.
- Convergencia: Richardson en N (3 mallas, razón 2) y en L (3 dominios), con el orden observado informado.
- Robustez del estimador: a y b con dos ventanas de ω y dos órdenes de polinomio. El *spread* entra en el presupuesto.

**Tolerancias (justificación física):** el resultado de LRCP es **exacto** para este modelo (BdG sin más
aproximaciones), así que la discrepancia solo puede venir de la numérica. Se fija:
- presupuesto numérico δ = Richardson (N) ⊕ dominio (L) ⊕ *spread* del estimador;
- PASS si |Δ| ≤ max(3δ, piso) **y** δ ≤ techo;
- **a:** piso 1e−3 relativo, techo 5e−3;
- **b** (coeficiente subdominante, más difícil de extraer): piso 5e−3, techo 2e−2.
Si δ supera el techo: **INCONCLUSO** (el instrumento no tiene resolución), nunca PASS.

**Clasificación:** PASS (a y b en las dos m_u) / PARCIAL (a pasa, b inconcluso) / FALLA (cualquier |Δ| > max(3δ, piso)
con δ ≤ techo) / INCONCLUSO. Una FALLA se diagnostica (convenciones de canal, normalización, discontinuidad de U)
**sin ajustar nada a posteriori**.

## 5. Persistencia (Control B)
`out/e1_lrcp_*.csv` (n_L(ω) por m_u, N, L; ajustes; integridad) y `out/e1_lrcp_manifest.json` (sha256 del driver,
de `fd_bdg_operator.py`, de este diseño; entorno/BLAS; parámetros; exit code).

## 6. Puntos abiertos para el referee
1. ¿Basta con m_u ∈ {0.5, 0.7}, o se exige también 0.3 (m_d = 11: aguas abajo muy diluido, ξ_d grande)?
2. ¿Tolerancias aceptables? Alternativa: fijarlas tras un *dry run* de resolución, pero **sin mirar a_th**.
3. ¿Comparar también la curva completa n_L(ω) hasta Ω? LRCP no da figura del *waterfall* (Fig. 6 es la delta),
   así que solo habría referencia analítica a baja ω.

4. **Comprobaciones previas pedidas por el referee (2026-11-16), a documentar en papel antes de implementar:**
   (a) correspondencia entre las convenciones de f, h (LRCP Ec. 44) y nuestra S_{u,d2} (orden y etiquetado de
   canales, signo de la fase, dependencia en ω^{−1/2}); (b) normalización de modos y flujo (corriente unidad en ambos);
   (c) definición exacta de a y b como coeficientes de ω·n_L(ω); (d) sensibilidad de a y b a la ventana infrarroja.

*Ventana ω y valores de N/L: se fijan tras la revisión, antes de implementar.*

---

## 7. Comprobaciones previas resueltas y parámetros fijados (2026-11-16, ANTES de implementar)
Decisión del usuario: E1 se hace en una sola sesión (límite: un día). Si a y b no cuadran, el resultado es
INCONCLUSO documentado y se cierra igualmente. El referee no ha visto este §7 antes de la ejecución; se le entrega después.

**7.1 Convenciones (lectura del original LRCP, §III.A–E y Ec. 44)**
- (a) **f, h frente a S_{u,d2}.** LRCP Ec. 44: S_{u,d2} = f/√ε_u + h·√ε_u + O(ε_u^{3/2}), con ε_u = ℏω/(mc_u²). En
  nuestras unidades (ℏ = m = g = n_u = 1 → c_u = 1) **ε_u = ω**. Entonces |S_{u,d2}|² = |f|²/ω + 2Re(f*h) + O(ω).
  a y b se definen a partir de |S|², así que **no dependen de la convención de fase** de los modos (una fase global o
  analítica en ω no los cambia). El canal: u|out es el saliente aguas arriba (norma +, v_g < 0) y d2|in el entrante
  aguas abajo de norma negativa. Coinciden con nuestro n_L (salida L de norma + al inyectar la entrada R de norma −).
  Mismo sentido del flujo (aguas arriba a la izquierda, V > 0).
- (b) **Normalización.** LRCP normaliza todos los modos reales a corriente |J_ℓ| = 1 (Ecs. 30, 35) y su S es
  η-unitaria (Ec. 43). La nuestra normaliza a |J| = 1 con J = (k+v)|U|² + (k−v)|V|². Las dos corrientes difieren a lo
  sumo en un factor constante común a todos los modos, y |S_{ν′ν}|² es un cociente de flujos: **es invariante**.
- (c) **Definición de a y b.** ω·n_L(ω) = a + b·ω + c·ω² + d·ω³; a_th = |f|² (A8), b_th = 2Re(f*h) (A14).
- (d) **Ventana infrarroja**, fijada antes de ejecutar y escalada con T_H(m_u) de la Ec. 57:
  - W1 (referencia): ω/T_H ∈ [0.1, 0.6], 12 puntos geométricos;
  - W2 (control): ω/T_H ∈ [0.05, 0.3], 12 puntos.
  Ajuste de referencia: W1, orden 3. *Spread* del estimador: máx |Δ| frente a {W1 orden 2, W2 orden 2, W2 orden 3}.
  Umbral Ω (fin del canal d2, calculado): 1.764 (m_u = 0.5) y 0.616 (m_u = 0.7); las ventanas quedan muy por debajo.

**7.2 Fondo e interfaz.** Rejilla con x = 0 en un nodo (N impar). U(0) = −U0/2: con este valor, la fila de la
interfaz aproxima el promedio de las ecuaciones a ambos lados (n″ salta en x = 0). Es una convención fijada a priori;
su efecto entra en el orden observado.

**7.3 Mallas y dominio.** dx ∈ {0.02, 0.01, 0.005} (razón 2) con L = 30; L ∈ {20, 30, 40} con dx = 0.01.
Valor central: Richardson en dx a L = 30 con el orden observado (si no es estimable, el valor de la malla fina).
δ = √(δ_N² + δ_L² + δ_est²):
- δ_N = |a_fina − a_Rich|;
- δ_L = máx |a(L) − a(30)|;
- δ_est = *spread* de 7.1(d).

**7.4 Precisión de I3 (orden compatible).** El esquema es de orden 2, pero la discontinuidad de U puede degradarlo a
orden 1 (§1). Se fija como compatible **p_obs ∈ [0.8, 2.4]**, además de que el defecto decrezca de forma monótona.
