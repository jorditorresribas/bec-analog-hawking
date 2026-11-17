# EXPEDIENTE TÉCNICO PARA AUDITORÍA CIENTÍFICA INDEPENDIENTE — `bec-gpe-solver`

**Propósito:** reconstruir el proyecto desde sus archivos y resultados reales para que un auditor
independiente (modelo de razonamiento distinto, y/o ChatGPT) pueda refutarlo contrastándolo con el
código. **No es una narración optimista del progreso.** Fase 1 de la auditoría adversarial de 3 fases
(Fase 1 = este documento; Fase 2 = auditoría Claude de mayor capacidad; Fase 3 = ChatGPT), propuesta
2026-11-16.

**Autor de la Fase 1:** Claude Code (la misma instancia que implementó el Bloque 2). **Sesgo declarado:**
construí B2.1–B2.4d; un auditor debe asumir que puedo haber heredado errores propios. **Dos modelos de
Claude NO son dos verificaciones independientes** (pueden compartir sesgos); por eso la Fase 3 externa.

**Método y límites de este expediente (honestos):**
- Elaborado por **inspección read-only** de archivos. **No se ejecutó ni modificó código** para redactarlo.
- Entorno congelado (`PROVENANCE.txt`, 2026-11-14): Python 3.9.6, numpy 2.0.2, scipy 1.13.1, sympy 1.14.0,
  Darwin x86_64. v1.0 con SHA-256 en `PROVENANCE.txt` / `FP_REVIEW_MANIFEST.txt`.
- **LIMITACIÓN DE TRAZABILIDAD CRÍTICA:** los **gates del Bloque 2 imprimen a stdout**; sus resultados
  numéricos **NO están persistidos en el repo**. Los números citados aquí para B2.1–B2.4d provienen de
  **ejecuciones de la sesión 2026-11-15/16**, capturadas en archivos de tarea **efímeros** del scratchpad
  (se borran). Son **reproducibles re-ejecutando el gate**, pero un auditor NO puede verificarlos contra un
  artefacto del repo. Esto se marca como **[stdout-sesión]** vs **[archivo-repo]** en todo el documento.

---

## 1. PREGUNTA CIENTÍFICA CENTRAL

¿Puede este proyecto calcular de forma **fiable** la ocupación Hawking `n_L(λ)` con **acoplo G fijo**, y
determinar si el comportamiento observado es **robusto frente a modificaciones de la dispersión UV**?

Subyacente (`NONLOCAL_SPEC.md §7`): el knob `build_op_lambda` (usado en trabajo previo) mostró cierta
robustez/colapso de `n_L`, **pero reescala G=gn/√λ** → no aísla la dispersión UV. El operador no local
`ε_λ(−i∂)` con **G fijo** sí la aísla. La pregunta: **¿el colapso de n_L venía de reescalar G, o persiste
con G fijo?** Responderla exige un **solver de scattering estabilizado** para el operador no local
(realización Padé) → medir `n_L(λ)`, λ∈[0.5,2].

**Estado de la respuesta: NO alcanzada.** Todo el Bloque 2 construye y valida el *instrumento* (solver);
la física `n_L(λ)` **no se ha medido**.

---

## 2. FORMULACIÓN FÍSICA Y MATEMÁTICA

### 2.1 Modelo (BEC 1D, gravedad análoga)
GPE 1D, ℏ=m=g=1. Fondo "waterfall" tanh: `n(x)=n_u+(n_d−n_u)·½(1+tanh(x/w))`, `v=j/n`, `c=√(gn)`; Mach
upstream 0.5, downstream 2; horizonte sónico donde c=v (x_H≈0.46). Parámetros: g=1, μ=1.125, j=0.5,
w=2, L=60 (dominio congelado x∈[−60,60], 8192 pts). **[archivo-repo: `out/waterfall_state.npz`]**

### 2.2 Dispersión de Bogoliubov y el knob de dispersión UV
`Ω²(k)=ε(k)[ε(k)+2G]`. BEC estándar: ε=k²/2, G=gn, c²=gn → `Ω²=c²k²+k⁴/4`.
Knob no local (G FIJO): `ε_λ(k)=−G+√(G²+c²k²+(λ/4)k⁴)` → `Ω²=c²k²+(λ/4)k⁴` ∀λ **[V, dif~1e−16, NONLOCAL_SPEC]**.
λ=1 ⇒ ε_1(k)=k²/2 EXACTO (c²=G) → operador local → **ancla de regresión** al benchmark FP.

### 2.3 Operador BdG y conservación
`L_± = ε_λ(−i∂) ∓ i(v∂+v'/2) + (G+Qp)`, `M=[[L_+,G],[−G,−L_−]]`. Norma `ρ=|u|²−|v|²`, métrica
`η=diag(1,−1)`. Conservación por **pseudo-hermiticidad** `M†η=ηM` (no hay corriente local J(x) simple para
operador no local). Tres niveles de control (global / balance regional / canales asintóticos), explícitos en
`NONLOCAL_SPEC §4`.

---

## 3. v1.0 — BENCHMARK FINAZZI–PARENTANI (CONGELADO)

**Afirmación acotada (`FP_BENCHMARK.md`):** para la familia tanh y las convenciones dadas, v1.0 reproduce
`T₀=κ̄/2π` con desviación **≤0.75%** en 5 anchuras w∈{2,3,4,6,8}, usando el observable FP `n_ω=|β_ω|²`
(ocupación del canal Hawking upstream `n_L`). **[archivo-repo: `reference/fp_benchmark_reference.csv`,
tabla en README_FP_BENCHMARK §2; T_up/kbar∈[0.9926,1.0050]]**.

- Solver v1.0: `mode_matching.py` (CONGELADO) + `bdg_smatrix.py` (S 3×3 con test S†ηS=η).
- Concordancia transfer vs mode-matching: β a 1e−4, N_ω a 1e−5 (declarado en README_FP §8).
- Comprobador por tolerancia `FP_check_reproduction.py` (RTOL=1e−4, ATOL=1e−12).
- **Limitaciones declaradas (honestas, FP_BENCHMARK §9):** familia de perfiles ≠ la de FP; d≈1.58 no
  universal; **ecuaciones de FP leídas del HTML de arXiv, pendiente cotejar con el PDF** (signos/convenios);
  régimen débilmente dispersivo. **El umbral 1e−4 del comprobador es de REGRESIÓN, NO cota del error físico.**

**Para el auditor:** v1.0 es el ancla. Si v1.0 fuera incorrecto (p.ej. error de convención FP), todo lo que
cuelga de la regresión λ=1 heredaría el error. El cotejo FP-PDF **no está hecho**.

---

## 4. RUMBO#2 — OPERADOR PADÉ CON G FIJO Y SU VALIDACIÓN

### 4.1 Realización Padé (local, orden finito)
El operador no local se aproxima por una racional en k² → operador diferencial local tratable con FD disperso
(`NONLOCAL_SPEC §6, ruta 2`). Realización (de `OPERATOR_VALIDATION_CLOSURE.md`):
`K_P = ½D + c₁ DAD + c₂ D²A³D²`, `A=(2·diag(G)+D)⁻¹`, `D=−∂²`, `c₁=(λ−1)/4`, `c₂=−(λ−1)²/16`.
`L_±=K_P±C+(G+Qp)`, `M_eff=[[L_+,G],[−G,−L_−]]`. Realización Dirichlet por nodos fantasma.

### 4.2 Validación del operador — CERRADA Y APROBADA (2026-11-15), **dentro de alcance documentado**
`OPERATOR_VALIDATION_CLOSURE.md` + `gates/GATES_REGISTRY.md`. Gates cerrados: P2 (analítico), P3.1 (A³ vs
Green, O(dx²)), P3.2-A (Schur/equivalencia aumentada), P3.2-B1 (simetría euclídea), P3.2-B2+B2c (identidad
de Green discreta→continuo O(dx²)), P3.2-B3 (pseudo-hermiticidad **PROYECTADA** sobre subespacio suave,
O(dx²)), Gates 1/2a (dispersión/regresión λ=1 a nivel operador), **P3.4-C2** (balance dinámico regional).
**[archivo-repo: los `gate_P3_*.py` existen; sus RESULTADOS numéricos no están persistidos — ver §9]**

**Lo que NO cierra (reservas explícitas, CLOSURE §3):** NO pseudo-herm matricial exacta ni autoadjunción
continua; **C1 (pseudo-unitaridad convergente de S) NO validada → RESERVADA** al solver estabilizado
(`STABILIZED_SOLVER_RESERVATION.md`); NO identidad SBP formal; NO ausencia universal de espurios; NO
`n_L(λ)`. El baseline λ=1 por transferencia (v1.0) tenía floor ~5e−5 sin convergencia en dx y sin meseta en
X0 — motivo de reservar C1.

**Para el auditor:** "validado dentro del alcance documentado" es una etiqueta recurrente. Escrutar si el
alcance declarado coincide con lo que los gates realmente prueban (P3.2-B3 es PROYECTADA, no exacta; P3.4-C2
usa estado de PRUEBA, no paquete BdG propagante).

---

## 5. BLOQUE 2 — SOLVER DE SCATTERING ESTABILIZADO (esta sesión)

**Objetivo (`STABILIZED_SOLVER_RESERVATION.md`):** construir un solver numéricamente estable que produzca la
matriz S del BVP (baseline λ=1 local → Padé λ≠1), sin la inestabilidad del transfer ingenuo, superando
criterios a priori (convergencia dx, meseta X0, S†η_out S=η_in con métricas IN/OUT separadas, control del
evanescente, sin espurios, regresión vs v1.0).

### 5.1 Arquitectura (BVP global aumentado) — diseño `STABILIZED_SOLVER_BLOCK1.md` (congelado)
- Estado `Y=(u,u',v,v')`, `Y'=A(x)Y`. **Esquema de caja** (trapezoidal) → 4N ecuaciones interiores.
- Clasificación modal POR ω: dirección `sign(v_g)` (NO sign k), norma `sign(Ω)`, evanescente por `Im k`
  dependiente del extremo. Verificado contra `bdg_channels.py`/`bdg_smatrix.py` de v1.0.
- **Matching por amplitudes aumentadas** `Y_b−Ξ_b a_b=0` (Ξ nunca se invierte). Restricciones: `a_grow=0`
  (evanescente creciente excluido), amplitudes IN prescritas.
- Sistema `(4N+12)×(4N+12)` **cuadrado**; unicidad ⟺ `ker(𝓜_aug)={0}`.
- **S∈ℂ³ˣ³** (ventana de Hawking 3×3: 1 IN izq + 2 IN der; 3 OUT). Guard de estructura completa de canales.
- Representación física: `ℓ=ξ_u=1/√2`, `Ỹ=(u,ℓu',v,ℓv')` (`STABILIZED_SOLVER_BLOCK2_2.md`).

### 5.2 Relaciones de corriente/norma (verificadas numéricamente, 18 modos)
`J=v_g(|u|²−|v|²)` **[V, ratio=1.0000, stdout-sesión]**; `sign(J)=sign(v_g)·sign(Ω)`; métrica de Krein
`η=sign(Ω)` (NO sign J). **v1.0 ya lo hace bien.**

---

## 6. TABLA DE GATES — BLOQUE 2 (estado de validación)

| Gate | Archivo | Qué valida | Resultado clave | Estado | Fuente |
|---|---|---|---|---|---|
| **B2.1** | `gate_B2_1_algebraic.py` | ensamblado 𝓜_aug cuadrado (4N+12), rango pleno (ker trivial), matching consistente, guard ventana 3×3 | rango pleno N=20…120; residuo 8e−15; guard FAIL fuera de ventana | CERRADO | [stdout-sesión] |
| **B2.2** | `gate_B2_2_conditioning.py` | condicionamiento: crudo/Ruiz/física, σ_min/rango, η_back, κ_out,b/M (adjunto), Δ_S | cond~N¹ (p≈0.97), η_back~1e−16, acuerdo repr ~1e−14, rango RESUELTO | CERRADO (calibración) | [stdout-sesión] |
| **B2.3** | `gate_B2_3_convergence.py` | convergencia en dx de S (X0=12 fijo) | **O(dx²)** (Richardson p→2.00, 2 freq); Δ_S(N=2000)=1e−3 | CERRADO | [stdout-sesión] |
| **B2.4** | `gate_B2_4_X0_convergence.py` | convergencia/condicionamiento en X0 (8…28) | R_cond≈1 (cond benigno); κ_out "explota" (←ARTEFACTO, ver §7); Δ_S O(1) (←FASE, ver §7) | **ABIERTO** | [stdout-sesión] |
| **B2.4b** | `gate_B2_4b_phase_sensitivity.py` | reanálisis fase + sensibilidad | moduli |S| convergen; S_fis converge con residuo O(dx²); adjunto roto a X0 grande | ABIERTO (provisional) | [stdout-sesión] |
| **B2.4c** | `gate_B2_4c_svd_verify.py` | SVD verdad-terreno de σ_min vs estimador | **estimador CORRECTO** (coincide SVD 5 cifras); cond~6.7e4 benigno a X0=28 | VERIFICADO | [stdout-sesión] |
| **B2.4d** | `gate_B2_4d_adjoint_check.py` | reproducción rutas adjunto | **`la.solve(Mᴴ)` DEFECTUOSA** (r_adj 5e6 a X0=28); `lu_solve(trans=2)` correcta; contigüidad REFUTADA | VERIFICADO | [stdout-sesión] |

**Operador (RUMBO#2):** P2, P3.1, P3.2-A/B1/B2/B2c/B3, Gates 1/2a, P3.4-C2 → CERRADOS (§4.2)
**[archivos `gate_P3_*.py`, `p2_analytic.py`, `gate_A3_kernel.py`; resultados no persistidos]**.

---

## 7. LA SAGA B2.4 — NARRATIVA HONESTA (incluye MIS errores, sin ocultarlos)

Esta sección documenta, por instrucción explícita (no ocultar errores previos ni sustituir conclusiones sin
explicar por qué cambiaron), la secuencia completa de B2.4, incluidas **dos falsas alarmas** del autor.

1. **Barrido B2.4:** Δ_S(X0) salió **O(1)** (S parecía NO converger) y `κ_out,M` **"explotó"** (2.7e5→5.4e13),
   `δ_adj`→6e6. El autor interpretó (alarma #1) que había una no-convergencia física + sensibilidad extrema.
2. **B2.4b (reanálisis):** (A) los **moduli |S_ij| CONVERGEN** → el Δ_S O(1) era **FASE** (S_ext rotaba con
   `e^{ikX0}`; Ξ es X0-independiente, la variación viene de la propagación interior). (B) `S_fis` reconstruida
   (fases exactas, signo corregido por ChatGPT: `S_fis=D_out⁻¹ S_ext D_in`) converge con residuo ~1e−3 que
   crece lineal en X0 → compatible con usar **k del continuo vs k discreto** (O(dx²)), **hipótesis no
   demostrada**. (C) prueba de perturbación: a X0=28 el adjunto `‖y‖=4.79e8`, identidad `s≠y†b` → el autor
   interpretó (alarma #2) que el **estimador de σ_min había fallado y cond explotaba e^{2κX0}** (C1 NO
   resuelto).
3. **B2.4c (SVD verdad-terreno):** **DESMONTÓ la alarma #2.** σ_min por SVD = **idéntico** al estimador
   (6.92e−5 a X0=28), cond~6.7e4 **benigno**. El estimador NO falló. `‖y‖` correcto (vía lu_solve) = 108.
4. **B2.4d (reproducción del adjunto):** la causa era un **bug de solver**: `la.solve(Mrep.conj().T, crep)`
   (ruta de `g2.kappa_out`) devuelve un `y` que **NO resuelve** el sistema a N grande (r_adj=5.44e6 a X0=28);
   `lu_solve(trans=2)` da el correcto (r_adj~1e−14). **La hipótesis de contigüidad del autor fue REFUTADA**
   (la copia contigua da idéntico resultado defectuoso). Micro-mecanismo exacto: **NO pinned**.

**Conclusión honesta de la saga:** (i) el condicionamiento de la matriz es **benigno** (SVD-confirmado), el
estimador fiable, el BVP **no** muestra el blow-up de cond del transfer en el rango probado; (ii) **toda la
columna κ_out,M/δ_adj de B2.4 es ESPURIA** (ruta adjunta defectuosa) → debe descartarse; (iii) el Δ_S O(1)
era fase; (iv) el autor **sobre-alarmó dos veces**. La disciplina de verificación cruzada (SVD + reproducción)
cazó tanto el artefacto de fase REAL como las dos falsas alarmas.

**Contradicción documental a resolver (para el auditor):** `g2.kappa_out` (en `gate_B2_2_conditioning.py`,
CERRADO/auditado) usa `la.solve(Mrep.conj().T, crep)` — la ruta ahora demostrada **no fiable a N grande**. Por
tanto el `κ_out` reportado en la calibración B2.2 (a N≤1200) y en B2.4 debe tratarse como **sospechoso**; a
N=1200 el error ya era ~1.9e−5 (pequeño) pero crece con N. **B2.2 está congelado; no se ha corregido.**

---

## 8. TAXONOMÍA DE AFIRMACIONES (rigurosa)

### Demostrado analíticamente
- Dispersión `Ω²=c²k²+(λ/4)k⁴` del knob no local; reducción exacta λ=1→ε_1=k²/2 (NONLOCAL_SPEC, [V]).
- Recuento BC del BVP y unicidad `ker(𝓜_aug)={0}` (BLOCK1).
- Relaciones `J=v_g(|u|²−|v|²)`, `η=sign(Ω)`.

### Verificado numéricamente (en los casos ensayados)
- Operador: P2–P3.4-C2 **dentro de alcance documentado** (pseudo-herm PROYECTADA, no exacta).
- B2.1 (ensamblado/rango), B2.2 (cond~N¹, η_back máquina), **B2.3 (O(dx²) por Richardson, 2 freq)**.
- B2.4c: estimador σ_min = SVD (cond benigno, X0≤28). B2.4d: ruta adjunta defectuosa identificada.
- A (moduli |S| convergen en X0).

### Compatible con los datos (no demostrado)
- `S_fis` converge físicamente en X0 salvo residuo O(dx²) por **k continuo vs discreto** (B2.4b-B): hipótesis
  plausible, **no verificada** (faltaría control dx=0.012 del residuo de fase + usar k discreto).
- "El BVP evita el blow-up de condicionamiento de C1": **favorable y SVD-respaldado en X0≤28**, no probado
  para X0 mayores ni como teorema.

### Hipótesis pendiente de comprobación
- Convergencia COMPLEJA de S en X0 (no solo moduli) con fase consistente con la discretización.
- Que `n_L(λ)` sea medible con fiabilidad (NO abordado).
- Micro-causa del bug `la.solve(Mᴴ)`.

### No resuelto
- **C1 (pseudo-unitaridad convergente S†η_out S=η_in)**: no medida de forma convergente/fiable. Era el objeto
  del solver; sigue abierta.
- **Pregunta física `n_L(λ)`**: no medida.
- Cotejo de ecuaciones FP vs PDF (v1.0).

---

## 9. RESULTADOS REPRODUCIBLES vs NO REPRODUCIBLES

- **Reproducibles con artefacto en repo:** benchmark FP (`FP_check_reproduction.py` vs
  `reference/fp_benchmark_reference.csv`); dispersión/escalas v1.0; S 3×3 (`bdg_smatrix.py`).
- **Reproducibles SOLO re-ejecutando (sin artefacto persistido):** TODOS los resultados de B2.1–B2.4d y de los
  gates P3.* del operador. El código está; **los números no**. Un auditor debe re-ejecutar para verificarlos.
  → **Riesgo de trazabilidad alto:** las conclusiones del Bloque 2 descansan en salidas efímeras.
- **Descartado como espurio (no usar):** columna κ_out,M/δ_adj de B2.4/B2.4b (bug adjunto).

---

## 10. RIESGOS PRIORIZADOS (matemáticos / numéricos / de proceso)

1. **[ALTO] C1 sin resolver.** El objetivo del solver (pseudo-unitaridad convergente de S) **no está medido
   fiablemente**. Sin ella, `n_L(λ)` no es confiable.
2. **[ALTO] Ruta adjunta defectuosa en código auditado.** `g2.kappa_out` usa `la.solve(Mᴴ)`; contamina toda
   métrica de sensibilidad basada en el adjunto. Hay que migrar a `lu_solve(trans=2)` y revalidar.
3. **[ALTO] Resultados del Bloque 2 no persistidos.** Trazabilidad frágil; re-ejecución obligatoria para auditar.
4. **[MEDIO] Convergencia física de S en X0 no demostrada** (solo moduli; S_fis con residuo O(dx²) hipotético).
5. **[MEDIO] Pseudo-hermiticidad del operador es PROYECTADA**, no exacta; el defecto matricial O(1) persiste en
   modos de rejilla. ¿Afecta a observables físicos?
6. **[MEDIO] Micro-causa del bug `la.solve(Mᴴ)` no identificada** → riesgo de que la ruta `lu_solve` comparta
   un modo de fallo latente a N aún mayores.
7. **[MEDIO] v1.0: ecuaciones FP sin cotejar con el PDF**; el ancla de regresión podría tener un error de
   convención no detectado.
8. **[BAJO] d≈1.58 y familia de perfiles ≠ FP** (limitación de alcance ya declarada).

---

## 11. PREGUNTAS QUE EL AUDITOR DEBE INTENTAR REFUTAR

1. ¿Es correcta la **clasificación de canales** (in/out por v_g, norma por Ω, evanescente por extremo) y la
   estructura 3×3? ¿Hay un signo o convención que invalide `S`?
2. ¿El **matching aumentado** y las BC (`a_grow=0` + IN prescritas) producen la S física correcta, o hay un
   error de fase/origen espacial además del ya detectado?
3. ¿Es **legítima** la reconstrucción de fase `S_fis=D_out⁻¹ S_ext D_in` con **k del continuo**? ¿El residuo
   O(dx²) es realmente k-discreto, o enmascara una no-convergencia física?
4. ¿Es `lu_solve(trans=2)` **fiable** para el adjunto a N grande, o comparte el fallo de `la.solve(Mᴴ)` a
   escalas mayores? (el micro-mecanismo no se conoce).
5. ¿La **pseudo-hermiticidad PROYECTADA** del operador basta para que `n_L(λ)` no sea un artefacto?
6. ¿Puede el estimador de σ_min **fallar** en regímenes no ensayados (lo contrario de lo que pasó: ¿y si falla
   dando σ_min grande cuando es pequeño)? La validación SVD fue solo X0∈{20,24,28}, out_sel=0.
7. ¿La reducción λ=1 → `build_MBdG` y la **regresión al benchmark FP** se ha ejecutado realmente para la
   realización Padé del solver (no solo a nivel operador)? (Gate P3.3 figura ⬜ no ejecutado en GATES_REGISTRY).
8. ¿Es correcto el **benchmark FP** mismo (convención de FP, κ̄ simétrica d=1.58)?

---

## 12. INVENTARIO DE ARCHIVOS RELEVANTES

**v1.0 congelado (física):** `mode_matching.py` (solver), `bdg_smatrix.py`, `bdg_channels.py`, `waterfall.py`,
`gpe.py`, `kappa_derive.py`, `FP_benchmark_table.py`, `FP_check_reproduction.py`, `FP_step{1..4}_*.py`;
docs `FP_BENCHMARK.md`, `README_FP_BENCHMARK.md`, `FP_SCALES/STEP2/STEP3/STEP4.md`, `REPRODUCE.md`;
datos `out/waterfall_state.npz`, `reference/fp_benchmark_reference.csv`; integridad `PROVENANCE.txt`,
`FP_REVIEW_MANIFEST.txt`.

**RUMBO#2 operador:** `residual_anatomy/NONLOCAL_SPEC.md`, `PADE_SPEC.md`, `OPERATOR_VALIDATION_CLOSURE.md`,
`STABILIZED_SOLVER_RESERVATION.md`, `README.md`, `P3_2*_DRAFT.md`, `P3_4_SMATRIX_DESIGN_DRAFT.md`;
gates `gates/p2_analytic.py`, `gate_A3_kernel.py`, `gate_P3_2A/B1/B2/B2c/B3*.py`, `gate_P3_4*.py`;
registro `gates/GATES_REGISTRY.md`.

**Bloque 2 solver (esta sesión):** diseños `STABILIZED_SOLVER_BLOCK1.md`, `BLOCK2_2.md`, `BLOCK2_4.md`,
`BLOCK2_LOG.md`; gates `gates/gate_B2_1_algebraic.py`, `gate_B2_2_conditioning.py`,
`gate_B2_3_convergence.py`, `gate_B2_4_X0_convergence.py`, `gate_B2_4b_phase_sensitivity.py`,
`gate_B2_4c_svd_verify.py`, `gate_B2_4d_adjoint_check.py`.

---

## 13. SÍNTESIS PARA EL AUDITOR

- **Lo firme:** v1.0 reproduce FP ≤0.75% (ancla); el operador Padé está validado *dentro de alcance
  documentado* (pseudo-herm proyectada); B2.1–B2.3 del solver son sólidos (ensamblado, condicionamiento
  benigno, **convergencia O(dx²) demostrada por Richardson**); el condicionamiento de la matriz **no explota**
  con X0 en el rango probado (SVD-confirmado); los moduli de S convergen en X0.
- **Lo frágil / abierto:** **C1 no resuelta** (objetivo central del solver); convergencia COMPLEJA de S en X0
  no demostrada (residuo de fase hipotético); **bug de solver adjunto** en código auditado que contamina las
  métricas de sensibilidad; **resultados del Bloque 2 no persistidos**; `n_L(λ)` **no medida**; regresión λ=1
  del solver Padé (P3.3) **no ejecutada**; cotejo FP-PDF pendiente.
- **Recomendación de la Fase 1:** NO declarar el solver validado ni la física resuelta. Antes de medir
  `n_L(λ)`: (a) migrar la ruta adjunta a `lu_solve(trans=2)` y revalidar sensibilidad; (b) cerrar la
  convergencia COMPLEJA de S en X0 (fase k-discreta) o declararla limitación; (c) ejecutar la regresión λ=1
  del solver Padé vs v1.0 (P3.3); (d) persistir resultados; (e) resolver C1 (pseudo-unitaridad convergente) con
  métricas IN/OUT separadas, que es el criterio a priori de la reserva.

*Fin del expediente Fase 1. No declara validado el solver ni la conclusión física. Preparado read-only,
2026-11-16.*
