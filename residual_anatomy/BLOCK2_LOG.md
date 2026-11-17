# BLOQUE 2 — Registro de gates del solver de scattering estabilizado (LOG VIVO)

Diseño en `STABILIZED_SOLVER_BLOCK1.md` (congelado). Este log es **vivo** (se amplía por gate).
Regla: no toca v1.0, ni los gates/docs del operador ya cerrados, ni `GATES_REGISTRY.md`.
Patrón: diseño en papel → revisión (ChatGPT) → implementar → auditar. "Validado" = dentro del
alcance documentado.

---

## B2.1 — Gate algebraico del BVP aumentado  ✅ CERRADO (2026-11-15)

**Archivo:** `gates/gate_B2_1_algebraic.py` (baseline λ=1 local, perfil CONGELADO, solo lectura).
**Revisión:** ChatGPT (estática, iterada) + verificación numérica propia. **Cerrado dentro de alcance.**

### Qué valida (estructura algebraica, SIN física)
- Ensamblado del BVP aumentado: interior por **esquema de caja** + matching por amplitudes
  aumentadas `Y_b − Ξ_b a_b = 0` (sin invertir Ξ) + 4 restricciones.
- **Cuadratura** `(4N+12)×(4N+12)` y **conteo** de restricciones (invariante `#in+#grow = 4`).
- **Clasificación modal por ω** con las guardas: dirección `sign(v_g)`, norma `sign(Ω)`,
  evanescente **dependiente del extremo** (izq crece si Im k>0; der crece si Im k<0), orden
  determinista `(Re k, Im k)`.
- **Rango numérico pleno** de 𝓜_aug (ker trivial ⇒ unicidad) en los casos ensayados
  (N=20,40,80,120; ω=0.05 y 0.10).
- **Matching consistente** (`Y_b=Ξ_b a_b` ~1e−15) + **prescripción explícita** de amplitudes
  entrantes (`a_in=s`, ~1e−16).
- **Guard de estructura completa de canales** (no solo #in): izq {1 in,1 out,1 dec,1 grow},
  der {2 in,2 out,0 dec,0 grow}. Solo con esa estructura se ejecuta el solve 3×3.

### Controles ejecutados
| caso | estructura | resultado |
|---|---|---|
| ω=0.10, N=60 | izq [in,out,dec,grow] · der [in,in,out,out] | PASS (res 7.7e−15) |
| ω=0.05, N=60 | idem | PASS |
| N=20…120 (ω=0.10) | rango pleno 92/172/332/492 | ker trivial, inv=4 |
| evan. derechos ω=0.5/0.8/1.2 | der gana evan; Im k<0→grow(excluir), Im k>0→dec | clasificación correcta |
| ω=0.80 (fuera de ventana) | izq/der 1+1 in | **FAIL(guard)**, no ejecuta solve |

### Lo que B2.1 NO valida (reservas explícitas)
1. Convergencia física del BVP, matriz **S**, **pseudo-unitaridad** `S†η_out S=η_in`.
2. **Ausencia general de modos espurios** (solo unicidad en los casos ensayados).
3. **Condicionamiento al crecer N:** rango pleno a N finito NO lo garantiza. Observado
   `σ_min/σ_max` DECRECIENTE con N (2.0e−3 @N=20 → 4.0e−4 @N=120) — **condicionamiento
   dependiente de escala, no pérdida de rango** → equilibrado/escalado **pendiente**.
4. λ≠1 (requiere ruta aumentada p,q1,q2,q3 con BC propias).
5. Continuidad de ramas modales al variar ω (orden determinista es por-ω, no seguimiento continuo).

---

## B2.2 — Estudio de condicionamiento de 𝓜_aug  ◐ DISEÑO CERRADO EN PAPEL (2026-11-15)

**Documento:** `STABILIZED_SOLVER_BLOCK2_2.md` (congelado, diseño; implementación NO autorizada).
**Revisión:** ChatGPT (iterada). Separación de ejes N (este gate) / X0 (gate aparte).

Decisiones congeladas: ℓ=ξ_u=1/√2 (healing upstream, indep. N); 2-norma; un adjunto por amplitud de
salida (3×3); selector bajo equilibrado `D_c†c`; crudo→Ruiz→físico; métricas `κ_{out,b}`, `κ_{out,M}`,
`η_back`, `Δ_S`; tolerancias PENDIENTES justificadas. Criterios C1–C5 fijados.

**Diseño + tolerancias:** cerrados en papel (constantes provisionales).

### Calibración ◑ CERRADA (instrumentación auditada, 2026-11-15)
**Gate:** `gates/gate_B2_2_conditioning.py` (importa el ensamblado auditado de B2.1; NO lo modifica).
**Auditoría ChatGPT:** 2 rondas (algebraica + instrumentación) → aprobada. Ejecución propia.

**Resultados (ω=0.10 y 0.05, X0=12, barrido N; ventana 3×3):**
- **Estabilidad algebraica excelente y N-estable:** `η_back~5e−17…1.5e−16` (máquina); control adjunto
  `δ_adj~1e−14` abs / `1e−13` rel; **acuerdo de las 3 representaciones (cruda/Ruiz/física) ~1e−14**;
  rango **RESUELTO** siempre (σ_min~1e−3 ≫ piso~5e−13).
- **Condicionamiento NO es artefacto de escala:** `drop` Ruiz ~2×, física 0.82 (no ayuda) →
  **no eliminado por los escalados ensayados**; `cond≈C·N^p` con **p≈0.97** (no satura). Las
  sensibilidades y errores observados indican que es **benigno para las amplitudes en los casos
  estudiados** (conclusión empírica, NO universal).
- **Sensibilidades (absolutas ⟂ relativas, no intercambiables):** relativas `κ_out,b~10²`,
  `κ_out,M~10⁴` (crecen ~N); **absolutas** `‖y‖~20` (a δb), `A_out,M~2e3` (a δ𝓜, cota 1er orden).
  Error forward de amplitudes `~1e−11` = **ESTIMACIÓN de 1er orden** (`κ_out,M·η_back`), **no** error
  forward medido directamente → amplitudes estables pese a cond~10³ en los casos estudiados.
- **Convergencia de S = limitada por DISCRETIZACIÓN:** orden Richardson (ratio dx=2) `p_obs→2.00`
  (1.97/1.99/1.99 en ω=0.10; íd. ω=0.05) → **compatible con O(dx²)**. Pero `Δ_S` es **O(1) en
  N=40–160**; ~1% a dx~0.05 (N~480), ~0.35% a dx~0.03 (N~800), aún bajando. El error algebraico
  (~1e−14) está **~10 órdenes por debajo** del de discretización.

**Constantes provisionales calibradas:** η_back ≤ ~1.5e−16 (umbral holgado); acuerdo repr ≤ ~3e−14;
rango RESUELTO (C₄~1, banda 10³); cond~N¹ (p≈0.97, C~25); κ_out,M~10⁴ creciendo ~N; **plató de Δ_S =
discretización** (a fijar en régimen fino N≳480).

**DISTINCIÓN (obligatoria):** calibración cerrada **≠** convergencia física de S **≠** cierre global de
B2.2. La parte algebraica/estabilidad está caracterizada; S **no** está convergido a precisión física.

**Próximo (NO iniciado, sin autorización):** barrido fino con la **caja como referencia de 2º orden**
(no cambiar esquema y resolución a la vez) antes de plantear un esquema de orden superior. No abrir el
siguiente gate sin luz verde.

---

## B2.3 — Barrido fino de convergencia de S  ✅ CERRADO (caracterización dx, 2026-11-15)

**Gate:** `gates/gate_B2_3_convergence.py` (reutiliza por import B2.1 `assemble`/`rhs_for` y B2.2
`physical_scaling`/`smatrix_and_solves`; NO los modifica). **Auditoría ChatGPT:** 2 rondas + control
N≈2000 → aprobada. **Esquema de caja de 2º orden, X0=12 FIJO** (sin cambiar esquema y resolución a la vez).

**Resultados (ω=0.10 y 0.05; S representación física):**
1. **Esquema caja 2º orden, X0=12 fijo.**
2. **Orden observado compatible con O(dx²):** Richardson (ratio dx=2) `p_Frob` 1.97→1.99→**1.999**
   (triple 500/1000/2000), `p_elem∈[1.97,2.00]`. Incluye control **N=2000**.
3. **Δ_S(2000)=1.01e−3** (dx=0.012), medido respecto a la referencia Richardson; `1e−2` alcanzado en
   N=800; `1e−3` ≈ N~2000 (medido, concordante con extrapolación N~2008); `1e−4` extrapolado y marcado
   **NO respaldado** (bajo la autoconsistencia de S_ref). Ratio Δ_S(1200)/Δ_S(2000)=2.77 ≈ (2000/1200)²=2.78.
4. **Autoconsistencia de S_ref = 2.7e−5** (pares finos 1000/2000 vs 600/1200; mejora desde 1.2e−4) —
   es **autoconsistencia entre extrapolaciones Richardson, NO cota de error al continuo**.
5. **Pendientes:** convergencia en **X0** (gate aparte, reservado) y **cierre global de B2.2**.

**NO demostrado:** que S_ref coincida con el límite continuo a <2.7e−5; que S_box(X0=12) = S física;
convergencia en X0; pseudo-unitaridad física; validez de n_L(λ).

**Próximo (NO iniciado, sin autorización):** gate de **ventana/convergencia en X0** (eje independiente).
El solve disperso queda como **optimización posterior** (a validar por equivalencia vs denso), no ahora.

---

## FASE DE CORRECCIÓN (tras auditoría adversarial de 3 fases, 2026-11-16)

Auditoría (`AUDIT_DOSSIER.md` Fase 1 + Opus Fase 2 + ChatGPT Fase 3, convergentes): **el Bloque 2 validó el
solver del operador LOCAL λ=1 (`A_at`), NO el operador Padé** — el instrumento para la física n_L(λ) no
existe. Las alarmas de B2.4 quedaron explicadas: Δ_S O(1) = artefacto de fase; "cond explota" = **bug del
solver adjunto** (`la.solve`/formar M^H). Plan en `CORRECTION_PLAN.md` (P1→P2→P4a→P3→P4nf→P5→P6).

### P1 — Certificar la ruta adjunta  ✅ CERRADO dentro de alcance (2026-11-16)
**Gate:** `gates/gate_P1_adjoint_cert.py`. **Artefactos persistidos (Control B):** `out/p1_adjoint_cert.csv`
(barrido X0=12/20/28, out_sel=0), `out/p1_selectors_X012.csv` (out_sel 0/1/2) + manifiestos (sha256, npz,
entorno, OpenBLAS 0.3.27).
- **R1 (LU trans-solve) ≡ R2 (QR) ≡ R3 (GMRES, precond. ILU adjunto)** a residuo máquina (1e−14–1e−13) en
  los 3 selectores (X0=12) y en X0=12/20/28 (out_sel=0); identidad `s=y†b` ~1e−16; cotas SVD cumplidas;
  entradas sin mutar.
- **Ruta defectuosa** (formar M^H explícitamente): residuo crece con N (1.9e−5 → 7.5 → 5.4e6); `r_fact`
  muestra que la degradación **ya aparece en la FACTORIZACIÓN** del M^H formado (idéntica en
  `la.solve`/LU-explícita/C/F → la contigüidad no la explica). **No se afirma** que esté descartada
  cualquier otra característica de `la.solve`; lo demostrado es que, en los controles realizados, el fallo
  ya está presente en la factorización. Micro-causa interna (OpenBLAS): **NO identificada** (legítimo).
- → **toda la columna κ_out/δ_adj de B2.2/B2.4 es ESPURIA** (ruta defectuosa); se descarta, no se recupera.
**Límite del cierre:** certificación en X0=12 (por selector) / 12,20,28 (out_sel=0), ω=0.10, la matriz
concreta del gate; **NO** prueba estabilidad universal para cualquier ω/malla/dominio.
**Siguiente autorizado:** P2 (re-referencia de fase). **NO recalcular κ_out** hasta resolver/verificar la
fase de S; cuando se autorice, será **serie nueva** (los κ_out antiguos no se recuperan).

---
*B2.1/B2.2/B2.3 cerrados (operador LOCAL). Auditoría 3 fases: el solver Padé NO existe; C1 y n_L abiertos.
P1 (ruta adjunta) CERRADO dentro de alcance. Siguiente: P2. Nada más abierto sin autorización.*
