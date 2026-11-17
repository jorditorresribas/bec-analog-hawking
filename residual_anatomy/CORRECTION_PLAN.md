# PLAN DE CORRECCIONES — solver de scattering para n_L(λ) con G fijo

**Estado:** plan de trabajo APROBADO (tras auditoría adversarial de 3 fases: expediente Fase 1 +
Opus Fase 2 + ChatGPT Fase 3, 2026-11-16). **La aprobación es del ORDEN DE TRABAJO, NO una validación
anticipada de los resultados que produzca.** No modifica `AUDIT_DOSSIER.md`, los gates, v1.0 ni el solver
hasta que cada paso pase su propio gate.

**Hallazgo estructural compartido por las tres auditorías:** el Bloque 2 (B2.1–B2.4) validó un solver del
operador **LOCAL λ=1** (`A_at`), **no el operador Padé K_P**. El instrumento para la pregunta física
**no existe todavía**. Buen condicionamiento del BVP local **NO** implica buen condicionamiento del BVP
Padé (cuyo `(2G+D)³` empeora con N, `PADE_SPEC §4.3`).

**Regla de evidencia:** nada se reutiliza como evidencia si se apoya en una ruta desacreditada
(p. ej. κ_out de B2.2/B2.4 vía `la.solve(Mᴴ)`). Cada gate persiste su resultado (Control B).

---

## SECUENCIA (con dependencias; P4 dividido: algebraica ANTES de P3)

```
P1 ─┐
    ├─► P4a ─► P3 ─► P4nf ─► P5 ─► P6
P2 ─┘        (Control A: Padé vs Weyl, antes de interpretar n_L)
Control B (manifiesto reproducible) = transversal a TODOS
```
1. **P1** Certificar la ruta adjunta.
2. **P2** Corregir la referencia de fase (en paralelo con P1).
3. **P4a** Fijar la convección skew-exacta del Padé (decisión algebraica ANTES de ensamblar).
4. **P3** Construir el solver Padé global (con la convección de P4a desde el inicio).
5. **P4nf** Validar pseudo-hermiticidad en el sector de scattering (numérica+física).
6. **P5** Regresión λ=1 (P3.3) + cierre de C1.
7. **P6** Extensión λ≠1 + medición de n_L(λ).

---

## P1 — Certificar la ruta adjunta  ✅ CERRADO dentro de alcance (2026-11-16)
**Cierre:** `gates/gate_P1_adjoint_cert.py`; artefactos `out/p1_adjoint_cert.csv`, `out/p1_selectors_X012.csv`
+ manifiestos. R1(LU trans)≡R2(QR)≡R3(GMRES precond. adjunto) a máquina en los 3 selectores (X0=12) y
X0=12/20/28 (out_sel=0); identidad `s=y†b`~1e−16; cotas SVD; sin mutación. Ruta defectuosa aislada a
**factorizar el M^H formado** (`r_fact` lo delata), micro-causa OpenBLAS NO identificada. κ_out antiguos
DESCARTADOS. **Límite:** X0=12/ω=0.10/matriz concreta; no universal. Detalle en `BLOCK2_LOG.md`.

**Problema (3 auditorías):** `la.solve(Mᴴ)` da adjunto incorrecto a N grande sin mecanismo conocido;
`lu_solve(trans=2)` y el estimador comparten rutina → sin 3ª vía independiente; κ_out de B2.2 contaminado.

- **Algebraica:** identidad `s=c†x` vs `s_adj=y†b` exacta; misma `M`,`b`,`c`,dtype; sin modificación in-situ.
- **Numérica — TRES rutas independientes** sobre la MISMA matriz/RHS:
  (i) `lu_solve(trans=2)`; (ii) **factorización independiente de Mᴴ construida explícitamente**;
  (iii) **GMRES** sobre `M†y=c` con **precondicionador documentado**.
  Registrar por ruta: `r_adj=‖M†y−c‖/‖c‖`, `r_primal=‖Mx−b‖/‖b‖` (**abs y rel**), y `s` vs `s_adj`.
  **GMRES se evalúa con el residuo del sistema ORIGINAL, no solo el precondicionado.** Investigar el
  **micro-mecanismo** del fallo de `la.solve` (F-contigüidad descartada; probar camino LAPACK, presión de
  memoria, dtype/overwrite). Barrido en N hasta donde la SVD es verdad-terreno.
- **Física:** ninguna. **Descartar** las sensibilidades antiguas; NO re-etiquetarlas "moderadas" sin tabla nueva.
- **Matices (ChatGPT):** una coincidencia entre dos soluciones NO basta si comparten un fallo → los
  residuos originales y la identidad primal-adjunta son los controles esenciales. **No exigir residuo
  máquina universal a GMRES**; fijar la tolerancia tras caracterizar las rutas directas. Para N mayores que
  la SVD no es "certificación": es **evidencia de estabilidad/consistencia**, no verdad-terreno independiente.
- **Gate:** residuo adjunto ~máquina en ≥2 rutas independientes **que coinciden entre sí**, en el rango de N
  accesible; más allá, evidencia de consistencia documentada como tal.

## P2 — Corregir la referencia de fase (des-rotado con propagador de caja)
**Problema:** el des-rotado usó `k` continuo; la caja propaga con `k_box=(2/dx)·arctan(k·dx/2)`
(derivación confirmada por las 3 auditorías contra `gate_B2_1`: `T_box=(I−dx/2A)⁻¹(I+dx/2A)`, autovalor `ik`).

- **Algebraica:** derivar el des-rotado **completo por canal y extremo**, con signos de entrada/salida
  explícitos; confirmar `D_out`,`D_in`. `k_box` es el autovalor homogéneo; **no sustituir k→k_box a ciegas**.
- **Numérica — distinguir TRES convergencias (ChatGPT):** (1) corrección algebraica de la fase;
  (2) convergencia en X0 a dx fijo; (3) convergencia en dx de la propia corrección. El gate **no** se limita
  a que desaparezca la tendencia lineal en X0: debe comprobar que los **nueve elementos complejos de
  `S_fis` convergen, con sus fases, al refinar dx**.
- **Física:** ninguna (consistencia numérica).
- **Gate:** `S_fis` compleja (9 elementos) converge en dx y en X0 con fase coherente; módulos siguen convergiendo.

## P4a — Convección skew-exacta del Padé (DECISIÓN algebraica, antes de P3)
**Problema:** `gate_B2_1`/`gate_P3_2B3` usan `-i(v·D1+½v')` (skew O(dx)); existe `½(v·D1+D1·v)` skew-exacta
en `gate0_pseudoherm.py` (ruta Weyl) **no trasladada** al Padé.

- **Algebraica:** incorporar `C_s=−i·½(v·D1+D1·v)` al `M_eff` Padé; verificar `C_s†=C_s` exacto
  (`D1†=−D1`, `v` real diag). **Escoger esta forma ANTES de ensamblar el solver** (evita implementar el
  operador dos veces). Verificar también el **tratamiento de bordes**.
- **Matiz (ChatGPT):** la skew-exacta elimina el defecto de **ese bloque**; NO demuestra pseudo-herm del
  `M` completo (con K_P, BC, métrica BdG), ni que conserve el continuo al orden esperado. Eso es P4nf.
- **Gate:** bloque convectivo skew-exacto a máquina en el ensamblado Padé.

## P3 — Solver Padé global (el instrumento que falta)
**Problema (crítico):** `gate_B2_1` usa `A_at` (local). El Padé `eps_P2 = s/2 + (λ−1)s²/(4(2G+s)) −
(λ−1)²s⁴/(16(2G+s)³)` (2º orden en δ, cubre k_* para λ∈[0.5,2]) **no está en el solver**.

- **Algebraica:** sistema aumentado disperso con auxiliares `p,q₁,q₂,q₃` y `(2G+D)` por band-solve, con
  **BC auxiliares explícitas** (la spec exige 2 por extremo); usar la convección de **P4a** desde el inicio;
  equivalencia aug≡denso.
- **Control de degeneración λ=1 (ChatGPT, nuevo):** en λ=1 los coeficientes `c₁,c₂` se anulan. Si el
  sistema aumentado conserva `p,q₁,q₂,q₃` sin contribución, **demostrar** que sus ecuaciones no dejan
  variables libres ni introducen singularidades: (i) rango del sistema aumentado en λ=1; (ii) límite λ→1
  continuo; (iii) ausencia de modos nulos/espurios de los auxiliares en ese límite.
- **Numérica:** repetir **B2.2–B2.4 sobre el Padé a λ=2**: ¿`cond` benigno o crece con N por `(2G+D)³`?
  → **secuencia de mallas (dato)**, no la advertencia de la spec. SVD verdad-terreno a N moderado.
  **Gate de espurios sobre el espectro del SISTEMA AUMENTADO** (no solo dispersión efectiva vs ideal):
  comprobar que las ramas auxiliares no introducen soluciones que satisfagan artificialmente las BC.
- **Física:** ninguna aún.
- **Gate:** aug Padé cuadrado, rango pleno, límite λ=1 regular, sin espurios auxiliares, cond caracterizado en N.

## P4nf — Pseudo-hermiticidad en el sector de scattering (numérica + física)
**Problema:** pseudo-herm de B3 solo **proyectada a baja k** (gaussianas suaves), posiblemente disjunta de
los canales físicos (k finito hasta k_*≈1.6).

- **Numérica — guardar TRES cantidades (ChatGPT):** (a) **defecto global** `M†η−ηM`; (b) **defecto
  proyectado** sobre los subespacios de scattering **entrante y saliente**; (c) **acoplamiento** del defecto
  entre el sector físico y su complemento. Primero **medir el contenido espectral** de los canales (qué k,
  qué amplitud a alta k). **La proyección usa la métrica BdG y la normalización modal apropiadas**, NO una
  base euclídea ortonormal. Barrido en N.
- **Física:** cuantificar si el defecto produce **violación medible de la norma en el observable n_L**.
  **NO asumir** que un defecto O(1) en rejilla invalida n_L si apenas acopla con los canales físicos
  (matiz ChatGPT: "disjunto"/"afecta necesariamente" es aún demasiado fuerte sin esta medición).
- **Gate:** defecto sobre el sector físico controlado y convergente, o declarado como límite con su magnitud
  y su efecto acotado sobre n_L.

## P5 — Regresión λ=1 (P3.3) + cierre de C1
**Problema:** P3.3 (⬜) y C1 (reservada) sin ejecutar; el baseline usa fondo interpolado (splines)+`np.gradient`
mientras v1.0 usa analítico cerrado → error adicional a controlar.

- **Algebraica:** confirmar que el solver Padé en λ=1 colapsa a `build_MBdG`.
- **Numérica:** **P3.3** reproducir β²=0.01008, T_up/(κ̄/2π)=1.005 (w=2) dentro de tolerancia, **controlando
  el error de interpolación del fondo** (o migrar a fondo analítico). **C1**: `S†η_out S=η_in` con métricas
  IN/OUT separadas, **convergente en dx y con meseta en X0** (criterios a priori de la reserva).
- **Física:** regresión al benchmark FP (ancla).
- **Gate:** P3.3 PASS vs v1.0 **y** C1 convergente. **Sin esto, NO se interpreta n_L.**

## P6 — Extensión λ≠1 y medición de n_L(λ)
**Problema:** `norm_current` y la base modal (`build_Xi` vía `A_homog`) son **específicas de λ=1**.

- **Algebraica:** base modal, `v_g` y corriente del Padé `J_∞=ε_λ'(k)(|U|²+|V|²)+v(|U|²−|V|²)` con
  `ε_λ'(k)=(c²k+(λ/2)k³)/(ε_λ+G)≠k` (NONLOCAL_SPEC ec.7–9); modos NO autovectores de `A_homog` local.
- **Numérica:** convergencia de n_L(λ) en dx, X0, N.
- **Física:** **medir n_L(λ), λ∈[0.5,2], G fijo → ¿persiste el colapso?** (pregunta central); comparar con el
  knob viejo (que reescala G).
- **CONDICIÓN DE ENTRADA ESTRICTA (ChatGPT):** NO basta con sistema cuadrado, rango pleno y estabilidad
  numérica. Debe existir evidencia de que (1) el operador implementado representa el **Padé previsto**
  (Control A), (2) bases modales/corrientes son las **adecuadas para λ≠1**, (3) la conservación de norma se
  comporta bien **en el sector físico** (P4nf) y C1 está cerrada (P5).
- **Gate:** n_L(λ) estable y reproducible con C1 satisfecha, antes de afirmar nada sobre el colapso.

---

## CONTROL TRANSVERSAL A — Padé vs Weyl (antes de interpretar n_L)
Superar los gates del solver NO demuestra que el Padé represente con precisión el operador no local de Weyl
que se quiere estudiar. Comparar: (i) dispersión homogénea en todo el k relevante; (ii) error de `eps_P2`
frente a Weyl por λ; (iii) convergencia al refinar malla; (iv) **comparar S y n_L, no solo los símbolos de
dispersión**. (`gate0_pseudoherm.py` tiene la referencia Weyl densa.)

## CONTROL TRANSVERSAL B — Manifiesto reproducible por ejecución
Cada gate importante produce un manifiesto (NO dejar resultados solo en stdout de terminal):

| Campo | Contenido |
|---|---|
| Código | Hash de cada archivo ejecutado + dependencias locales |
| Entorno | Python, NumPy, SciPy, BLAS/LAPACK, plataforma |
| Parámetros | λ, ω, X0, dx, N, G, perfil, condiciones de contorno |
| Resultado | CSV de datos COMPLETOS, no solo resúmenes |
| Diagnóstico | Residuos, rango, condicionamiento, tiempos |
| Reproducción | Comando exacto + código de salida |

---

## SEPARACIÓN DE TIPOS DE PRUEBA
- **Algebraica:** P1, P2, P4a, P3 (ensamblado/degeneración), P6 (bases/corrientes).
- **Numérica:** P1 (rutas/residuos), P2 (convergencia S compleja), P3 (cond/espurios), P4nf (defecto), P5 (C1), P6 (convergencia n_L), Control A.
- **Física:** P5 (regresión FP), P4nf (violación de norma en n_L), P6 (n_L(λ) / colapso), Control A (S/n_L Padé vs Weyl).

**NO se interpreta física hasta P5 (C1) cerrada y P6 cumple su condición de entrada estricta.**

*Plan congelado 2026-11-16 tras 3 auditorías convergentes. No valida resultados; fija el orden de trabajo.*
