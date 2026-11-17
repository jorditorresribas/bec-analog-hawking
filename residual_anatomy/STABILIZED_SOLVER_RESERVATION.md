# ACTA DE RESERVA — Solver de scattering ESTABILIZADO (fase posterior)

**Fecha:** 2026-11-15. **Estado: RESERVADO, NO diseñado, NO implementado.** Este documento fija el
OBJETIVO y los CRITERIOS DE ACEPTACIÓN a priori. **NO contiene arquitectura** (se elige con datos, en su
fase); **NO fija tolerancias numéricas nuevas** (se establecerán antes de implementar, justificadas con la
precisión y los controles disponibles). No toca v1.0 ni el operador ya validado.

---

## 1. Motivación

La validación del operador Padé con G fijo (`OPERATOR_VALIDATION_CLOSURE.md`) quedó cerrada **salvo la
conservación de scattering convergente (C1)**. El baseline λ=1 por transferencia (de v1.0) mostró:
- **floor** ~5e−5 en `max|S†η_out S − η_in|` sin convergencia demostrada en dx;
- **sin meseta en X0** (ventana de estabilidad + blow-up del modo evanescente creciente ~e^{+κX0}).

Los controles A/B (réplica exacta de `B.smatrix`; invariancia a la escala evanescente) y la sensibilidad
al umbral modal **acotan el origen a la PROPAGACIÓN por transferencia**, no al operador. Por tanto, para
certificar C1 con rigor se necesita un **solver de scattering estabilizado** — tarea de nivel instrumento,
separada de la validación del operador.

**Doble valor (por eso NO se abandona):**
1. **Broche de oro de la fase numérica:** certifica que el operador produce una matriz S fiable y conserva
   la norma de Krein esperada.
2. **Instrumento imprescindible de la fase física:** es lo que permitirá medir `n_L(λ)` con G fijo (la
   pregunta central aún ABIERTA), extendiéndose al operador Padé para λ≠1.

---

## 2. Objetivo

Construir un solver de scattering **numéricamente estable** que produzca la matriz S del problema BdG
(primero baseline λ=1 local, luego operador Padé con G fijo), sin la inestabilidad del transfer ingenuo,
y que **supere los criterios de aceptación de abajo**. La arquitectura (recursión de matriz S,
propagación de subespacios, matching directo, Riccati u otra) **se decide en su fase, no aquí**.

---

## 3. Criterios de aceptación (a priori, NO se moverán a posteriori)

| Criterio | Requisito |
|---|---|
| **Convergencia espacial** | S estable al refinar `dx` (pendiente/meseta registrada, no un floor) |
| **Estabilidad asintótica** | **meseta** de S respecto a `X0` (sin blow-up; el modo creciente controlado) |
| **Conservación** | `S† η_out S = η_in` (pseudo-unitaridad con métricas IN/OUT SEPARADAS) |
| **Métricas** | η_in y η_out correctamente separadas y asignadas por canal |
| **Modos evanescentes** | control del modo creciente y del **condicionamiento** (sin dependencia de su escala) |
| **Espectro** | ausencia de **soluciones espurias** (auxiliares/polos no físicos clasificados como canal) |
| **Regresión** | coincidencia con v1.0 (`B.smatrix`) **en la ventana donde la transferencia es estable** |

**Tolerancias numéricas:** NO se fijan aquí. Se establecerán antes de implementar, justificándolas con la
precisión alcanzable y los controles de convergencia (igual que en cada gate previo).

---

## 4. Orden de desarrollo (cuando se autorice la fase)

1. **Diseño en papel** de la arquitectura candidata + controles (revisión antes de implementar).
2. **Baseline λ=1 (local):** implementar, comparar con v1.0 en su ventana estable, exigir los criterios de
   §3 (convergencia dx, meseta X0, pseudo-unitaridad). Prueba de infraestructura.
3. **Extensión al operador Padé (λ≠1):** sistema aumentado disperso (auxiliares p,q1,q2,q3) con sus
   **gates de BC propios** (derivar sistema → justificar BC → homogéneo reproduce raíces/modos → sin modos
   espurios → scattering), raíces del **sistema aumentado completo** (no de la dispersión ideal; posibles
   polos/modos auxiliares), `vg` de la rama espectral real.
4. **Entonces:** conservación de scattering para λ≠1 + medir `n_L(λ)`, G fijo, λ∈[0.5,2] → **¿persiste el
   colapso?** (pregunta física central).

---

## 5. Qué NO es este documento
- NO es un diseño ni una arquitectura.
- NO autoriza implementación.
- NO forma parte de la validación del operador (que se cierra sin él).
- NO fija umbrales numéricos nuevos.

**Siguiente acción sobre este documento:** abrir la fase solo cuando el usuario lo decida; entonces,
diseño en papel + revisión antes de cualquier código.
