# Cierre de la VALIDACIÓN DEL OPERADOR Padé con G fijo (RUMBO #2)

**Fecha:** 2026-11-15. **Alcance:** validación del OPERADOR (no del solver de scattering ni de la física).
**Regla:** v1.0 FP (benchmark congelado) y `bdg_smatrix.py` NO se modifican; se usan como referencia.
Revisión externa en paralelo (ChatGPT); cada gate auditado antes de cerrarse.

> **Qué es este documento:** el cierre de la fase "¿el operador Padé `K_P` con `G=g·n(x)` FIJO está bien
> construido y tiene las propiedades que debe (simetría, Green, pseudo-hermiticidad, conservación a nivel
> de operador)?". **NO** cierra la conservación de scattering convergente (C1), que queda **reservada**
> para la fase del solver estabilizado (ver `STABILIZED_SOLVER_RESERVATION.md`).

Operador: `L_± = K_P ± C + (G+Qp)`, `M_eff = [[L_+, G],[−G, −L_−]]`,
`K_P = ½D + c₁ DAD + c₂ D²A³D²`, `A=(2·diag(G)+D)⁻¹`, `D=−∂²`, `c₁=(λ−1)/4`, `c₂=−(λ−1)²/16`.
Realización: **Dirichlet por nodos fantasma** (plana, euclíd-simétrica).

---

## 1. Resultados por gate (todos CERRADOS, con alcance explícito)

| Gate | Qué valida | Resultado | Alcance / limitación |
|---|---|---|---|
| **Paso 2 analítico** (P2.1–2.7) | dominios, Green, BC, residuo G', flujo δ_J, positividad 2G+D | ✅ bajo hipótesis declaradas | NO maximalidad del dominio; NO autoadjunción del operador inhomogéneo |
| **P3.1** A³ vs Green | band-solve A³ = símbolo 1/(k²+κ²)³ en el interior | ✅ interior O(dx²) (pend 2.00) | Robin es aprox. O(1/L); caja O(L²e^{−κL}) |
| **P3.2-A** Schur | sistema aumentado disperso ≡ K_P | ✅ aug≡seq≡denso ~1e−12 ∀N,L,BC | SOLO equivalencia algebraica |
| **P3.2-B1** simetría | K_P euclíd-autoadjunto (discretización plana) | ✅ 4e−13 | Robin/Dirichlet rompen la simetría; C herm. solo O(dx) (rel) |
| **P3.2-B2 + B2c** Green discreta + continuo | identidad de Green discreta; convergencia al concomitante continuo | ✅ (opción B) B2.3 cerrado por B2c: BF_disc→BF_cont O(dx²); F_matriz→BF_cont O(dx) | B2.2 = consistencia numérica, **NO** identidad SBP formal; G=1, región [−6,6] |
| **Gate 1** homogéneo | dispersión Ω²=c²k²+λk⁴/4, modos, corrientes, λ=1 exacto | ✅ máquina + O(dx²) | fondo homogéneo |
| **Gate 2a** regresión operador λ=1 | operador no local(λ=1) → build_MBdG | ✅ O(dx²) | nivel operador (NO scattering) |
| **P3.2-B3** pseudo-herm M_eff | `R=M†𝒢−𝒢M` proyectado sobre subespacio suave | ✅ O(dx²) (pend 2.00), proj_C domina, proj_KP~máquina | **PROYECTADA, NO exacta**; defecto matricial O(1) persiste en modos de rejilla; realización densa se condiciona a N grande |
| **P3.4-C2** balance dinámico regional | `dQ_R/dt = i⟨Ψ,[M,P_R]Ψ⟩_𝒢 + iΨ†(M†𝒢−𝒢M)P_R Ψ` durante evolución `expm(−iMt)` | ✅ balance cierra; O(dt²) temporal; ΔQ_R/flujo convergen con N; ∫defecto→0 O(dx²); recon~máquina | **estado inicial de PRUEBA** (no modo BdG propagante); identidad algebraica + balance discreto, NO física de Hawking; realización λ=1 local |

**Nota clave:** B3 y P3.4-C2 son evidencia de conservación **a nivel de operador** e **independientes del
método de transferencia** (B3 opera sobre la matriz; C2 evoluciona M con `expm`). No heredan el floor de C1.

---

## 2. Lo que SÍ queda validado (dentro del alcance documentado)

- La **construcción** del operador Padé con G fijo y su **equivalencia aumentada** (P3.2-A).
- Su **simetría** euclídea en la discretización plana **y para la realización algebraica ensayada (B1),
  con las limitaciones de contorno documentadas** (Robin/Dirichlet la rompen; no es independiente de la
  realización), y la **identidad de Green discreta** con convergencia O(dx²) al concomitante continuo
  (B2+B2c).
- Su **pseudo-hermiticidad PROYECTADA** sobre el subespacio suave, O(dx²), dominada por el término
  convectivo C (B3).
- Su **conservación dinámica regional** (balance con término de defecto), con convergencia espacial y
  temporal controladas (P3.4-C2).
- La **regresión λ=1** a nivel de operador (= build_MBdG) y la **dispersión/modos** homogéneos (Gates 1, 2a).
- **Robustez de la clasificación modal** al umbral |Im k|<thr (invariante 1e−6…1e−12; S−Sref=0.0).

---

## 3. Lo que NO se afirma (limitaciones y reservas explícitas)

1. **NO** se demuestra pseudo-hermiticidad **matricial exacta** ni autoadjunción en el **dominio continuo**
   (B3 es proyectada; el defecto matricial O(1) persiste en modos de rejilla).
2. **NO** se valida la **pseudo-unitaridad convergente de la matriz S (C1)**. El baseline λ=1 por
   transferencia tiene un **floor ~5e−5** (sin convergencia demostrada en dx) y **sin meseta en X0**
   (ventana de estabilidad + blow-up del modo evanescente creciente). Los controles A/B (réplica exacta de
   `B.smatrix` e invariancia a la escala evanescente) y la sensibilidad al umbral **descartan varias
   fuentes posibles de error** (réplica del matching, normalización de los evanescentes, clasificación
   modal). Los resultados **apuntan a la propagación por transferencia como principal sospechosa** del
   floor y del blow-up, **pero no demuestran que sea su causa exclusiva**. **C1 queda RESERVADO al solver
   estabilizado** (`STABILIZED_SOLVER_RESERVATION.md`).
3. **NO** se extiende el resultado a la realización **Robin**, a **otras condiciones de contorno**, ni al
   operador **completo fuera de la realización Dirichlet-fantasma ensayada**.
4. **NO** se afirma identidad **SBP formal** (B2.2 es consistencia numérica).
5. **NO** se usa un **paquete BdG propagante** en C2 (estado de prueba); el paquete físico queda para la
   fase de scattering.
6. La realización **densa** `D²A³D²` se condiciona a N grande (afecta a normas crudas, no al subespacio
   suave ni al balance regional); la realización en **banda** es tarea de la fase posterior.
7. **NO** se afirma **ausencia universal de modos espurios** en cualquier realización (la clasificación de
   raíces físicas vs auxiliares/polos es un control del **solver de scattering**, fase posterior).
8. **NO** se responde todavía la pregunta física `n_L(λ)` (fase posterior; depende del solver Padé).

---

## 4. Veredicto

**El operador Padé con G fijo queda VALIDADO dentro del alcance documentado:** construcción, simetría,
identidad de Green (→ continuo O(dx²)), pseudo-hermiticidad proyectada (O(dx²)) y conservación dinámica
regional, todo a nivel de operador e independiente del método de transferencia, más las anclas λ=1.

**"Validado dentro del alcance documentado" NO significa** (y nunca debe citarse sin esta coletilla):
(a) pseudo-hermiticidad **matricial exacta**; (b) autoadjunción del operador **continuo** para condiciones
de contorno generales; (c) **conservación de scattering convergente** (C1); (d) validez física de
**n_L(λ)**; (e) **ausencia universal de modos espurios** en cualquier realización. Lo demostrado es
pseudo-hermiticidad **proyectada** sobre el subespacio suave y una **identidad dinámica discreta** (con el
término de defecto incluido) — resultados relevantes, pero más acotados que una demostración general de
consistencia física del modelo.

**Reserva formal:** la **conservación de scattering convergente (C1)** NO está validada y se traslada a la
fase del **solver estabilizado**. El cierre de la validación del operador **no la presupone**.

Próxima fase (reservada, NO iniciada): `STABILIZED_SOLVER_RESERVATION.md`.
