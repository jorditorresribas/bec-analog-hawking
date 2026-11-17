# P3.2 — BdG aumentado disperso y pseudo-hermiticidad (BORRADOR, NO CONGELAR)

**Estado:** borrador matemático en papel. NO implementa solver. Marca explícita de lo PROBADO vs
HIPÓTESIS. Cautela central (ChatGPT): la métrica no se elige para forzar `M†Σ=ΣM`; se deriva de las
variables auxiliares y la equivalencia con `K_P`, con términos de borde. NO toca v1.0.

## 1. Operador físico BdG
Estado físico `(u,v)`, métrica `η=diag(1,−1)`, norma `Q=∫(|u|²−|v|²)dx`. Con `D=−∂²`,
`A=(2G+D)⁻¹`, convectivo `C=−i(v∂_x+v'/2)` y el cinético Padé
`K_P = ½D + c₁ DAD + c₂ D²A³D²` (c₁=(λ−1)/4, c₂=−(λ−1)²/16):

    L_± = K_P ± C + (G+Qp) ,     M = [[ L_+ , G ],[ −G , −L_− ]] ,     M (u,v)ᵀ = ω (u,v)ᵀ   (1)

**`C` es formalmente HERMÍTICO** (no antihermítico): `A_v=v∂+v'/2` es antihermítico (`A_v†=−A_v`,
salvo borde), luego `C=−iA_v` cumple `C†=iA_v†=−iA_v=C`. [N verificado: `|C−C†|≈borde`, `|C+C†|=O(1)`.]
**Signo del convectivo (contrastado con `build_op_lambda`):** el símbolo de `L_+` es `tk²+vk+G`
(`_symbol_lam`: `Lp=tk²+vk+G`); como `C e^{ikx}=+vk·e^{ikx}` [N verificado: 0.2407 vs vk=0.2410],
el convectivo de `L_+` es `+C` → **`L_± = K_P ± C`** (`L_+=K_P+C`, `L_−=K_P−C`). Idéntico a
`build_op_lambda` salvo que el cinético local se sustituye por `K_P` (§NONLOCAL_SPEC 3).
Como `K_P`, `C` y `G+Qp` son formalmente hermíticos, `L_±` lo son (en su dominio, salvo borde), y
`M†η=ηM` formal (ver §6b, con la salvedad de los términos de borde).

## 2. Realización aumentada (dispersa)
`K_P` actúa sobre `u` (en `L_+`) y sobre `v` (en `L_−`) → auxiliares POR CAMPO:

    p_a = A D a ,  q1_a = A D² a ,  q2_a = A q1_a ,  q3_a = A q2_a ,   a ∈ {u,v}
    K_P a = ½ D a + c₁ D p_a + c₂ D² q3_a

Estado aumentado (10 campos):  Ψ = (u, v, p_u,q1_u,q2_u,q3_u, p_v,q1_v,q2_v,q3_v)ᵀ.

**Filas dinámicas (llevan ω):**

    ω u = ½Du + c₁ D p_u + c₂ D² q3_u + C u + (G+Qp)u + G v
    ω v = −G u − [ ½Dv + c₁ D p_v + c₂ D² q3_v − C v + (G+Qp)v ]

**Filas de restricción (algebraicas, SIN ω; una por auxiliar, para a∈{u,v}):**

    (2G+D) p_a − D a = 0 ;  (2G+D) q1_a − D² a = 0 ;  (2G+D) q2_a − q1_a = 0 ;  (2G+D) q3_a − q2_a = 0

Problema generalizado  𝕄 Ψ = ω 𝔼 Ψ,  con `𝔼 = diag(1,1,0,…,0)` (solo `u,v` dinámicos).
**Estructura dispersa (corrección ChatGPT):** las restricciones son dispersas de BANDA — `(2G+D)` es
tridiagonal, pero los términos `D²a` (en `q1_a` y en el `c₂ D²q3` dinámico) son **pentadiagonales**.
El sistema aumentado completo conserva una estructura dispersa de banda (ancho fijo), no toda tridiagonal.

## 3. Identidad algebraica con K_P  [PROBADO — es la construcción]
Eliminar los auxiliares = complemento de Schur del bloque auxiliar = aplicar `A`. Como `(2G+D)` es
invertible **en su dominio y BC concretos** (P2.2; la invertibilidad NO es independiente del dominio),
la eliminación es única y **exacta**: recupera (1) con el `K_P` definido por ESE mismo `A`. No hay
aproximación en la equivalencia; la única aproximación es Padé (ε_P≈ε_λ), ajena a esta reducción.

**TRES equivalencias distintas (NO confundir; la algebraica NO implica las otras dos):**
- (i) **algebraica** — Schur reproduce `K_P` como expresión. PROBADO aquí (§3).
- (ii) **de dominios** — las BC auxiliares reproducen el dominio del resolvente `A` que define `K_P`.
  Es cuestión de las BC (§4) + la prueba `A³` (P3.1), NO se sigue de (i).
- (iii) **pseudo-herm efectiva** — `M_eff` cumple la identidad con `η`, contando el término de borde.
  Es §6b (continuo) + gate discreto (§8), NO se sigue de (i) ni de (ii).
En scattering abierto el operador tiene identidad de Green CON flujo de borde: NO se afirma
`M_eff†η=ηM_eff` como identidad de operador en cualquier dominio sin la cancelación de esos términos
(solo se da en el bulk; el borde ES el flujo físico, §6b).

## 4. Condiciones de contorno  [heredadas de §9/P3.1, con matices]
- Auxiliares `p_a,q_i_a`: cada solve `(2G+D)·=·` con BC AUTOADJUNTA. **Matiz (ChatGPT):**
  (a) **Dirichlet** en caja finita define un problema autoadjunto para `(2G+D)`, pero **NO es
  exactamente el inverso de recta completa** — la equivalencia exige el límite de caja (convergencia en
  L, medida en P3.1). (b) **Robin real** representa la cola evanescente asintótica **solo cuando el
  perfil ya es suficientemente constante** en el extremo; NO es "transparente exacta" para un perfil
  inhomogéneo arbitrario. Anulan los términos de borde INTERNOS (autoadj. de `L,L³`) bajo esas salvedades.
- Campos físicos `u,v`: BC de scattering (radiación / matching modal en los extremos), NO homogéneas.

## 5. Modos espurios — enunciado ACOTADO  [parte probada, parte a comprobar por convergencia]
**Corrección de alcance (ChatGPT):** para `G(x)` INHOMOGÉNEO, las soluciones homogéneas de
`(2G(x)−∂²)φ=0` NO son exponenciales exactas; solo en las **regiones asintóticas** (`G(x)→G_±`)
aparecen colas `e^{±κ_± x}` con `κ_±=√(2G_±)`. La invertibilidad de `(2G+D)` por sí sola **NO
garantiza** ausencia de modos espurios ni independencia de la BC de caja.

Lo que SÍ se puede afirmar (y no más):
> En las regiones homogéneas asintóticas, el bloque auxiliar `(2G+D)` no introduce raíces
> PROPAGANTES reales adicionales (único polo en `k²=−2G_±`, imaginario). Si el resolvente se
> implementa con el MISMO dominio y BC que definen `K_P`, la eliminación de Schur recupera el problema
> físico. Los posibles efectos de las BC auxiliares de caja **deben comprobarse por convergencia**
> (en N, L, BC), no se dan por garantizados (residuo de borde `O(e^{−κL})` medido en P2.6).

## 6. Pseudo-hermiticidad — la métrica NO se impone  [parte PROBADA, parte ABIERTA]

### 6a. La métrica física es DEGENERADA en el espacio aumentado
La norma física `Q=∫(|u|²−|v|²)` NO involucra a los auxiliares (son esclavos, determinados por `u,v`).
Luego la métrica física sobre el espacio aumentado es **degenerada**:

    Σ_phys = diag(η, 0_aux) = diag(1,−1,0,…,0)

**No** es una métrica no degenerada para un enunciado estándar `𝕄†Σ=Σ𝕄`. Inventar una `Σ` no
degenerada para que la pseudo-herm "salga" sería exactamente el error que advierte ChatGPT.

### 6b. Ruta RIGUROSA (recomendada): η-pseudo-herm del operador EFECTIVO  [PROBADO en continuo]
La pseudo-hermiticidad es propiedad del operador EFECTIVO `M_eff` tras eliminar auxiliares (Schur) =
el `M` físico (1) con `K_P`. Como (demostrado en §9):
- `K_P` es SIMÉTRICO salvo el flujo de borde físico (términos externos = `J_∞`, no espurios);
- `C=−i(v∂+v'/2)` es formalmente **HERMÍTICO** (§1), salvo borde; `G+Qp` real.
⟹ `L_±` formalmente hermíticos ⟹ `η M_eff − M_eff† η = 0` en el BULK (identidad formal interior),
con término de borde = flujo físico de canales (cuya conservación = pseudo-unit. de S, P3.4).
**Tres enunciados distintos, NO intercambiables:** (1) identidad formal en el interior; (2) identidad
de Green COMPLETA con los términos de borde explícitos (escrita y verificada en `PADE_SPEC §9.3`:
`½𝒲(w,u)+c₁[𝒲(w,p_u)+𝒲(p_w,u)]+c₂[𝒲₂(w,q3_u)+𝒲₂(q3_w,u)]`, con flujo en §9.7); (3) pseudo-unit.
de S `S†ηS=η`. El gate P3.2-B (§8) mide el DEFECTO INTERIOR `max‖M_eff†η−ηM_eff‖/‖M_eff‖` en región
central — **diagnóstico del bulk**, que debe acompañarse de la identidad de Green (§9.3) y del flujo
(P3.4), no sustituirlos. NO requiere elegir ninguna `Σ`: usa la `η` física y el `M_eff` real.

### 6c. Ruta de `Σ` en el espacio completo  [ABIERTO — hipótesis, NO resuelto]
Si se quisiera una `Σ` NO degenerada sobre los 10 campos con `𝕄†Σ=Σ𝕄` (+término de borde), debería
DERIVARSE de que los `q` son imágenes de `A` (p.ej. pesos `A`-dependientes en los bloques auxiliares),
no postularse. No está demostrado que tal `Σ` exista sin circularidad. **Se deja como cuestión abierta
y NO se usa para el gate.** La ruta 6b la evita por completo.

## 7. Qué está PROBADO vs HIPÓTESIS (explícito)
| Resultado | Estado |
|---|---|
| Equivalencia algebraica aumentado ↔ `K_P` (Schur) | **PROBADO** (construcción; invertibilidad en su dominio) |
| Sin raíces PROPAGANTES espurias en regiones asintóticas | **PROBADO** (§5, acotado) |
| Equivalencia de DOMINIOS (BC auxiliares ≡ dominio de `A`) para `G(x)` | **PARCIAL**: P3.1 solo probó `G=1` homogéneo; inhomogéneo por convergencia (PENDIENTE) |
| BC auxiliares autoadjuntas en caja (internos se anulan) | **PROBADO** con salvedades §4 (Dirichlet≠recta completa exacta; Robin solo asintótico) |
| `K_P` simétrico en el BULK salvo flujo físico de borde | **PROBADO en continuo** (§9, P2.3) |
| η-pseudo-herm del `M_eff` físico (ruta 6b) | **bulk formal OK**; Green completo + flujo + discreto = gate P3.2-B/P3.4 (PENDIENTE) |
| `Σ` no degenerada en espacio completo (ruta 6c) | **ABIERTO** (no imponer) |
| Maximalidad/autoadj. del operador inhomogéneo | **ABIERTO** (§9.6) |

**Estructura de problema generalizado (precisión):** el sistema aumentado es `𝕄Ψ=ω𝔼Ψ` con `𝔼`
SINGULAR (parte esencial de la estructura). NO confundir la pseudo-herm del operador REDUCIDO `M_eff`
(con `η`) con una pseudo-herm de la matriz aumentada `𝕄` (que no la tiene con una métrica no degenerada;
ruta 6c abierta).

## 8. Diseño del gate P3.2 — DOS pruebas independientes (corrección ChatGPT)
La contradicción de la versión previa ("medir sobre el sistema con TODOS los auxiliares" vs "medir el
`M_eff` reducido con `η`") se resuelve **separando** en dos pruebas; NO se mide pseudo-herm sobre la
matriz aumentada `𝕄` (no hay `Σ` no degenerada; ruta 6c abierta):

- **P3.2-A — sistema AUMENTADO completo** (`u,v,p,q₁,q₂,q₃` para a∈{u,v}): comprobar las RESTRICCIONES
  auxiliares, el residual algebraico de Schur y la equivalencia de la eliminación. Objeto: que las
  filas de restricción + la eliminación reproducen `K_P` numéricamente (incl. convergencia en N,L,BC).
  NO es una prueba de pseudo-hermiticidad.
- **P3.2-B — operador REDUCIDO `M_eff`** (con `η`): medir el DEFECTO INTERIOR de pseudo-herm
  `max‖M_eff†η−ηM_eff‖/‖M_eff‖` en región central. **Etiqueta obligatoria:** es un **diagnóstico del
  bulk**, NO una demostración de pseudo-herm para BC de scattering abiertas; debe acompañarse de la
  identidad de Green completa (`PADE_SPEC §9.3`) y de la comprobación del flujo. Convergencia en N,L,BC.

- **P3.4 (aparte): pseudo-unitaridad de S** `S†ηS=η` — es la comprobación de scattering abierto, no se
  mezcla con P3.2-B.

Precaución transversal: NO confundir pseudo-herm (infraestructura) con corrección física (Gate 2b/P3.3).

## Referencias
- `PADE_SPEC.md §9` — dominios, Green, BC, flujo (Paso 2 cerrado).
- `gates/GATES_REGISTRY.md` — P2.1–2.7 (PASS), P3.1 (PASS), P3.2–3.5 (pendientes).
- `NONLOCAL_SPEC.md §3–4` — operador físico y métrica η.
