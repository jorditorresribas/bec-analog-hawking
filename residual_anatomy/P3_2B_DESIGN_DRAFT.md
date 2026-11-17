# P3.2-B — DISEÑO (paper) η-pseudo-herm de M_eff (v2, CORREGIDO) — NO CONGELAR

**Estado:** DISEÑO en papel, NO implementación. v2 tras la revisión de ChatGPT: corrige el orden del
producto tensorial y REFUTA la premisa de que `W` hacía autoadjunto `K_P`. Todo marcado [N] está
verificado numéricamente (2026-11-15). NO toca v1.0. Depende de P3.1, P3.2-A (APROBADOS).

## B0. Orden de los grados de libertad y métrica  [CORREGIDO]
Almacenamiento REAL (verificado en `build_MBdG`): **por BLOQUES**, `u=[0:N]`, `v=[N:2N]`. Por tanto la
métrica física NO es `W⊗η` (orden intercalado) sino:

    𝒢 = diag(W, −W) = η ⊗ W ,   W = diag(½,1,…,1,½) (peso trapezoidal de ∫dx)   [CORREGIDO]

`Q = ∫(|u|²−|v|²)dx ≈ Ψ†𝒢Ψ`. B0 verifica el orden y construye `𝒢` acorde.

## B1. Simetría ponderada de los bloques — el hallazgo que CORRIGE el diseño  [VERIFICADO]
**Mi premisa anterior era falsa.** Que `W` simetrice `A` NO implica que simetrice `K_P`: los factores
externos `D`, `D²` tienen su propia simetría. Medido (defecto relativo, N=120):

    operador     euclídeo(rel)   W(rel)
    D (plano)      0.00           0.25      ← D es EUCLÍDEO-s.a., NO W-s.a.
    A (Robin)      0.38           3e-16     ← A-Robin es W-s.a., NO euclídeo
    K_P (Robin)    0.51           0.22      ← ¡NI euclídeo NI W!  (mezcla D eucl + A-W)
    K_P (Dirichlet-identidad) 0.29 0.21     ← tampoco

**Conclusión:** con la discretización MIXTA (D plano + A con BC Robin/Dirichlet-identidad), `K_P` no es
autoadjunto en ninguna métrica simple. La `W`-simetría de `A` no se transfiere.

## B1b. RESOLUCIÓN (verificada): discretización PLANA simétrica + L grande  [VERIFICADO]
Si `A = (2G + D_plano)⁻¹` con `D_plano` SIN filas de BC (**truncamiento finito de la discretización de
recta** — NO demostrado equivalente a Neumann; una matriz plana sin modificar bordes puede ser simétrica
sin representar exactamente Neumann), entonces `2G+D_plano` es euclídeo-simétrico → `A` euclídeo-s.a. →
`K_P` **euclídeo-autoadjunto**:

    K_P eucl-sim(rel) = 4.4e-13 (N=120,L=12) ,  3.9e-13 (N=400,L=40)   [N, reproduce P2.3]

(Estos valores son el RESULTADO de esa prueba concreta con sus N,L,dx,κ; NO una cota general — el gate
los registra con N,L,dx,κ y región.)

La Robin/Dirichlet-identidad de P3.1 (útiles para reducir el ERROR de borde de A³) **rompen** la
autoadjunción de `K_P`. **Para la medida de pseudo-herm (P3.2-B) se usa la discretización PLANA
simétrica a L grande**, donde `K_P` es euclídeo-s.a. sobre TODA la matriz (sin necesidad de restringir
al interior, evitando el problema de §B3b). En el bulk `W≈I`, así que euclídeo ≈ métrica `𝒢`.

**CONTROL DE CONSISTENCIA (tensión P3.1↔P3.2-B) — RESUELTO [N]:** autoadjunta ≠ correcta, así que se
re-pasó el test Green de P3.1 a la PLANA: `A³` plana vs `e^{−κ|x−y|}/2κ` reproduce el resolvente de recta
completa en el interior a **6.7e−7** (= piso FD O(dx²), igual que P3.1) para L=10,20,40. ⟹ la plana es
**autoadjunta Y correcta** (reproduce A³), no solo simétrica. La tensión queda cerrada.

## B2. Identidad de Green DISCRETA (control independiente)
Verificar `⟨w, K_P u⟩ − ⟨K_P w, u⟩` (producto de la MISMA discretización plana) contra los términos de
borde de `PADE_SPEC §9.3` (`½𝒲+c₁[…]+c₂[𝒲₂…]`). Con la plana simétrica a L grande y campos localizados,
el defecto → 0 (términos de borde `O(L²e^{−κL})`). Para ondas planas, el borde = flujo `ε_P'(k)|U|²`
(§9.7). NO basta comparar expresiones continuas: misma malla y mismo producto escalar.
**DOS comprobaciones DIFERENCIADAS (ChatGPT), no la misma:** (1) **simetría matricial** `K_P†=K_P` (si la
plana es simétrica, esto da ~máquina — verifica la simetría discreta, NADA más); (2) **correspondencia
con Green**: el término de borde discreto CONVERGE al término continuo previsto de §9.3 para los campos de
prueba (esto sí prueba que la discretización reproduce cada término de la identidad). No confundirlas.

## B3. Pseudo-hermiticidad de M_eff en la métrica física
`M_eff = [[L_+, G],[−G, −L_−]]`, `L_± = K_P ± C + (G+Qp)`, `𝒢 = diag(W,−W)`. Condición:

    M_eff† 𝒢 = 𝒢 M_eff   ⟺   L_+† W = W L_+ ,  L_−† W = W L_− ,  (G diagonal real: bloques compatibles)

Con la discretización plana (W≈I bulk), equivale a `L_±` euclídeo-s.a. en el bulk. **Cautela (ChatGPT):
la autoadjunción euclídea de `K_P` NO implica la ponderada de `L_±`** — con `W` trapezoidal la diferencia
en los extremos persiste; el gate la MIDE, no la descarta por definición. Factores a verificar por
SEPARADO (Control 0), NO asumir: `K_P` (B1b: OK plana euclídeo); `G+Qp` (diagonal real, trivial);
`C=−i(v∂+v'/2)` — **hallazgo B1 [N]:** la forma de `build_MBdG` v1.0 (`diag(v)D1+½diag(v')`) es hermítica
solo a **O(dx)** en norma relativa de operador (NO de borde, NO O(dx²); `max|C−C†|~O(v')` abs pero
`max|C|~O(v/dx)`). ⟹ el defecto de pseudo-herm de `M_eff` tendrá un **piso O(dx) dominado por `C`**, que
CONVERGE linealmente. La forma exactamente antisimétrica `−i/2(vD1+D1v)` da `C` hermítica a máquina [N]
pero difiere de v1.0 a O(dx) (dentro del error de discretización). B3 usará la forma de v1.0 por
consistencia con el benchmark y reportará la convergencia O(dx); opcionalmente contrastará con la exacta.

### B3b. Defecto de bulk vs borde — con la cautela de ChatGPT
Calcular el defecto matricial COMPLETO `R = M_eff†𝒢 − 𝒢M_eff` y DESPUÉS separar bloques
interior-interior, interior-borde, borde-interior (NO construir el operador restringido al interior y
adjuntarlo: `interior(prod) ≠ prod(interior)` para operadores no locales — verificado que con Robin el
"interior" de K_P da defecto espurio por acoplo con el borde). Con la plana simétrica el defecto es
~máquina en TODA la matriz, así que la separación es limpia; la Robin se deja solo como diagnóstico.

## Controles (secuencia revisada de ChatGPT)
- **B0:** orden de DOF + construcción de `𝒢=diag(W,−W)`. [hecho en diseño; a codificar]
- **B1 (Control 0, OBLIGATORIO):** simetría euclídea y `W` de `D`, `A`, `K_P`, `C` POR SEPARADO, para
  plana / Robin / Dirichlet-interior. Localiza defectos. [VERIFICADO el núcleo; a consolidar en gate]
- **B2:** identidad de Green discreta (misma malla y producto). [a codificar]
- **B3:** pseudo-herm de `M_eff` en `𝒢`, defecto de bulk → 0, convergencia N,L. [a codificar]
- **P3.4 (aparte):** `S†ηS=η`.

## PROBADO vs PENDIENTE (de este diseño v2)
| Punto | Estado |
|---|---|
| Orden por bloques → `𝒢=diag(W,−W)` | **VERIFICADO** (corrige `W⊗η`) |
| `W` simetriza A-Robin pero NO K_P | **VERIFICADO** (refuta premisa v1) |
| Plana simétrica + L grande → K_P euclídeo-s.a. | **VERIFICADO** (4e-13, reproduce P2.3) |
| `interior(prod) ≠ prod(interior)` (Robin) | **VERIFICADO** (defecto espurio) |
| Pseudo-herm de M_eff completo (con C, borde) | **HIPÓTESIS** — gate B3, pendiente |
| Identidad de Green discreta | **HIPÓTESIS** — gate B2, pendiente |

## Lo que este diseño NO reclama
- NO autoadjunción con maximalidad de dominio (abierto). NO equivalencia de dominios inhomogéneos.
- NO scattering (P3.3 + P3.4). La métrica/discretización se ELIGE por consistencia (plana simétrica),
  verificando que da K_P autoadjunto — NO para forzar el resultado.

## Referencias
- `PADE_SPEC §9.3, §9.7`; `P3_2_AUGMENTED_DRAFT §6b`; `gates/GATES_REGISTRY.md`.
