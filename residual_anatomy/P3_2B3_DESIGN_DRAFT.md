# P3.2-B3 — DISEÑO (papel) η-pseudo-hermiticidad de M_eff en 𝒢=diag(W,−W) — NO CONGELAR, NO IMPLEMENTADO

**Estado:** DISEÑO en papel. NO hay implementación; NO toca v1.0 ni los gates aprobados
(`gate_P3_2B2_green.py`, `gate_P3_2B2c_continuum.py`). Depende de P3.1, P3.2-A, P3.2-B1 (APROBADOS) y
del cierre documental de B2 (opción B) y B2c (APROBADO). Redactado a petición de ChatGPT con sus 5
puntos. Marcas: **[PROBADO]** analítico/numérico aquí · **[HIPÓTESIS]** a verificar en el gate ·
**[ABIERTO]** sin respuesta.

Regla heredada (cautela viva #3): la métrica **NO se elige para forzar** `M†𝒢=𝒢M`; es la métrica
FÍSICA (norma BdG Q). Se mide el defecto TOTAL que salga, sin imponer de antemano la ley O(dx) de C.

---

## B3.0 — Espacio, producto interno y métrica  (punto 1 de ChatGPT)

Estado aumentado NO; aquí el operador FÍSICO BdG 2N (bloques `u=[0:N]`, `v=[N:2N]`, como `build_MBdG`).
Operador (de `P3_2_AUGMENTED_DRAFT.md §1`, con el cinético Padé `K_P` en vez del local):

    L_± = K_P ± C + (G + Qp) ,     M = [[ L_+ , G ],[ −G , −L_− ]]

- `K_P = ½D + c₁ DAD + c₂ D²A³D²`,  c₁=(λ−1)/4, c₂=−(λ−1)²/16, A=(2G+D)⁻¹, D=−∂².
- `C = −i(v∂_x + v'/2)` convectivo, **formalmente hermítico** (`C†=C`), símbolo +vk en L_+.
- `G+Qp` multiplicador real diagonal.

**Métrica física** (de `P3_2B_DESIGN_DRAFT.md §B0`, almacenamiento por bloques):

    𝒢 = diag(W, −W) = η ⊗ W ,   η=diag(1,−1),   W = diag(½,1,…,1,½)  (peso trapezoidal de ∫dx)

Producto interno (INDEFINIDO, espacio de Krein):  ⟨Φ,Ψ⟩_𝒢 = Φ†𝒢Ψ.  Norma BdG discreta:
`Q(Ψ) = Ψ†𝒢Ψ ≈ ∫(|u|²−|v|²)dx`. **[PROBADO]** `𝒢†=𝒢`, `𝒢` real, bloque-diagonal.

**Realización ensayada (declarada explícitamente, punto 3 de ChatGPT):**
- **Discretización de D/A/K_P:** matriz simétrica de **Dirichlet por nodos fantasma** (la misma de B2c /
  B1b-plana), NO la realización Robin de §9.2 (que rompe `K_Pᵀ=K_P`). Con ella `K_P†=K_P` euclídeo.
- **Perfil de G(x):** dos ensayos, reportados por separado — **(i) homogéneo G=1** (baseline algebraico,
  C con v constante → v'=0 trivial) y **(ii) inhomogéneo** (fondo waterfall de v1.0, v(x) y G(x)
  variables → C≠0 no trivial con v'≠0; es el caso con física de horizonte). El perfil se IMPONE (no se
  re-resuelve), como en NONLOCAL_SPEC.
- **PRE-REQUISITO (gate propio, antes de la descomposición de R):** verificar `‖K_Pᵀ − K_P‖/‖K_P‖ ≤ 1e−12`
  PARA ESA realización y esos N,L (reproduce B1b: 4e−13). Si falla → la descomposición de B3.1-2 no
  aplica (K_P no euclíd-s.a.) → PARAR y re-abrir B1 (criterio R2). **[HIPÓTESIS hasta ejecutarlo.]**

El objetivo de B3 es la **η-pseudo-hermiticidad** de M en esta métrica:

    M es 𝒢-pseudo-Hermítico  ⟺  M†𝒢 = 𝒢M  ⟺  R := M†𝒢 − 𝒢M = 0.

**Alcance (corregido tras ChatGPT, no sobre-afirmar):** en el continuo, la estructura FORMAL del operador
sugiere conservación de la norma de Krein Q **bajo las hipótesis adecuadas de dominio y condiciones de
contorno** (en particular para los términos NO locales de K_P). Este gate **NO demuestra autoadjunción ni
conservación exacta en el continuo**; B2 tampoco la demostró. B3 mide el defecto DISCRETO `R≠0` de una
realización concreta (abajo); la conservación de S se verá en P3.4.

---

## B3.1 — Defecto total de pseudo-hermiticidad, descompuesto por bloques  (punto 2)

**NO se presupone que R escale como O(dx).** Se calcula el bloque completo y se descompone.

**Derivación explícita (desde la definición, con † = transpuesta conjugada; C complejo, C†=C; G real
diagonal → G†=G):**

    M  = [[ L_+ ,  G ],[ −G , −L_− ]] ,   𝒢 = [[ W , 0 ],[ 0 , −W ]]
    M† = [[ L_+† , −G ],[ G , −L_−† ]]                         (porque (−G)†=−G, G†=G)
    M†𝒢 = [[ L_+† ,−G ],[ G ,−L_−† ]]·[[W,0],[0,−W]]
        = [[ L_+† W ,  (−G)(−W) ],[ G W , (−L_−†)(−W) ]]  = [[ L_+† W ,  G W ],[ G W ,  L_−† W ]]
    𝒢M  = [[W,0],[0,−W]]·[[ L_+ , G ],[ −G , −L_− ]]
        = [[ W L_+ , W G ],[ (−W)(−G) , (−W)(−L_−) ]]    = [[ W L_+ , W G ],[ W G ,  W L_− ]]

    ⇒  R = M†𝒢 − 𝒢M = [[ L_+†W − W L_+ ,   GW − WG ],
                        [ GW − WG       ,   L_−†W − W L_− ]]

**El signo del bloque inferior es `+(L_−†W − W L_−)`**, NO `−`: las dos negaciones (`−L_−` en M y el
`−W` del bloque inferior de 𝒢) se **cancelan**. **[PROBADO]** numéricamente (N=6, C hermítico complejo,
K_P euclíd-sim real, G,G+Qp diagonales): los 4 bloques de `R` coinciden con esta fórmula a ~máquina, y
el bloque (2,2) reproduce `L_−†W − W L_−` con signo +.

**Bloques fuera de diagonal:**  `GW − WG = 0` **EXACTO** (G y W diagonales → conmutan; verificado 0.0).
**[PROBADO]**

**Bloques diagonales**, con `L_± = K_P ± C + (G+Qp)` y linealidad:

    L_±†W − W L_± = (K_P†W − W K_P)  ±  (C†W − W C)  +  ((G+Qp)W − W(G+Qp))

1. **(G+Qp):**  `(G+Qp)W − W(G+Qp) = 0` **EXACTO** (multiplicador diagonal real). **[PROBADO]**
2. **K_P:**  `K_P†W − W K_P`. Con la **discretización PLANA simétrica** de B1b, `K_P†=K_P` euclídeo
   (`K_Pᵀ=K_P`, medido 4e−13 en B1b). Entonces `K_P†W − W K_P = [K_P, W]`, con entradas
   **`[K_P,W]_{ij} = (K_P)_{ij}(W_j − W_i)`**. Como `W_j−W_i = 0` salvo cuando EXACTAMENTE uno de i,j es
   un extremo (W=½) y el otro interior (W=1), hay que **distinguir dos afirmaciones** (precisión pedida
   por ChatGPT):
   - **(2a) [PROBADO]** El **bloque INTERIOR×INTERIOR** (ambos índices interiores) se anula
     **EXACTAMENTE** (`W_i=W_j=1` ⇒ diferencia 0). Verificado max=0.0 (N=6).
   - **(2b) [PROBADO]** Las entradas que **conectan un nodo de borde con el interior NO son nulas** en
     general: `[K_P,W]_{0,j} = (K_P)_{0,j}(1−½) = ½(K_P)_{0,j} ≠ 0` (verificado: 0.2766, 0.1285, …).
     Decaen como `(K_P)_{ij} ~ e^{−κ|i−j|}` (colas de A³) al alejarse del extremo, pero NO son
     idénticamente cero.
   ⇒ el defecto de K_P está **localizado en las dos "cruces" de borde** (fila+columna de cada extremo),
   con amplitud decreciente hacia el interior; NO es estrictamente "2 nodos diagonales".
   **CLAVE (corrección ChatGPT):** cada entrada de la cruz tiene SIEMPRE **al menos un índice de borde**
   (fila 0 o fila N−1), por lejos que llegue su cola. Por tanto el **bloque interior×interior es cero
   EXACTO con independencia de la distancia a los extremos**; NO hace falta ningún margen `m·dx≳1/κ`
   para excluir la cola: basta definir `R_bulk` por índices **estrictamente interiores** (B3.2) y las
   cruces caen enteras en `R_borde`. **[HIPÓTESIS]** a re-verificar (pre-requisito B3.0): que `K_P†=K_P`
   euclídeo se mantenga al ensamblar M (si no → re-abrir B1, criterio R2).
3. **C:**  `C†W − W C`. Es el **único término que puede dar defecto de BULK**. En v1.0 la discretización
   de C es hermítica **solo a O(dx)** (B1b: `C` O(dx) euclídeo). ⇒ se **PREDICE** un **piso O(dx)** en
   el bulk, procedente de C. **NO se impone**: se mide el R total y se ajusta su orden. **[ABIERTO]** el
   valor real del piso y su constante.

**Consecuencia de diseño:** el defecto de bulk de M está dominado por C. K_P y (G+Qp) no contribuyen al
bulk (exacto / ~máquina); el borde es un término separado que se reporta aparte, no se mezcla con el
bulk. Esto respeta la cautela #1 (C herm. solo a O(dx) = piso de B3) y la #3 (no forzar la métrica).

**Verificación estructural [PROBADO]** (bloques, N=6, K_P euclíd-sim REAL, **C hermítico COMPLEJO**
(C†=C), G,G+Qp diagonales reales): los 4 bloques de R coinciden con la fórmula a ~máquina; `GW−WG`=0
exacto; `(G+Qp)` parte=0 exacto; parte K_P = [K_P,W] con interior×interior=0 exacto y cruces de borde no
nulas; parte C = O(1) estructural (el valor O(dx) es propiedad de la discretización REAL de C de v1.0,
a medir, no de esta prueba con C aleatorio).

---

## B3.2 — Normas y normalizaciones  (punto 3)

Se reportan, SIN mezclar regiones:

- **Defecto absoluto por bloque:** `‖R_bloque‖` en (a) norma espectral ‖·‖₂ y (b) máximo entrada a
  entrada ‖·‖_max. Las dos, por la cautela de B1 (la simetría depende de la métrica; no colapsar a un
  único número).
- **Defecto relativo:** `‖R‖ / ‖𝒢M‖` (misma norma), para hacerlo adimensional y comparable entre N.
- **Separación BULK vs BORDE (definición EXACTA en la matriz 2N×2N, corregida por ChatGPT):**
  - **Índices de borde** (4, por los dos bloques): `B = {0, N−1, N, 2N−1}` (donde W=½ en cada bloque).
  - **Índices interiores:** `I = {1,…,N−2} ∪ {N+1,…,2N−2}`.  `R_bulk := R[I, I]`; `R_borde` = el resto
    (toda fila o columna con índice en B) — incluye las cruces de `[K_P,W]` enteras, sin margen extra.
  - **Conservar la separación u/v al tomar normas:** reportar `R_bulk` por sub-bloques `uu` (=`[1:N−1,
    1:N−1]`) y `vv` (=`[N+1:2N−1, N+1:2N−1]`), NO solo la norma global de `R_bulk` (una norma única puede
    ocultar qué bloque diagonal aporta el defecto).
  - Por B3.1-2 la contribución de K_P a `R_bulk` es **cero exacto** (algebraico, condicionado al
    pre-requisito `K_P†=K_P`); **queda `R_bulk = R_bulk[C]` SÓLO si se cumplen las hipótesis** (K_P
    euclíd-sim, W=I interior, G+Qp real diagonal) — es consecuencia ALGEBRAICA, NO resultado
    experimental: se CONFIRMA ejecutando la descomposición por términos, no se presupone.
  Se reportan por separado. **El criterio de aceptación mira `R_bulk`**; `R_borde` se documenta pero NO
  entra en la ley de convergencia del bulk.
  - Margen extra `m` (OPCIONAL, solo DIAGNÓSTICO de sensibilidad espacial): puede estudiarse ampliar el
    margen más allá de los nodos de borde para ver si el piso cambia, pero **NO es necesario para excluir
    `[K_P,W]`** (ya excluido por la definición estrictamente interior). Un cambio al ampliar `m` NO se
    atribuye automáticamente a colas de K_P (no las hay en el interior); se interpreta como estructura
    espacial de C cerca del borde.
- **Por términos:** además del R total, se reporta la descomposición medida `R = R[K_P] ± R[C] +
  R[G+Qp]` para confirmar empíricamente que el bulk = R[C] (y que R[K_P], R[G+Qp] son ~máquina en bulk).

Normalización de campos: NINGUNA sobre M (es el operador, no un modo). El defecto es una propiedad de
la matriz M y la métrica 𝒢, independiente de estados. (La normalización de modos a flujo unidad —el bug
de la ocupación en op='lambda'— es de P3.3/P3.4, NO de B3.)

---

## B3.3 — Ensayos de convergencia  (punto 4)

Barridos (todos con la realización Dirichlet-fantasma de B3.0, y los dos perfiles de G declarados allí:
(i) homogéneo G=1 y (ii) inhomogéneo waterfall; previo a cada uno, el pre-requisito `K_Pᵀ=K_P` de B3.0):

- **En N (dx→0), L fijo:** medir `‖R_bulk‖(dx)` y ajustar pendiente log-log sobre ≥4 puntos (como B2c,
  sin descartar puntos). Reportar `‖R_bulk‖` de K_P y (G+Qp) por separado (deben quedar a nivel de
  máquina, sin pendiente). **Papel de los dos perfiles (cautela ChatGPT):**
  - **(i) G=1 homogéneo con v CONSTANTE = CONTROL NEGATIVO.** Aquí v'=0 y C puede ser hermítico también
    en la discretización ⇒ `R_bulk` podría ser **~0 (sin pendiente)**. **NO se exige pendiente 1 en este
    ensayo**; su papel es confirmar que, sin fuente de defecto, el bulk es ~máquina (si NO lo es, hay un
    bug en el ensamblado, no física).
  - **(ii) waterfall inhomogéneo = ENSAYO PRINCIPAL.** Aquí v(x) varía (v'≠0) y la discretización de C
    introduce el defecto que se quiere medir ⇒ es DONDE se evalúa la predicción **pendiente ≈ 1 (O(dx))**.
- **En L (dominio), N/L fijo (dx fijo):** comprobar que R_bulk NO depende de L (el defecto es local),
  y que R_borde queda confinado en los extremos.
- **Condicionamiento:** registrar `cond(2G+D)` y el N a partir del cual la A³ DENSA se corrompe (como en
  B2.3, N≲3200); por encima, usar realización en banda. No confundir piso físico O(dx) con ruido de
  condicionamiento.

---

## B3.4 — Criterios de aceptación / rechazo (ANTES de ejecutar)  (punto 5)

Fijados a priori para no ajustar el criterio al resultado:

**ACEPTAR B3** si TODO:
- (A1) `R` fuera de diagonal = 0 a ~máquina (≤1e−12 rel). [predicho exacto]
- (A2) `R_bulk[K_P]` y `R_bulk[G+Qp]` ≤ 1e−10 rel, sin pendiente en dx (planos). [predicho ~máquina]
- (A3) **[REFORMULADO en DOS criterios tras ejecutar, ChatGPT 2026-11-15 — ver nota de resultados]:**
  - **(A3-matriz)** defecto de `R_bulk` en norma **max/espectral**: se REGISTRA, NO se exige que desaparezca.
    Resultado: **O(1) constante** (artefacto de rejilla: `C†−C ~ v'` solo actúa en modos k~1/dx).
  - **(A3-suave)** defecto proyectado al **subespacio SUAVE ortonormalizado (QR)** — paquetes de espectro
    concentrado a baja k, localizados en el horizonte: se exige **O(dx²)** (NO O(dx): C es hermítico en el
    continuo y usa FD centrada O(dx²)). En el **control (i) G=1/v cte**, `~máquina` (sin pendiente). El
    `R_bulk ≈ R_bulk[C]` se confirma por la descomposición por términos (A2), no se presupone.
- (A4) [reformulado, ChatGPT — para que NO sea cierto por definición] DOS condiciones:
  (A4a) `R_bulk` **converge y es estable al aumentar L** con dx y perfil físico fijos (el defecto es
  local, no crece con el dominio);  (A4b) la contribución de `C` en las **primeras filas interiores** se
  caracteriza por separado y **no domina artificialmente** el `R_bulk` (que el piso no venga solo de 2–3
  nodos junto al borde). NO se exige que `R_borde` sea pequeño (las cruces de K_P pueden tener amplitud
  finita).

**RECHAZAR / PARAR** si:
- (R1) `R_bulk` NO converge, o pendiente < 0.5 (defecto O(1) en bulk) → la realización de K_P/C en el
  operador completo introduce un error no controlado; revisar discretización antes de seguir.
- (R2) `R_bulk[K_P]` o `R_bulk[G+Qp]` NO son ~máquina → la 𝒢-simetría supuesta (B1b) no se transfiere
  al operador completo; re-abrir B1.
- (R3) el defecto depende de L en bulk → no es local; la interpretación como "piso de borde + O(dx)
  interior" es falsa.
- (R5) al redefinir arbitrariamente la frontera interior/borde el piso de `R_bulk` cambia de forma
  incoherente con una estructura O(dx) de C cerca del borde (p.ej. salta de orden) → revisar el stencil
  de C en los primeros nodos interiores antes de concluir. (NO se atribuye a colas de K_P: por B3.1-2 no
  hay K_P en el interior estricto.)

**Diagnóstico (NO criterio de bulk), corrección de ChatGPT sobre el antiguo R4:** la comparación de `𝒢`
físico (W) con el euclídeo `diag(I,−I)` **NO distingue nada en el bulk**: como `W=I` en el interior,
`R_bulk^(W) = R_bulk^(I)` **idénticamente** (no es predicción asintótica; es álgebra). ⇒ el bulk es
**automáticamente independiente de la métrica**, lo que de hecho **refuerza** que el piso O(dx) medido no
es artefacto de la elección de métrica (cautela #3 satisfecha de oficio en el interior). La comparación
W vs I SÍ puede diferir en el **defecto total y de borde** (allí W≠I): se reporta como diagnóstico de
borde, no como prueba sobre el bulk.

**Qué NO prueba B3** (alcance, condición opción B): no demuestra autoadjunción (dominio maximal); no se
extrapola a otras BC, a otros perfiles fuera de los ensayados, ni se apoya en ninguna identidad SBP (que
no existe). B3 es un gate numérico de pseudo-hermiticidad del operador completo, con sus propios
controles, como exigió ChatGPT.

---

## Orden posterior (sin cambios)
B3 (este) → P3.4 (conservación de S, `S†ηS=η`, 3 columnas) → Gate 2b/P3.3 (regresión λ=1 vs benchmark
congelado) → medir n_L(λ), G fijo, λ∈[0.5,2]: **¿persiste el colapso?** (pregunta física, ABIERTA).

**Pendiente antes de IMPLEMENTAR B3:** visto bueno de ChatGPT a este diseño (en especial B3.1 — la
descomposición por bloques y el papel de C— y los criterios B3.4). No arrancar código hasta entonces.

---

## ✅ RESULTADO — B3 CERRADO/APROBADO (ChatGPT, 2026-11-15; alcance limitado)
Implementado en `gates/gate_P3_2B3_pseudohermiticity.py`. **Formulación admisible para el registro y la
documentación (redacción ChatGPT):**

> Para la realización discreta Dirichlet-fantasma estudiada, el defecto de pseudo-hermiticidad proyectado
> sobre el subespacio de paquetes suaves de baja frecuencia converge como **O(dx²)**. La descomposición
> por términos atribuye ese comportamiento al operador convectivo **C**; la contribución proyectada de
> **K_P** permanece a nivel de precisión numérica. El defecto matricial completo **NO desaparece** y
> persiste en modos de rejilla. El resultado es un **diagnóstico de pseudo-hermiticidad PROYECTADA, NO una
> demostración de pseudo-hermiticidad exacta del operador completo ni de autoadjunción en el dominio
> continuo.**

Números: A1 fuera-diag=0 exacto ✓; A2(G+Qp)=0 exacto ✓; A2(K_P) limitado por condicionamiento denso en
norma cruda pero proj_KP~máquina (1e−12); A3-matriz=O(1) (caracterizado); A3-suave (QR)=O(dx²) pend 2.00,
proj_C pend 2.00 y proj_TOTAL≡proj_C (recon_err=0); A4a estable en L. Tareas independientes pendientes:
(1) **P3.4** conservación de S; (2) realización en **BANDA** (dispersa) antes de escalar a N grande, con
sus propios controles de equivalencia/estabilidad/pseudo-herm.
