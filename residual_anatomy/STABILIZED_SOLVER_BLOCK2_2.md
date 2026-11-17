# B2.2 — Diseño del estudio de CONDICIONAMIENTO de 𝓜_aug (CERRADO EN PAPEL)

**Fase:** solver de scattering estabilizado (`STABILIZED_SOLVER_RESERVATION.md`), bloque 2.
**Fecha:** 2026-11-15. **Estado: DISEÑO CERRADO EN PAPEL** (revisión cruzada ChatGPT).
**NO autoriza implementación.** Tolerancias numéricas: **pendientes justificadas** (se fijan antes de
escribir el gate). No toca v1.0, ni los gates/docs del operador cerrados, ni `GATES_REGISTRY.md`,
ni `STABILIZED_SOLVER_BLOCK1.md` (congelado), ni el gate B2.1 (cerrado).

> **Cuestión central que B2.2 debe resolver:** demostrar que el solver calcula **amplitudes de
> scattering estables y reproducibles al refinar la malla**, *incluso si* el sistema algebraico
> completo está mal condicionado. No es "conseguir un número de condición bonito".

---

## 0. Objetivo y separación de ejes N / X0

Dos fuentes de condicionamiento, **ejes distintos**:
- **Eje N** (dx→0, **X0 fijo**): discretización. **Es lo que B2.2 estudia.**
- **Eje X0** (tamaño de dominio + cola evanescente `~e^{2κX0}`): **gate de ventana de X0, separado.**

**Matiz:** que `e^{2κX0}` sea constante al variar N **no** lo hace inocuo — multiplica el
condicionamiento base por una constante (posiblemente grande) y encarece el solve. Por eso B2.2
**fija y documenta X0**; solo afirma que su *variación* es otro eje, no que su influencia desaparezca.

## 1. Escalas de bloques de 𝓜_aug (ESTRUCTURAL)

| Bloque | Entradas | Escala |
|---|---|---|
| Interior (caja) `(I∓(dx/2)A)` | identidad O(1) ± O(dx·‖A‖) | fila O(1) |
| Matching `Y_b−Ξ_b a_b` | identidad O(1) y columnas de Ξ | Ξ mezcla O(1) (u,v) con **O(k)** (u′=iku, v′=ikv) |
| Restricción | un único 1.0 | fila O(1); `b` solo puebla estas filas |

`b` solo en filas de restricción ⇒ la sensibilidad relevante vive en **columnas concretas de 𝓜⁻¹**
(motiva §4). Candidato a artefacto de escala: derivadas O(k) vs amplitudes O(1).

## 2. Equilibrado: crudo → Ruiz → físico (los tres, en ese orden)

1. **Cruda:** referencia.
2. **Ruiz:** equilibrado iterativo de dos lados; diagnóstico genérico del desequilibrio algebraico.
3. **Físico:** no-dimensionalización coherente aplicada a **incógnitas Y ECUACIONES** (no solo a las
   columnas de derivada), vía
   $$\widetilde Y=(u,\ \ell\,u',\ v,\ \ell\,v')^{\mathsf T}$$

**Escala física fija (DECISIÓN):**
$$\boxed{\ \ell=\xi_u=\frac{1}{\sqrt{2\,g\,n_u}}=\frac{1}{\sqrt2\,c_u}\ \stackrel{g=n_u=1}{=}\ \frac{1}{\sqrt2}\approx0.7071\ }$$
Longitud de healing **asintótica upstream**, de los parámetros nominales del fondo. **Independiente de
N, dx y X0.** Se decide ANTES del barrido. (`1/κ` NO se usa como escala primaria: depende del modo y de
ω; sirve para la cola evanescente / eje X0.)

Aplicación coherente: transformar también las **filas de matching** y mantener bien definida la relación
entre coeficientes modales y campos — no basta reescalar columnas de derivada.

**Matiz:** `cond(𝓜)/cond(D_r𝓜D_c)` es **descriptivo, no invariante** del artefacto — el equilibrado
puede cambiar cond sin mejorar la precisión de los observables. **Juez real = §4 y §6**, no §2.

## 3. Condicionamiento intrínseco vs pérdida de rango (sobre matrices EQUILIBRADAS)

- Ley de crecimiento en N: ajustar `cond(equilibrado)≈C·N^p`, reportar `p`; ¿**satura** hacia constante
  fijada por X0 (benigno) o crece sostenido (patología)?
- **Rango:** `σ_min` frente a `n_dof·ε_mach·σ_max` es **referencia de rango numérico, NO frontera
  universal** (depende de algoritmo, precisión y escalado). Por tanto: calcular `σ_min` **sobre las
  matrices equilibradas** y **documentar la convención** del umbral. Discriminador: `σ_min` muchos
  órdenes sobre el piso = condicionamiento (fenómeno 2); hacia el piso = pérdida de rango (fenómeno 3).
- **Alcance:** el test `ker=∅` de B2.1 es un **test de rango a precisión finita**, no prueba de unicidad
  del problema **continuo**.

## 4. Sensibilidad de las amplitudes vía adjunto (2-norma; b y M separadas)

Normas: **espectral/2-norma** en todo (`‖𝓜‖₂=σ_max`). Sea `𝓜x=b`, salida escalar `s=c†x`.
Adjunto `𝓜†y=c` ⇒ `s=c†𝓜⁻¹b=y†b`.

**Perturbación de b (aditiva):**
$$|\delta s|\le\|y\|_2\|\delta b\|_2\ (\text{abs});\qquad
\boxed{\ \kappa_{\mathrm{out},b}=\frac{\|y\|_2\,\|b\|_2}{|s|}\ }\ (s\neq0,\ \text{rel})$$

**Perturbación de la matriz 𝓜 (1er orden):** de `δx≈−𝓜⁻¹δ𝓜\,x`, `δs=c†δx`:
$$\delta s\simeq-\,y^\dagger\,\delta\mathcal M\,x\ \Rightarrow\ |\delta s|\le\|y\|_2\|\delta\mathcal M\|_2\|x\|_2;\qquad
\boxed{\ \kappa_{\mathrm{out},M}=\frac{\|y\|_2\,\|\mathcal M\|_2\,\|x\|_2}{|s|}\ }\ (s\neq0)$$
(se omiten términos de orden superior en δ𝓜.)

**Backward error (coherente):**
$$\boxed{\ \eta_{\mathrm{back}}=\frac{\|\mathcal Mx-b\|_2}{\|\mathcal M\|_2\|x\|_2+\|b\|_2}\ }$$
Pequeño ⇒ backward-stable, **NO** garantiza error forward pequeño.

**Advertencias:** (i) `s→0` infla la condición relativa aunque el error absoluto sea pequeño → reportar
**absoluta Y relativa**; (ii) **mantener separadas** δb y δ𝓜. `κ_out` suele ser `≪ cond(𝓜)`.

**Multiplicidad (DECISIÓN):** hay **3 salidas × 3 problemas (RHS)**. Se calcula la sensibilidad para
**cada amplitud de salida y cada RHS**, no para una representativa. **Cada salida tiene su propio
adjunto `y`.**

**Selector bajo equilibrado (precisión):** si `𝓜̃=D_r𝓜D_c` con `x=D_c z`, entonces `s=c†x=(D_c†c)†z`
⇒ en variables nuevas el selector es `D_c†c` y el adjunto `𝓜̃†w=D_c†c`. Se mide así la salida
**física** correcta, no una transformada.

## 5. Backward / forward error y refinamiento

- `η_back` (§4) como estabilidad backward.
- **Refinamiento iterativo** con la **misma** matriz discretizada; comprobar si reduce el error
  estimado. **No** confundir con convergencia al continuo (esa = comparar distintos N, §6).

## 6. Consistencia de S: crudo/equilibrado y entre mallas (protocolo)

**A. Crudo vs equilibrado** (mismo N, X0, ω):
- Resolver con las 3 representaciones (cruda, Ruiz, física).
- **Recuperar siempre las amplitudes en la convención física original** (deshacer transformaciones).
- Comparar S tras deshacer: **diferencia absoluta máxima por elemento** + diferencia relativa Frobenius.
- Elementos de S casi nulos → informar **diferencia absoluta**, no solo relativa.

**B. Convergencia en N** (X0, perfil, ω, normalizaciones fijos; mismas ω en todos los refinamientos):
- Comparar `S_N(ω)` entre mallas consecutivas y respecto de la más fina.
- Absolutas y relativas, con el mismo trato de elementos próximos a cero.
- Métrica relativa matricial:
  $$\boxed{\ \Delta_S(N_a,N_b)=\frac{\|S_{N_a}-S_{N_b}\|_F}{\max(1,\ \|S_{N_b}\|_F)}\ }$$
  (evita que un S pequeño produzca comparación relativa artificialmente enorme; **no** sustituye el
  examen elemento a elemento.)

**Pseudo-unitaridad de S = comprobación DISTINTA** (fase posterior): que dos soluciones coincidan no
demuestra que ambas sean físicamente correctas.

## 7. Criterios de aceptación a priori (sin fijar tolerancias)

- **(C1)** El equilibrado elimina el desequilibrio algebraico (`drop` caracterizado) — descriptivo,
  confirmado por §4/§6.
- **(C2)** `cond(equilibrado)` satura o crece suave/polinómicamente en N a X0 fijo.
- **(C3)** `σ_min(equilibrado)` órdenes sobre el piso (convención documentada) → sin pérdida de rango.
- **(C4)** `κ_{out,b}`, `κ_{out,M}`, `η_back` pequeños y **estables en N**; S coincide crudo-vs-
  equilibrado y converge en N (`Δ_S`) → **amplitudes estables aunque cond(𝓜) grande**. ← CENTRAL.
- **(C5)** Tolerancias numéricas: se fijan **antes de implementar**, justificadas con la precisión
  alcanzable (como cada gate previo).
- **Ejes:** B2.2 fija y documenta X0, barre N; el eje X0 (`e^{2κX0}`) es del gate de ventana de X0.

---

## Decisiones congeladas

| Elemento | Decisión |
|---|---|
| Escala física ℓ | `ξ_u = 1/√(2 g n_u) = 1/√2 ≈ 0.7071` (nominal upstream, indep. de N) |
| Variable adimensional | `Ỹ=(u, ℓu′, v, ℓv′)`, aplicada a incógnitas **y** ecuaciones |
| Normas | 2-norma (espectral) en todas las métricas |
| Adjunto | **uno por amplitud de salida** (3 salidas × 3 RHS) |
| Selector bajo equilibrado | `D_c†c` (medir la salida física correcta) |
| Comparación crudo/equilibrado | mismo problema, amplitudes recuperadas en convención física |
| Convergencia de S | barrido en N, X0 fijo; métrica `Δ_S` + examen elemento a elemento |
| Elementos de S ≈ 0 | comparación absoluta |
| Tolerancias | pendientes justificadas, antes de implementar |

## 8. TOLERANCIAS (CONGELADAS: forma/orden; constantes PROVISIONALES)

**Marco:** `n_dof·ε_mach` es **escala de referencia teórica, NO cota universal** (`γ` LU depende de
factores/escalado/aritmética; las anclas de B2.1 son un caso concreto). Los umbrales algebraicos se
**calibran empíricamente** con el backward error **observado** en los solves de control.

### Nivel 1 — ALGEBRAICAS (aceptación, calibración empírica)
En **variables transformadas** `Ỹ=(u,ℓu′,v,ℓv′)` para no mezclar escalas campo/derivada.

| Métrica | Qué mide | Umbral | Tipo·alcance |
|---|---|---|---|
| `η_back=‖𝓜x−b‖₂/(‖𝓜‖₂‖x‖₂+‖b‖₂)` | backward stability (NO precisión forward) | **calibrado** del observado; `n_dof·ε`=referencia | aceptación, (ω,RHS)+máx |
| `r_match,b=‖Ỹ_b−Ξ̃_b a_b‖₂` **(+ versión relativa, denominador documentado)** | matching en vars transformadas | calibrado del observado | aceptación, por borde |
| `\|a_in−s\|` | prescripción entrantes (≠matching) | calibrado (ref `n_dof·ε·\|s\|`) | aceptación, por in |
| `r_grow=‖Ξ̃_grow a_grow‖₂/max(‖Ỹ_b‖₂, escala_ref)` | **contaminación del CAMPO por el modo excluido** (no el coef., que depende de la convención modal) | calibrado (≈cero-máquina) | aceptación, por grow |

### Nivel 2 — ESTABILIDAD (diagnóstico + alarma)
Sensibilidades: conservar **absolutas y relativas**; **sin umbral de aceptación sobre `κ_out` aislado**.
$$|\delta s|_b\le\|y\|_2\|\delta b\|_2,\quad |\delta s|_M\lesssim\|y\|_2\|\delta\mathcal M\|_2\|x\|_2\ (\text{1er orden}),\quad
\kappa_{\mathrm{out},b}=\tfrac{\|y\|_2\|b\|_2}{|s|},\ \kappa_{\mathrm{out},M}=\tfrac{\|y\|_2\|\mathcal M\|_2\|x\|_2}{|s|}$$

| Métrica | Detecta / no concluye | Umbral | Tipo |
|---|---|---|---|
| `\|δs\|` abs + `κ_out` rel (adjunto por salida) | amplificación por salida; `s→0` infla la rel → informar **absoluta** | informar abs+rel + **evolución en N** por (ω,salida); alarma si crece más que polinómico | diagnóstico+alarma |
| acuerdo S: cruda/Ruiz/física | una discrepancia revela **inconsistencia** entre representaciones; por sí sola NO identifica cuál falla ni prueba que alguna sea inestable (podría ser error al transformar/recuperar amplitudes/ensamblar) | **cota combinada** `\|δs\|≲‖y‖(‖δb‖+‖δM‖‖x‖)`; **sin estimación fiable de δ𝓜,δb → diagnóstico empírico SIN umbral teórico `κ_out·ε`** | diagnóstico empírico |
| refinamiento iterativo (1 paso) | cambio pequeño **NO** prueba piso | registrar `η_back` **antes/después** + `‖δx‖_rel` + error forward de salidas (si disponible); **diagnóstico algebraico, no convergencia física** | diagnóstico |

**Rango (3 zonas, sobre matriz equilibrada):** *resuelto* (`σ_min` holgado sobre umbral) / *alerta*
(se aproxima → repetir con más precisión / diagnóstico independiente) / *no resuelto* (dentro de
`C₄·n_dof·ε·σ_max`, convención documentada). La 3ª zona **NO** es prueba de singularidad matemática ni de
pérdida de unicidad del problema **continuo**.

### Nivel 3 — CONVERGENCIA (aceptación real)
**El plató de `Δ_S` NO se exige `=κ_out·ε`** (la diferencia entre mallas contiene error de
discretización+dominio+solve+normalización; puede dominar `O(dx^p)` por encima del piso algebraico).
1. Exigir **tendencia de convergencia consistente** y compatible con régimen asintótico — **NO
   decrecimiento estrictamente monótono** (se admiten pequeñas oscilaciones entre mallas).
2. Si aparece plató, **investigar** compatibilidad con errores algebraicos/sensibilidades; `κ_out·ε` =
   referencia de interpretación, no predicción.
3. **Separar `p≈0` (saturación) de `p>0` (polinómico):** `p>0` **NO es automáticamente benigno** —
   depende de que `κ_out`/amplitudes sigan estables.
4. Estimar orden de S con **≥3 mallas (preferible 4)**.
5. Separar (no intercambiables): por elemento `S_ij` (≈0 → absoluta) · Frobenius `Δ_S` · estabilidad
   entre representaciones.

`‖S†η_out S−η_in‖` (pseudo-unitaridad) = comprobación **DISTINTA**, nunca criterio de convergencia.

### Tabla C1–C5
| Criterio | Métrica | Tolerancia | Estado |
|---|---|---|---|
| C1 artefacto escala | `drop=cond_cruda/cond_equil` | descriptivo, sin pass | diagnóstico |
| C2 cond intrínseco | `cond(equil)≈C·N^p` | `p≈0` satura / `p>0` solo si sensibilidades estables | provisional |
| C3 rango | `σ_min(equil)` 3 zonas | resuelto/alerta/no-resuelto; 3ª≠singularidad | forma ya; `C₄` prov. |
| C4 amplitudes estables | `\|δs\|`+`κ_out`(abs+rel), `η_back`, `Δ_S`, acuerdo repr. | η_back calibrado; sensib. informadas+N-trend; `Δ_S` tendencia→plató interpretado | mayoría provisional |
| C5 tolerancias | — | calibrar todo tras 1er barrido (≥3–4 mallas) | antes de implementar |

**Justificable AHORA (forma/orden):** `n_dof·ε`=referencia (no cota); cotas de sensibilidad
**absolutas**; rango en **3 zonas**; aceptación de observables por **tendencia de convergencia** (no por
residuo, no por plató predicho); matching y contaminación en **variables transformadas**.
**Provisional (1er barrido):** umbrales algebraicos calibrados, `p`, `C₄`, magnitud/N-trend de `κ_out`,
plató de `Δ_S` y su interpretación, métricas antes/después del refinamiento.

## Qué NO valida B2.2

Convergencia física del BVP, corrección física de S, **pseudo-unitaridad** `S†η_out S=η_in`, ausencia
general de modos espurios, ni el eje X0. Rango pleno a N finito no demuestra buen condicionamiento al
crecer N — eso es justo lo que este estudio mide.

---
*Diseño en papel cerrado. Implementación (gate B2.2) NO autorizada: requiere fijar tolerancias y luz
verde del usuario. No modifica v1.0, gates/docs cerrados, BLOCK1 ni el gate B2.1.*
