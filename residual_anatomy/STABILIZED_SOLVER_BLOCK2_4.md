# B2.4 — Diseño del gate de CONVERGENCIA en X0 (CERRADO EN PAPEL)

**Fase:** solver de scattering estabilizado, bloque 2. **Fecha:** 2026-11-15.
**Estado: DISEÑO CERRADO EN PAPEL** (ChatGPT + revisión propia con datos del perfil congelado).
**NO autoriza implementación.** No toca v1.0, validación del operador, BLOCK1, B2.1, B2.2, B2.3 ni
`GATES_REGISTRY.md`. Reutiliza por import la instrumentación auditada de B2.1 y B2.2 sin modificarla.

> **Pregunta central:** ¿el BVP global presenta una **ventana de X0** en la que S es estable y las
> amplitudes de salida están bien condicionadas, **sin el crecimiento exponencial** `e^{2κX0}` que
> perjudicaba al método de transferencia (causa acotada de C1 en v1.0)?

---

## 0. Objetivo y aislamiento del eje

Determinar si `S_box(X0)` tiene una **ventana fiable** al aumentar el dominio, con **discretización fija**.
NO se busca aún demostrar que `S_box(X0)` = S física; se busca si el problema global está suficientemente
desacoplado de los bordes y si existe una ventana numéricamente fiable.

**Tres fenómenos a separar:** (1) **error de frontera** (X0 pequeño → extremos no asintóticos);
(2) **condicionamiento** (X0 grande → representación de los modos evanescentes); (3) **error de
discretización** (ya caracterizado en B2.3 — O(dx²) —, a mantener controlado, NO re-estudiar).

**Se mantienen FIJOS:** esquema de caja 2º orden · representación física de B2.2 (ℓ=1/√2) ·
clasificación de canales y normalización de B2.1 · dx.

## 1. Barrido (DECISIÓN congelada)

| Parámetro | Valor |
|---|---|
| X0 | 8, 10, 12, 14, 16, 20, 24, 28 |
| dx principal | 0.02  (N = 2·X0/dx = 100·X0 → N=800…2800) |
| Control fino dx | 0.012, en X0 ∈ {12, 24} |
| Frecuencias | 0.05 y 0.10 |

Dominio congelado `x∈[−60,60]` (8192 pts, grid≈0.0147) → X0≤28 holgadamente dentro.
**Error de frontera ya medido (|perfil−asint| en ±X0):** ~7e−4 (X0=8), ~1e−5 (12), ~2e−7 (16), ~5e−9 (20)
→ decae ~e^{−X0}; **despreciable en X0≳14**, por lo que los X0 grandes **aíslan el condicionamiento**.

**Coste (nota de implementación):** X0=28 → N=2800 → ndof≈11212; la matriz densa ~2GB y `svdvals` es
O(n³)/O(n²) memoria. Para X0 grandes, `cond`/`σ_min` pueden estimarse con un **estimador de condición**
(o σ_min por iteración inversa) en vez de SVD completa, **validado contra SVD completa a X0 moderado**.
Las sensibilidades (adjunto), η_back y δ_adj son baratas.

## 2. El evanescente (corrección respecto al planteamiento inicial)

En el BVP, la amplitud del evanescente creciente hacia el exterior está **IMPUESTA a cero**
(`a_L,grow=0`): es una **restricción**, NO una variable libre que observar crecer con X0. Por tanto
`r_grow≡0` (métrica de B2.2) **no sirve** aquí. **No** interpretar `a_grow≈0` como prueba de buen
condicionamiento: la restricción lo fuerza.

Lo que SÍ se vigila: (a) **residuo** de esa restricción; (b) contribución de los **evanescentes
permitidos** (decrecientes) en las fronteras; (c) estabilidad de sus coeficientes al cambiar la
normalización del autovector (el coeficiente de un evanescente **no** es invariante); (d) independencia
de columnas de las bases modales (`σ_min(Ξ_b)`).

## 3. Métricas congeladas (cuatro grupos) + diagnósticos específicos de X0

| Grupo | Métricas |
|---|---|
| **Frontera** | `E_borde(X0)` por ambos extremos (ver abajo) |
| **Condicionamiento** | `cond₂`, `σ_min/σ_max`, zona de rango (3 zonas de B2.2), **`R_cond(X0)`** |
| **Sensibilidad** | `κ_out,b`, `κ_out,M` (abs+rel), `η_back`, `δ_adj` (instrumentación de B2.2) |
| **Dispersión** | `Δ_S` Frobenius + por elemento (≈0 → absoluta) + amplitudes individuales |

**E_borde(X0)** (diagnóstico directo de frontera, independiente de S):
$$E_{\rm borde}(\pm X_0)=\frac{\|A_{\rm at}(\pm X_0)-A_{\rm homog}(\pm\infty)\|_2}{\|A_{\rm homog}(\pm\infty)\|_2}$$
(2-norma, **normalización fija** por `‖A_homog‖` para comparar entre X0; se conserva también el diagnóstico
directo |perfil−asint| si es más simple). Separa el error de frontera del comportamiento de S.

**R_cond(X0)** (desenredar N de X0):
$$R_{\rm cond}(X_0)=\frac{{\rm cond}(\mathcal M(X_0))}{{\rm cond}(\mathcal M(12))}\left(\frac{N(12)}{N(X_0)}\right)^{p},\quad p=0.97\ (\text{de B2.2})$$
- `R_cond≈const` → el crecimiento es compatible con el de la malla (N).
- `R_cond` creciente → **deterioro adicional asociado a X0**.
- `R_cond` decreciente → X0 mejora el condicionamiento relativo a esa referencia.
**Caveat (límite #2):** `p` se midió **refinando dx a X0 fijo** (B2.2); aquí N crece **extendiendo el
dominio** — variaciones de N físicamente distintas → `R_cond` es **comparación heurística con B2.2, NO
ley universal**.

**Pseudo-unitaridad** `‖S†η_out S−η_in‖`: diagnóstico **separado**, nunca criterio de convergencia.

## 4. Convergencia de S_box(X0)

$$\Delta_S(X_0^{(a)},X_0^{(b)})=\frac{\|S(X_0^{(a)})-S(X_0^{(b)})\|_F}{\max(1,\|S(X_0^{(b)})\|_F)}$$
+ diferencia máxima por elemento + absolutas para elementos ≈0 + evolución de cada amplitud.
**NO** pre-elegir X0=28 como referencia: primero ver si aparece **meseta**; si la hay, `S_ref` se define
de la **región estable observada**, no del X0 mayor.

**Hipótesis de cancelación (C, NO demostrada):** a dx fijo, `Δ_S(X0_a,X0_b)` resta dos S con error de
discretización similar (transición común) → el error de dx **se cancela en gran parte**, haciendo
resolubles diferencias de X0 ~1e−5 por debajo del piso absoluto (~2.8e−3 a dx=0.02). **No** garantizada
(cambiar X0 cambia BC y nº de nodos) → el **control dx=0.012** comprueba la robustez.

## 5. Interpretación del condicionamiento (D — NO dicotómico)

- **~polinómico** → evidencia favorable de que se evitó el mecanismo exponencial de C1.
- **exponencial** → desfavorable; comprobar si procede del evanescente, de la formulación o del escalado.
- **irregular** → no clasificable todavía.
Incluso si es polinómico, **comprobar que `κ_out` sigue controlado** (cond ≠ precisión de amplitudes, B2.2).

## 6. Qué NO concluye B2.4 (aunque aparezca meseta)

- Que `S_ref` sea el límite continuo exacto. · Que `S_box(X0)` = S física. · Pseudo-unitaridad física.
- `n_L(λ)`. · Ausencia de modos espurios fuera de la ventana de frecuencias estudiada.

## 7. Tres límites explícitos (del diseño congelado)

1. `E_borde` **no sustituye** a la convergencia de S.
2. La normalización por `N^p` es **comparación con B2.2**, no ley universal (ver caveat §3).
3. Un condicionamiento favorable **no demuestra** que S sea físicamente correcta.

## 8. Estructura del gate (NO implementado)

`residual_anatomy/gates/gate_B2_4_X0_convergence.py`: único archivo nuevo; importa B2.1 y B2.2 sin
modificarlos; no toca B2.3, v1.0 ni la validación del operador; barrido de X0 con dx fijo; **sin** solve
disperso ni cambio de arquitectura. Umbrales de Δ_S: **no** universales todavía (primero observar
borde / meseta / condicionamiento). Tolerancias aplicables de B2.2 a revisar antes de fijar criterios.

---
*Diseño en papel cerrado. Implementación (gate B2.4) NO autorizada: requiere luz verde del usuario.
No modifica v1.0, gates/docs cerrados, BLOCK1, B2.1/2/3.*
