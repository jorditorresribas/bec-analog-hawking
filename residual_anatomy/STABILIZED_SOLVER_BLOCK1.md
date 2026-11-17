# BLOQUE 1 — Condiciones modales en los extremos (diseño en papel, CERRADO)

**Fase:** diseño del **solver de scattering estabilizado** (ver `STABILIZED_SOLVER_RESERVATION.md`).
**Fecha de cierre en papel:** 2026-11-15. **Estado: APROBADO EN PAPEL** (revisión cruzada ChatGPT +
verificación numérica propia). **NO autoriza implementación.** El bloque 2 (primer gate algebraico con N
pequeño) permanece CERRADO hasta autorización expresa del usuario.

> **Qué es este documento:** la derivación, revisada y corregida, de las **condiciones de contorno
> modales** del BVP global que sustituirá al transfer ingenuo de v1.0 (`bdg_smatrix.py`), que es la causa
> acotada del floor ~5e−5 / blow-up en X0 de C1. **No** es el solver. **No** toca v1.0 ni los gates/docs
> ya cerrados. Todo lo etiquetado PROBADO está verificado numéricamente o es lectura directa del código
> congelado; HIPÓTESIS = decisión de diseño a revisar en su gate; ABIERTO = pendiente.

---

## 1.0 Convención de estado (PROBADO — lectura de `bdg_smatrix.py`)

Vector de 1er orden **Y = (u, u′, v, v′)ᵀ**, con `Y′ = A(x) Y`:

```
A = [[0,      1,       0,     0      ],
     [a_uu,  −2i V_f,  2gn,   0      ],
     [0,      0,       0,     1      ],
     [2gn,    0,       b_vv,  2i V_f ]]
a_uu = V_f² + 4gn + 2V − 2μ − 2ω ;  b_vv = V_f² + 4gn + 2V − 2μ + 2ω
```

`V_f` = velocidad de flujo (≠ amplitud de Bogoliubov v). **Decisión (HIPÓTESIS fijada):** se adopta este
orden (u, u′, v, v′), no (u, v, u′, v′); es permutación trivial, pero fijarla evita bugs de reordenación.

Parámetros del fondo congelado (`out/waterfall_state.npz`): g=1, μ=1.125, j=0.5; upstream (izq)
**subsónico** (Mach≈0.5), downstream (der) **supersónico** (Mach≈2.0).

---

## 1.1 Base modal homogénea (PROBADO)

Dispersión verificada (residuo del cuártico ~1e−12):

> **Ω² = c² k² + ¼ k⁴,  Ω ≡ ω − V_f k,  c² = gn**
> → `k⁴/4 + (c²−V_f²)k² + 2ωV_f k − ω² = 0` : 4 raíces k_m por región y por ω.

**Autovector analítico** (de la ecuación de Bogoliubov uniforme; cruza con `eig(A)`):

> **r_m ≡ v_m/u_m = (Ω_m − ½k_m² − gn) / (gn)**
> **ξ_m = u_m · (1, i k_m, r_m, i k_m r_m)ᵀ**

Útil como control cruzado cerca de umbrales, donde `eig` numérica degrada.

---

## 1.2 Clasificación — por ω, con las cuatro guardas (PROBADO vs `bdg_channels.py`)

Para **cada ω** (nunca presuponer conteo):

| Magnitud | Fórmula | Uso |
|---|---|---|
| Propagante | `\|Im k\| < tol·(1+\|Re k\|)` | si no → evanescente |
| Norma de Krein | `sign(Ω)`, Ω = ω − V_f k | signo de η por canal |
| **Dirección** | `v_g = V_f + (c²k + k³/2)/Ω` → `sign(v_g)` | in/out (**NO** sign(k)) |
| Evanescente | `Im(k)` | izq: `Im k > 0` ⇒ crece en x→−∞ ⇒ **excluir** |

Asignación in/out (confirmada en el código):
- **Izquierda x=−X0 (subsónica):** IN ↔ v_g>0 (hacia el horizonte); OUT ↔ v_g<0.
- **Derecha x=+X0 (supersónica):** OUT ↔ v_g>0 (se aleja); IN ↔ v_g<0.

**ABIERTO:** que la región supersónica tenga 0 evanescentes es cierto *dentro de la ventana de Hawking*,
pero **no se presupone**: se comprueba ω a ω; si aparece un evanescente creciente hacia +∞, se le impone
amplitud 0 y el conteo se re-evalúa (ver invariante en §1.4bis.D).

---

## 1.3 Dirección y acotación del evanescente excluido (PROBADO)

El modo excluido en la izquierda es el **acotado-violante cuando x→−∞**: `e^{ikx}` con `Im k > 0` da
`e^{−(Im k)x} → +∞` al ir a −∞. Simétricamente, cualquier evanescente creciente hacia +∞ en la derecha se
excluiría. "Creciente" se entiende SIEMPRE como no acotado hacia el infinito del extremo correspondiente.

---

## 1.4 Matching por amplitudes aumentadas, SIN formar Ξ⁻¹ (HIPÓTESIS de diseño)

Descomposición en cada borde: `Y_b = Ξ_b · a_b`, con `Ξ_b = [ξ_1 … ξ_4]` (4×4) y columnas **ordenadas por
el clasificador**:

- Izq: `Ξ_L = [ξ_in | ξ_out | ξ_dec | ξ_grow]`, `a_L = (a_in, a_out, a_dec, a_grow)`.
- Der (en ventana): `Ξ_R = [ξ_in1 | ξ_in2 | ξ_out1 | ξ_out2]`, `a_R = (a_in1, a_in2, a_out1, a_out2)`.

Las amplitudes `a_L, a_R ∈ ℂ⁴` entran como **incógnitas** y el matching como **filas lineales**
`Y_b − Ξ_b a_b = 0`; `Ξ_b` es un bloque 4×4 denso que **nunca se invierte**. QR/SVD se usan solo para
(a) base modal ortonormalizada y bien condicionada, y (b) detectar coalescencia (`σ_min(Ξ_b)`).

**Matiz honesto:** el matching aumentado **no reduce** `cond(Ξ_b)`; elimina la amplificación artificial de
fabricar filas de BC con una inversa explícita (`R⁻¹Qᴴ`) y organiza el conteo. En la **coalescencia modal
genuina** el matching es mal-puesto para cualquier formulación → esa ω es **frontera de ventana**, se
marca, no se fuerza. Coincide con la divergencia de la normalización |J|=1 (§1.5) y con `σ_min(Ξ_b)→0`.

---

## 1.4bis Ensamblado y unicidad del sistema aumentado (PROBADO — recuento; condición de unicidad)

### A. Interior — esquema de caja (trapezoidal), recuento limpio

Nodos `x_0=−X0 … x_N=+X0` (N+1 nodos), `Y_j ∈ ℂ⁴`. Una ecuación vectorial por intervalo:

> `(I − (dx/2)A_{j+1}) Y_{j+1} − (I + (dx/2)A_j) Y_j = 0`,  j=0…N−1

→ **4N ecuaciones** interiores; **4(N+1)** incógnitas de campo; déficit = 4 → exactamente 4 BC.
(La discretización concreta es HIPÓTESIS del bloque 2; la caja se elige por el recuento exacto.)

### B. Filas de matching (8) y C. restricciones (4)

Matching: `Y_0 − Ξ_L a_L = 0` (4), `Y_N − Ξ_R a_R = 0` (4).
Restricciones (selectores canónicos `e_m^⊤`):

> `e_4^⊤ a_L = 0` (grow izq anulado) · `e_1^⊤ a_L = s_L` (in izq) · `e_1^⊤ a_R = s_{R1}` · `e_2^⊤ a_R = s_{R2}`

Para los 3 problemas de scattering, `(s_L, s_{R1}, s_{R2})` recorre `e_1, e_2, e_3` de ℂ³.

### D. Recuento → CUADRADO

| | Ecuaciones | Incógnitas |
|---|---|---|
| Interior (caja) | 4N | — |
| Campo Y_j | — | 4(N+1) |
| Matching izq+der | 8 | — |
| Amplitudes a_L,a_R | — | 8 |
| Restricciones | 4 | — |
| **Total** | **4N+12** | **4N+12** |

**Invariante de conteo (debe cumplirlo el clasificador por ω):**

> **#(incoming, ambos extremos) + #(evanescentes crecientes hacia el exterior, ambos extremos) = 4**

En ventana: izq {1 in + 1 grow} = 2, der {2 in + 0} = 2 → 4 ✓ (verificado: izq 2 propagantes + 2
evanescentes; der 4 propagantes + 0). Si no suma 4 → umbral/mal-puesto → **bandera**, reconstruir, no
parchear.

### E. Extracción de S

Las salientes son **incógnitas resueltas**: `a_out^L=(a_L)_2`, `a_out1^R=(a_R)_3`, `a_out2^R=(a_R)_4`.
Resolviendo con `(s_L,s_{R1},s_{R2}) = e_1,e_2,e_3` → las 3 columnas de **S ∈ ℂ³ˣ³** directamente, sin
proyección posterior. Requiere columnas propagantes normalizadas a |J|=1. Pseudo-unitaridad luego con
`η = sign(Ω)` por canal.

### F. Rango y unicidad (condición N&S y jerarquía)

> **Unicidad del sistema discreto ⟺ ker(𝓜_aug) = {0}**  (única condición necesaria **y** suficiente)

| Nivel | Enunciado | Rol |
|---|---|---|
| **Estructural (necesaria, no suficiente)** | `Ξ_L, Ξ_R` regulares | un vector nulo de `Ξ_b` da amplitud modal que no altera el campo de borde ⇒ pérdida de unicidad del sistema aumentado tal como está formulado |
| **N&S de unicidad** | `ker(𝓜_aug) = {0}` | homogéneo solo trivial |
| **Numérica (suficiente, NO necesaria)** | `I−(dx/2)A_{j+1}` invertibles | comodidad para la recurrencia; `Y_{j+1}` aparece en DOS filas → bloque-columna `[B_{j+1};−C_{j+1}]` puede tener rango 4 aunque `B_{j+1}` sea singular ⇒ la global puede ser regular con un bloque aislado singular |

**Pérdida de regularidad modal:** la dependencia lineal de columnas de `Ξ_b` **puede** aparecer en
coalescencias de raíces —incluidos umbrales con `v_g→0`—, se detecta con `σ_min(Ξ_b)`, **sin presuponer**
que toda singularidad provenga de ese mecanismo.

**ker(𝓜_aug)={0} vs σ_min pequeño (precisión):** `ker(𝓜_aug)={0}` es **unicidad exacta** del sistema
discreto; un **valor singular mínimo pequeño** es **casi-singularidad numérica**, cuya interpretación
**depende del escalado**. No son lo mismo y no se deben intercambiar.

**Homogéneo no trivial:** indica **no-unicidad / casi-singularidad** del BVP truncado. **NO** se etiqueta
como "modo espurio" sin controles adicionales (dependencia con X0, dx, consistencia asintótica) que
distingan modo ligado real, resonancia, artefacto de caja o degeneración numérica.

**Diagnósticos — qué prueban y qué NO:**

| Diagnóstico | Qué prueba | Qué NO prueba |
|---|---|---|
| `σ_min(Ξ_b)` | independencia de la base modal | nada del sistema global |
| Pivotes LU | alerta de condicionamiento (depende de escalado/orden) | rango |
| `cond(𝓜_aug)` estimado | condicionamiento global | unicidad exacta |
| Residuo `‖𝓜_aug x − b‖` | que el solve algebraico cerró | **ni rango ni unicidad** |
| `ker(𝓜_aug)` (homogéneo) | **unicidad** del sistema discreto | clasificación física de la causa |

Regla operativa: residuo pequeño **no** implica rango completo; pivotes razonables **no** sustituyen el
análisis de núcleo; `σ_min` pequeño es casi-singularidad dependiente de escala, no `ker≠{0}`.

---

## 1.5 Corriente, velocidad de grupo y norma (PROBADO — verificado numéricamente)

Fórmula de corriente de v1.0:

> `J = Im(u*u′) + V_f|u|² + Im(v*v′) − V_f|v|²`

**Verificación numérica (18 modos propagantes: izq+der, ω=0.05/0.1/0.2; estado congelado):**

> **J = v_g (|u|² − |v|²)** — ratio J/predicción = **1.0000** en las 18 filas.
> **sign(J) = sign(v_g)·sign(Ω)** — se cumple en las 18 filas.
> **sign(J) ≠ sign(Ω)** en general — p. ej. der, ω=0.05, k=+2.13: sign(J)=+1, sign(Ω)=−1 (modo de norma negativa, canal Hawking).

Relaciones cerradas:

> **ρ_K ≡ |u|² − |v|²  (densidad de norma de Bogoliubov),  sign(ρ_K) = sign(Ω)**
> **J = v_g · ρ_K  ⇒  sign(J) = sign(v_g) · sign(Ω)**
> **η_canal = sign(Ω)**  (métrica de Krein, **NO** sign(J));  **|J|=1** solo normaliza el flujo.

**Nota:** v1.0 ya guarda `sign(J)` y `norm=sign(Ω)` **por separado** (correcto). El error de conflación
`η=sign(J)=sign(Ω)` fue solo de una redacción previa de este bloque, ya corregido; la referencia v1.0 es
consistente.

---

## Qué queda para el bloque 2 (NO autorizado)

Primer gate **algebraico con N pequeño**, SIN observables Hawking: ensamblado de 𝓜_aug, cuadratura
(4N+12), `σ_min(Ξ_b)`, `ker(𝓜_aug)` y consistencia del matching. Después: discretización interior
definitiva, derivadas en los extremos, **ventana de X0** (diagnóstico que detecta "demasiado corto = no
asintótico" y "demasiado largo = subdesbordamiento del evan decreciente ~e^{−2κX0}"), y previsión de la
ruta aumentada (p,q1,q2,q3 + BC propias) para λ≠1.

**Pregunta física central (ABIERTA):** con G fijo y λ∈[0.5,2], ¿persiste el colapso de n_L(λ)?

---
*No modifica v1.0 ni gates/docs cerrados. Revisión cruzada ChatGPT (2026-11-15) + verificación numérica
propia de §1.5. Bloque 2 cerrado hasta autorización.*
