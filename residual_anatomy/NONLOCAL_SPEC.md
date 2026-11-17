# Especificación teórica — operador no local de dispersión UV con acoplo fijo

**Estado:** especificación teórica, SIN código. Define el modelo antes de implementarlo
(recomendación de la sesión 2026-11-15: primero el modelo matemático, luego programar).
No toca el benchmark FP v1.0 ni `residual_gate.py`.

**Motivación.** El knob `build_op_lambda` (usado en CORRELATIONS §7–8) genera una deformación
que conserva la velocidad del sonido c pero, para lograrlo, reescala el acoplo de Bogoliubov
G=gn/√λ y la masa t=√λ/2 (ligados por 2tG=c²). Es una deformación COHERENTE pero NO aísla la
dispersión UV: modifica conjuntamente dispersión, acoplo efectivo y escala de cruce (ver
`README.md`, veredicto del diagnóstico). Este documento define un knob alternativo que varía
**solo la forma UV de la relación cinética** manteniendo **G=gn fijo**.

Las fórmulas marcadas [V] están verificadas numéricamente a precisión de máquina (2026-11-15).

---

## 1. Qué significa "variar solo la dispersión UV"

La dispersión de Bogoliubov homogénea factoriza como

    Ω²(k) = ε(k) · [ε(k) + 2G]                                    (1)

donde ε(k) es la **relación cinética** (energía de partícula libre) y G el **acoplo** (químico).
En el BEC estándar (ħ=m=1): ε(k)=k²/2, G=gn, c²=gn, dando Ω²=c²k²+k⁴/4.

"Variar solo la dispersión UV" = deformar la cola de Ω(k) a k grande **sin** tocar:
- la velocidad del sonido c (comportamiento IR: Ω→c k, ε→k²/2),
- el acoplo G=gn (la composición partícula-hueco y el gap),
- el fondo (densidad n(x), flujo v(x), perfil del horizonte, μ).

La deformación se introduce en la **relación cinética** ε(k), manteniendo G fijo:

    ε_λ(k) = −G + √( G² + c²k² + (λ/4)k⁴ )     [V]                 (2)

Sustituyendo (2) en (1) se recupera exactamente la dispersión objetivo **para todo λ**:

    Ω²(k) = c²k² + (λ/4)k⁴     [V, dif ~1e-16]                     (3)

A diferencia de `build_op_lambda`, aquí **G=gn no se reescala**: la deformación vive
íntegramente en ε(k), que es lo que representa la física cinética/UV.

---

## 2. Magnitudes invariantes (lo que NO cambia con λ)

| Magnitud | Valor | Nota |
|---|---|---|
| Fondo n(x), v(x)=j/n | perfil tanh actual | idéntico a `bg()` |
| Horizonte x_H, μ, j | actuales | idénticos |
| c(x) = √(g n(x)) | velocidad del sonido local | idéntica |
| **G(x) = g n(x)** | **acoplo de Bogoliubov** | **FIJO — diferencia clave con build_op_lambda** |
| Presión cuántica Qp del fondo | λ-indep (fondo impuesto) | ver §6 |

**Decisión de modelado (importante):** el fondo n(x),v(x) se trata como **externamente impuesto**
(V_ext ajustado para sostenerlo), y solo se modifica el operador de FLUCTUACIÓN BdG. No se
re-resuelve el estado estacionario. Esto permite una comparación limpia a fondo fijo, igual que
hace `build_op_lambda` (Qp λ-indep). Alternativa (re-resolver el fondo con ε_λ) queda descartada
para la comparación controlada.

Composición partícula-hueco en este modelo (rama norma +):

    |V/U|² = (ε_λ + G − Ω) / (ε_λ + G + Ω) = (√(G²+Ω²) − Ω)/(√(G²+Ω²) + Ω)     (4)

que depende solo de Ω/G. En el IR (k→0, Ω→0): |V/U|→1 para todo λ (límite fonónico universal,
ya verificado en el diagnóstico). La diferencia con `build_op_lambda` está en k finito: aquí G
es fijo, allí G=gn/√λ.

---

## 3. El operador

**Símbolo homogéneo** (marco del fluido, v=0):

    M_sym(k) = [[ ε_λ(k) + G ,      G      ],
                [    −G       , −(ε_λ(k)+G) ]]                     (5)

con autovalores ±Ω(k). Añadiendo el flujo v (marco lab), el operador completo sobre (u,v):

    L_± = ε_λ(−i∂_x)  ∓  i( v ∂_x + v'/2 )  +  ( G + Qp )
    M   = [[ L_+ , G ],[ −G , −L_− ]]                             (6)

donde ε_λ(−i∂_x) es el **operador no local** obtenido promoviendo (2). El término convectivo
−i(v∂+v'/2) y el diagonal G+Qp son idénticos a `build_op_lambda`; la ÚNICA diferencia es que
el cinético −t∂² (local) se sustituye por ε_λ(−i∂) (no local) y G no se reescala.

**No localidad:** ε_λ(k) es no polinómico (raíz) para λ≠1 → ε_λ(−i∂) es un operador integral
(no diferencial) en espacio real. En fondo homogéneo se define por Fourier. En fondo inhomogéneo
(c(x), G(x) variables) requiere una prescripción de ordenamiento (§6).

---

## 4. Conservación de norma y normalización de modos

La norma conservada es ρ = |u|² − |v|² (métrica η=diag(1,−1)). **Distinción crítica para un
operador no local** (ChatGPT, correcta): hay que separar DOS cosas que en el caso local coinciden
pero aquí no.

### 4a. Flujo asintótico de ondas planas (normalización de canales)

En las regiones homogéneas asintóticas el campo es superposición de ondas planas, y el flujo de
norma por modo está bien definido:

    J_∞ = ε_λ'(k)(|U|² + |V|²) + v(|U|² − |V|²)     (onda plana)      (7)
    ε_λ'(k) = (c²k + (λ/2)k³) / (ε_λ(k) + G)        [V, dif ~1e-11]   (8)

Consistencia λ=1: ε_1'(k)=k → J_∞=(k+v)|U|²+(k−v)|V|² = corriente validada de `modes_full`. [V]
Velocidad de grupo (idéntica en ambos modelos, Ω(k) es la misma):

    v_g = v + dΩ/dk = v + (c²k + (λ/2)k³)/Ω                          (9)

Esta J_∞ sirve SOLO para normalizar los canales de dispersión en los extremos (flujo unidad
|J_∞|=1, UNA sola vez, convención de `modes_full`; 1 cuanto ↔ |J_∞|=1). **NO es una corriente
local exacta en el interior inhomogéneo** — ver 4b.

> **Advertencia (obligatoria en el manuscrito):** la expresión (7) caracteriza el flujo de modos
> planos en regiones homogéneas. No debe interpretarse como una corriente local J(x) exacta en el
> fondo inhomogéneo sin derivación adicional. Un operador no local tiene, en general, un flujo que
> involucra el campo en distintos puntos (núcleo / representación pseudodiferencial).

### 4b. Conservación de norma en el interior (global + balance por región)

Para ε_λ(−i∂) no polinómico (λ≠1) NO existe una ley de continuidad local ∂_tρ+∂_xJ(x)=0 con J(x)
sencilla. La conservación se establece a nivel GLOBAL por **pseudo-hermiticidad** del operador:

    Q = ∫(|u|²−|v|²)dx conservada  ⟺  M†η = ηM  ⟺  L_± autoadjuntos.

Condiciones suficientes (que la implementación DEBE cumplir), **válidas bajo hipótesis de dominio,
condiciones de contorno y realización discreta apropiadas** (matiz C, ChatGPT — no basta
enumerarlas; ver la precaución de discretización en §4d):
- ε_λ(−i∂) autoadjunto → cuantización de **Weyl** (símbolo real, §6);
- término convectivo −i(v∂+v'/2) simetrizado → antihermítico → ∓i·(antiherm) hermítico;
- G(x)+Qp(x) reales.

### 4c. Identidad de balance por región (resuelve el "flujo interior", ChatGPT)

Para una región R con proyector P_R (χ_R=1 en R, 0 fuera) y norma regional Q_R=⟨Ψ,P_RΨ⟩_η, con
i∂_tΨ=MΨ y M pseudo-hermítico (M†η=ηM):

    dQ_R/dt = i⟨Ψ, [M,P_R] Ψ⟩_η          [derivada verificada; V, dif 1.6e-4]   (10)

El negativo, F_R = −dQ_R/dt = −i⟨Ψ,[M,P_R]Ψ⟩_η, es el **flujo neto saliente INTEGRADO** sobre la
región R, NO una densidad de flujo local sobre una superficie (matiz A, ChatGPT: la identidad da
la variación de norma regional, no una corriente puntual). Con P_R=I, [M,I]=0 → dQ/dt=0
(conservación global). La identidad **no presupone operador local ni corriente J(x) sencilla**.

Estructura no local: con M = [[K+B+A, G],[−G, −K+B−A]] (K cinético no local autoadjunto,
B=−i(v∂+v'/2), A,G multiplicación real), A y G conmutan con P_R; B da contribución de borde; K da
el conmutador de núcleo

    [K,P_R](x,y) = K(x,y)·[χ_R(y) − χ_R(x)]                                     (11)

≠0 solo cuando x,y están en lados OPUESTOS de la frontera de R. **Interpretación (propiedad de la
representación, NO trayectoria física; matiz B, ChatGPT):** el núcleo conecta puntos a ambos lados
de la frontera → hay una contribución DIRECTA entre ellos en el balance regional. Esto NO implica
una trayectoria física que salte los puntos intermedios; es una propiedad de la representación no
local del operador, que una corriente local no puede capturar.

### 4d. Los tres niveles de conservación (controles COMPLEMENTARIOS, no redundantes)

Están relacionados matemáticamente pero NO son pruebas numéricas redundantes (ChatGPT): el nivel
global (M†η=ηM) testea el operador; el nivel entre-regiones testea ADEMÁS la implementación del
proyector P_R y del conmutador no local [M,P_R]; el asintótico testea la normalización de canales y
el balance de flujo de S. Pasar uno no implica pasar los otros.

| Nivel | Condición | Gate numérico |
|---|---|---|
| Global | M†η=ηM | max‖M†η−ηM‖/‖M‖ pequeño (discreto, ver abajo) |
| Entre regiones | identidad del conmutador (10) | dQ_R/dt directo == i⟨Ψ,[M,P_R]Ψ⟩_η |
| Canales asintóticos | J_∞ (7) + pseudo-unit de S | S†ηS=η |

**Precaución (ChatGPT):** para un estado ESTACIONARIO de frecuencia definida la densidad de norma
no cambia → balance de región = 0 trivial. Para probar (10) hay que usar un **paquete de ondas /
superposición de frecuencias**, no un modo estacionario.

**Precaución de discretización (ChatGPT, precisión 1):** que el símbolo continuo de Weyl sea real
NO garantiza M†η=ηM en la caja finita; dependen del dominio, las BC y la realización discreta. Hay
que verificar la pseudo-hermiticidad DISCRETA. El valor ~5e-5 medido es **diagnóstico de la
discretización LOCAL actual** (N, L, FD, norma del residuo max‖M†η−ηM‖/‖M‖), **NO una tolerancia
universal** que deba imponerse al operador no local: la métrica, la norma del residuo y el tamaño de
matriz deben fijarse y declararse para que la comparación sea reproducible. El criterio J(x)=const
NO aplica.

---

## 5. Límite λ=1 (ancla de consistencia)

    ε_1(k) = −G + √(G² + c²k² + k⁴/4).  Con c²=G=gn:
    G² + G k² + k⁴/4 = (G + k²/2)²  →  ε_1(k) = k²/2   EXACTO   [V, dif=0]

Por tanto en λ=1 el operador no local (6) se reduce **exactamente** al operador local
`build_MBdG` validado (ε_1 es polinómico → sin ambigüedad de ordenamiento), y debe reproducir
β²=0.01008 y T_upstream/(κ̄/2π)=1.005 (w=2). **Gate de regresión obligatorio** antes de usar λ≠1.

---

## 6. El problema del ordenamiento (fondo inhomogéneo)

ε_λ(−i∂_x) está bien definido en fondo homogéneo. Con c(x),G(x) variables, el símbolo
ε_λ(k;x)=−G(x)+√(G(x)²+c(x)²k²+(λ/4)k⁴) depende de x, y promoverlo a operador exige elegir el
ordenamiento de k↔x. Distintas prescripciones difieren en términos con derivadas de c(x),G(x),
**relevantes cerca del horizonte** (donde los gradientes son grandes). Requisitos que la
prescripción DEBE cumplir:

- (a) preservar la corriente de norma conservada (7) (pseudo-hermiticidad respecto a η);
- (b) reproducir la dispersión asintótica correcta lejos del horizonte;
- (c) reducirse a `build_MBdG` en λ=1 (automático: ε_1 polinómico, sin ambigüedad).

Propiedad clave de Weyl (ChatGPT): como c²(x)=G(x), el símbolo en λ=1 es ε_1(x,k)=k²/2
**independiente de x**, luego Op^W[ε_1]=−∂_x²/2 (local) exactamente, para cualquier perfil. El
ancla de regresión λ=1 sobrevive al ordenamiento. Para λ≠1 el símbolo sí depende de x (vía c(x),
G(x)) y el ordenamiento es genuino.

**Dos rutas de implementación (decisión para la fase de código):**

1. **Pseudo-espectral / kernel de convolución.** Aplicar ε_λ(−i∂) en Fourier. Con fondo
   inhomogéneo, prescripción de Weyl (ordenamiento simétrico): núcleo K(x,x') construido por
   simetrización. Exacto en la dispersión, pero núcleo denso → coste O(N²) y BC no triviales.

2. **Aproximación local de orden finito (Padé/racional).** Aproximar √(G²+c²k²+(λ/4)k⁴) por una
   función racional de k² → operador diferencial LOCAL de orden finito (p.ej. 4–6), tratable con
   FD disperso como el actual. Error controlable por el orden de Padé; recupera el caso exacto
   en el límite de orden alto. Sacrifica exactitud UV lejana por tratabilidad. **Probablemente
   la ruta práctica** para empezar; validar el error de truncamiento contra el espectro.

La elección de ordenamiento NO es libre: debe fijarse ANTES y verificarse con los gates (§7).

---

## 7. Plan de validación (gates antes de cualquier física)

En orden, cada uno puerta dura (si falla, parar):

0. **Pseudo-hermiticidad discreta del operador NO LOCAL — PRIMER gate, ANTES de cualquier
   observable (ChatGPT).** Construir M no local discreto y medir max‖M†η−ηM‖/‖M‖ con métrica y
   norma declaradas. Objetivo: que una desviación en n_L(λ) NO pueda confundirse con un artefacto
   de discretización. Si la realización no es pseudo-hermítica al nivel exigido, NO calcular n_L.
   **No limitarse a un único tamaño de malla (ChatGPT):** estudiar la CONVERGENCIA del residuo con
   N, L y las condiciones de contorno; comprobar que disminuye al refinar Y que su comportamiento es
   compatible con el orden de convergencia esperado. Un residuo pequeño en una única discretización
   NO demuestra pseudo-hermiticidad en el límite continuo ni ausencia de errores sistemáticos (p.ej.
   de consistencia u ordenamiento, que esta prueba puede no captar).
1. **Homogéneo.** Dispersión del operador discreto = (3) a precisión de malla; modos de (5)
   reproducen (4) y el flujo asintótico (7) se normaliza bien.
2. **Regresión λ=1.** El operador no local reducido a λ=1 reproduce `build_MBdG` y el benchmark
   (β²=0.01008, T_upstream/(κ̄/2π)=1.005, w=2) dentro de tolerancia.
3. **Conservación (tres niveles, §4d; NO J(x)=const).** Para λ≠1: (i) pseudo-hermiticidad DISCRETA
   M†η=ηM; (ii) identidad de balance por región (10) con PAQUETE DE ONDAS (dQ_R/dt directo ==
   i⟨Ψ,[M,P_R]Ψ⟩_η); (iii) pseudo-unitaridad S†ηS=η de la matriz S completa. El criterio local
   J(x)=const NO aplica al operador no local.
4. **Control de modelo (opcional, Padé).** Convergencia del observable con el orden de Padé.

### Tres requisitos (ChatGPT) — formulación teórica CERRADA; queda verificar la discretización

| Requisito | Qué queda demostrado | Estado |
|---|---|---|
| Flujo asintótico | Fórmula (7) de J_∞ y relación con v_g, para normalizar canales (§4a) | cerrado (teórico) [V] |
| Norma global | Pseudo-hermiticidad del BdG (Weyl + real + simetrizado) → Q conservada (§4b) | cerrado (teórico) |
| Flujo en el interior | Identidad de balance (10) dQ_R/dt=i⟨Ψ,[M,P_R]Ψ⟩_η con conmutador de núcleo (11) (§4c) | **CERRADO (teórico) [V]** |

La verificación de los tres a nivel de la DISCRETIZACIÓN concreta (y la pseudo-hermiticidad discreta)
pertenece a la fase de implementación (gates de §7), no a la formulación teórica. La identidad (10)
sirve igual para pseudo-espectral y para Padé → no obliga a decidir la ruta todavía.

Solo tras 1–3 se puede medir n_L(λ) / T_upstream(λ) con este knob y compararlo con el de
`build_op_lambda`. La comparación entre ambos knobs responde la pregunta de fondo: **¿el colapso
de n_L venía de reescalar G, o persiste con G fijo?**

---

## 8. Decisiones abiertas (requieren elección antes de programar)

1. **Ruta de implementación:** pseudo-espectral exacta vs Padé local de orden finito (§6).
2. **Ordenamiento:** Weyl simétrico por defecto; confirmar que cumple (a)–(c) de §6.
3. **Caveat de comparación entre modelos (matiz ChatGPT):** comparar n_L entre dos operadores
   distintos (build_op_lambda vs no local) arrastra una convención de qué se fija como base.
   Aquí la base es: fondo, c, G, flujo y perfil del horizonte idénticos; se compara la respuesta
   de la ocupación a cambiar ÚNICAMENTE la forma UV de ε(k). Debe declararse explícitamente al
   reportar, porque la pseudo-unitaridad no la resuelve.
4. **Alcance:** ¿barrido en w grande (hidrodinámico, deformación suave) además de w=2? El eje w
   sigue siendo el de "suavidad"; este knob redefine el eje de "dispersión".

---

## Referencias internas
- `README.md` (residual_anatomy) — diagnóstico que motivó esta espec.
- `CORRELATIONS.md §7–8` — robustez trans-planckiana con el knob build_op_lambda (a reinterpretar).
- `FP_BENCHMARK.md` — benchmark λ=1, intacto (ancla de regresión §5/§7).
