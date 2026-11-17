# Aproximación racional (Padé) del operador no local — derivación y ordenamiento

**Estado:** derivación teórica, SIN solver todavía (tarea fijada por ChatGPT antes de codificar).
Complementa `NONLOCAL_SPEC.md`. **Weyl (denso) se CONSERVA como referencia matemática**; Padé es la
implementación DISPERSA para scattering (gates 2b–4). No son automáticamente el mismo modelo: su
equivalencia en la banda hay que DEMOSTRARLA (ver §5). Fórmulas [V] verificadas numéricamente.

## 1. Objetivo
El operador no local de Weyl (NONLOCAL_SPEC §3) es denso (kernel N×N) → scattering O(N³) + cola del
kernel truncada en los bordes. Impracticable para n_L(ω) en toda la banda. Se busca una **relación
cinética racional** ε_Padé(k) ≈ ε_λ(k) que, al promover k²→−∂², dé un operador **diferencial/racional
LOCAL** (disperso de banda), reutilizable en `mode_matching`.

## 2. Derivación del Padé [2/1]
Con c²=G (=gn) y s=k², el radicando factoriza:

    G² + c²s + (λ/4)s² = (G + s/2)² + δ s² ,   δ = (λ−1)/4.              (1)  [V]

A primer orden en δ (√(u²+δs²) ≈ u + δs²/(2u), u=G+s/2):

    ε_λ(s) = −G + √((G+s/2)² + δs²)  ≈  s/2 + (λ−1) s² / (4(2G+s))  ≡ ε_P(s)   (2)

Es un **Padé [2/1] en s=k²** (numerador grado 2 / denominador grado 1). Propiedades [V]:
- **Ancla λ=1 EXACTA:** δ=0 → ε_P(s)=s/2 (dif 0.0 a máquina), recupera Bogoliubov estándar.
- **Dispersión:** Ω_P² = ε_P(ε_P+2G) ≈ c²k²+(λ/4)k⁴ (error ≈ 2× el de ε_P).

## 3. Polos y modos espurios
Único polo en **2G+s=0 ⇒ s=k²=−2G** (k imaginario) → **fuera del eje físico k∈ℝ**. No hay polos para
s=k²≥0 → **sin modos propagantes espurios** en la banda física. (A s grande ε_P ~ (λ+1)s/4, mientras
ε_λ ~ √λ·s/2; el Padé [2/1] no captura la cola UV lejana, pero esa región es k≫banda.)

## 4. AUDITORÍA (gate_pade_audit.py, 4 pruebas de ChatGPT) y decisión de orden
**Tolerancia fijada:** error de dispersión Ω² < 1% en la banda física (hasta k_*≈1.6, donde ocurre la
conversión), para no contaminar observables al nivel del % que se busca.

- **[2] Positividad/invertibilidad de 2G(x)+D:** min autovalor 2.0 (homog), 0.80 (horizonte) > 0 →
  positivo definido, invertible. ✓
- **[3] Espectro del cinético Padé:** sin ramas complejas espurias (max|Im|=0), Re≥0, homog y horizonte. ✓
- **[4] Padé vs Weyl inhomogéneo (acción):** funciones suaves difieren ~3-4e-3 en el horizonte (λ=2),
  comparable a la línea base de discretización FD-vs-espectral (λ=1: ~1e-3); localizada en el horizonte
  ~1.2e-2. El ordenamiento NO introduce discrepancia catastrófica (el conmutador §5 no domina); está al
  nivel del error de la aproximación.
- **[1] Precisión — DECIDE EL ORDEN.** Ω² al 1% se mantiene hasta:

| λ | **[2/1]** (1er orden δ) | **[4/3]** (2º orden δ) | ancla λ=1 |
|---|---|---|---|
| 0.5 | k=1.63 | toda la banda | exacta |
| 2.0 | k=1.23 (✗<k_*) | **k=2.18 (✓)** | exacta |
| 4.0 | k=0.76 | k=1.13 (✗<k_*, δ grande) | exacta |

**DECISIÓN:** la aproximación de 1er orden en δ es insuficiente para λ≥2 (1% solo hasta k<k_*). Se usa
la **aproximación racional de 2º orden en δ** (truncada, error O(δ³) — NO es un Padé [4/3] estándar en el
sentido de cociente único de polinomios grado 4/3; es una SUMA de términos racionales de la expansión):

    eps_P2(s) = s/2 + (λ−1)s²/(4(2G+s)) − (λ−1)²s⁴/(16(2G+s)³) + O(δ³)   [V]      (2b)

Derivación (√(u²+δs²)=u+½δs²/u−⅛δ²s⁴/u³+…, u=G+s/2): el ½ y el ⅛ son necesarios (su omisión en un paso
intermedio da factores 2 y 4 mal; el resultado (2b) es el correcto, verificado numéricamente).
Cumple el 1% hasta k≈2.18 para λ∈[0.5,2] (cubre k_*), con ancla λ=1 exacta. **λ=4 (δ=0.75) aún no
lo cumple → queda como extensión con orden mayor o Padé óptimo en banda.** El programa central
(¿persiste el colapso de n_L con G fijo?) se ataca en λ∈[0.5,2], donde el knob viejo ya mostraba el
colapso (18.66→0.44).

Ordenamiento hermítico del término de 2º orden: `−(λ−1)²/16·D²(2G(x)+D)⁻³D²` (forma D²A³D², A=(2G+D)⁻¹
hermítico). Implementación dispersa: 3 resoluciones de banda de (2G+D). Coste a MEDIR (§7).

**AUDITORÍA FASE 2 (gate_pade_audit2.py, pruebas 2-4 de ChatGPT) — PASA:**
- (2) Ω_P² vs c²k²+λk⁴/4 DIRECTO: 1% hasta k_* para λ∈[0.5,2] (λ=0.5 maxerr 7.6e-4; λ=2 3.3e-3 a k≤1.6).
- (3) `D²(2G+D)⁻³D²`: directa vs 3 band-solves coinciden 1e-13→1e-9 (N=200→800), adjunto ~1e-9
  (hermítico). Caveat: error crece con N por condicionamiento de (2G+D)³, vigilar en producción.
- (4) Padé2 vs Weyl, separando malla (λ=1) de estructural: **error estructural ~0.5%** (suave y
  localizado), <1%; en el caso localizado domina la discretización FD-vs-espectral, no el estructural.
  → despreciable frente al colapso buscado (~99%), marginal para efectos finos del %.
**Nomenclatura:** (2b) es aprox. racional de 2º orden en δ (O(δ³)), NO un Padé [4/3] estándar.

## 5. Ordenamiento inhomogéneo (hermítico) — y el caveat de ChatGPT
Con G=G(x)=gn(x), D=−∂². El término s²/(2G+s) promovido debe ser HERMÍTICO. La forma
`D²(2G(x)+D)⁻¹` NO lo es (D² y (2G+D)⁻¹ no conmutan para G(x) variable). La forma **simétrica** sí:

    ε_P_op = D/2 + (λ−1)/4 · D (2G(x)+D)⁻¹ D                              (3)

Hermítica: (D A D)† = D A† D con A=(2G+D)⁻¹ hermítico (2G(x)+D es hermítico: G real mult. + D hermítico).
En homogéneo (3) = (2) exacto. 

**CAVEAT CENTRAL (ChatGPT):** (3) — Padé con coeficientes x — **NO es idéntico al Weyl del símbolo
inhomogéneo** ε_λ(x,k): difieren por conmutadores [D,G(x)] ~ G'(x), **relevantes en el horizonte**
(gradiente grande). La equivalencia en la banda **debe demostrarse numéricamente** comparando (3) con
el operador de Weyl en el fondo con horizonte — NO basta con que coincidan las dispersiones homogéneas.

## 6. Criterios de validación del Padé (ChatGPT) — antes de usarlo para temperaturas
1. **Dispersión homogénea:** error de ε_P(k) y Ω_P(k) en TODA la banda de k, no solo IR (tabla §4).
2. **Polos:** inspeccionar todos los polos, descartar espurios en la banda física (§3: único polo en k im.).
3. **Pseudo-hermiticidad DISCRETA** del operador (3) + bloques BdG (no solo el símbolo homogéneo).
4. **Convergencia:** variar orden de Padé, malla (N) y caja (L).
5. **Regresión λ=1:** recuperar el operador local Bogoliubov y los observables del benchmark congelado.
6. **Comparación Padé vs Weyl en fondo INHOMOGÉNEO** (no solo dispersiones homogéneas) — acotar la
   diferencia por conmutador del §5 en la banda.
Además **re-pasar gates 0–1 sobre el operador Padé concreto** (no heredan de Weyl).

## 7. Implementación dispersa y coste
(3) se aplica resolviendo el sistema de banda (2G(x)+D)φ = Dψ (tridiagonal) y luego Dφ → todo DISPERSO.
En `mode_matching` se incorpora con una variable auxiliar φ (sistema aumentado disperso). **No dar por
hecho O(N) (ChatGPT):** medir el coste REAL por frecuencia (depende del matching, los bloques y las BC)
antes de comprometer estimaciones.

## 8. Plan
(a) elegir orden de Padé por §4/§6; (b) construir (3) disperso; (c) re-pasar gates 0–1 sobre Padé;
(d) gate adicional §5/§6.6 (Padé vs Weyl inhomogéneo); (e) regresión λ=1 de observables (gate 2b);
(f) gates 3–4 y medir n_L(λ) con G fijo. Weyl permanece como referencia cruzada en todo el proceso.

## 9. Realización auxiliar en caja: dominios, Green y condiciones de borde (Paso 2)

**ETIQUETA DE ALCANCE:** esta sección es el **CIERRE ANALÍTICO de la realización auxiliar en caja,
bajo hipótesis de borde**. NO valida el solver de scattering. La validez para física (scattering) la
decide el **Paso 3** (pseudo-herm discreta del BdG completo, prueba `A³` sobre fuentes, regresión λ=1).
Derivado en papel y verificado numéricamente (símbolos [V]/[N]); revisión cruzada con ChatGPT (paralelo).

### 9.1 K_P y la realización secuencial
Operador cinético Padé 2º orden (D=−∂², A=(2G(x)+D)⁻¹, c₁=(λ−1)/4, c₂=−(λ−1)²/16):

    K_P = ½D + c₁·D A D − |c₂|·D² A³ D²  =  ½D + c₁ D A D + c₂ D² A³ D²

Se realiza con **auxiliares esclavas** resueltas en banda (tridiagonal), NO invirtiendo A densa:

    (2G+D) p  = D u ;   (2G+D) q₁ = D²u ;   (2G+D) q₂ = q₁ ;   (2G+D) q₃ = q₂
    K_P u = ½ D u + c₁ D p + c₂ D² q₃

**Realización (a), elegida:** A³ = tres solves SECUENCIALES de (2G+D), cada uno con BC autoadjunta.
Entonces A³ es potencia de operador autoadjunto acotado → autoadjunto automáticamente.

### 9.2 Dominio y autoadjunción de L, A, A³  [cerrado, Robin real]
L = 2G(x)−∂² en [x_L,x_R], dominio D(L)={φ∈H²: (∂+κ_R)φ(x_R)=0, (∂−κ_L)φ(x_L)=0}, κ=√(2G) real
(G const en los extremos). (i) **Simétrico**: 𝒲(ψ,φ)=[ψ̄'φ−ψ̄φ']=0 bajo Robin real. (ii) **Positivo**:
⟨φ,Lφ⟩=∫(2G|φ|²+|φ'|²)+κ_R|φ(x_R)|²+κ_L|φ(x_L)|² ≥ 2G_min‖φ‖²>0.
El **signo +κ** (modo decreciente) es el que preserva positividad e invertibilidad. ⟹ A=L⁻¹ y A³
acotados, autoadjuntos, positivos.
**[N] Verificación del mínimo, con el perfil IDENTIFICADO:**
- Perfil **homogéneo G=1** (caja test): min autoval (2G+D) = **2.003** (= 2G=2 + D≥0). Es el valor que
  figuraba antes sin etiquetar — corresponde a ESTE test, no al perfil completo.
- Perfil **tanh real** (w=2, horizon): n_d=(jcur/(2√g))^{2/3}=**0.397** (NO 0.25) → 2G_min=2n_d=**0.794**;
  min autoval (2G(x)+D) medido = **0.800** (>2G_min por la contribución de D). >0 ⟹ invertible también
  en el perfil con horizonte. La cota 2G_min‖φ‖² usa el 2G_min del perfil concreto (0.794 aquí).

### 9.3 Identidades de Green de DAD y D²A³D²  [derivadas]
Con p=ADu, q₃=A³D²u, y usando que L,L³ son autoadjuntos sobre las auxiliares (∈D(L)), los términos
INTERNOS (cruces entre q₁,q₂,q₃) se cancelan y quedan SOLO los Wronskianos de las D externas:

    ⟨w,K_P u⟩−⟨K_P w,u⟩ = ½𝒲(w,u) + c₁[𝒲(w,p_u)+𝒲(p_w,u)] + c₂[𝒲₂(w,q₃ᵤ)+𝒲₂(q₃_w,u)]
    𝒲(a,b)  = [ā'b − āb']
    𝒲₂(a,b) = [ā b''' − ā'b'' + ā''b' − ā'''b]   (forma de borde de D²=∂⁴)

A autoadjunto NO implica DAD, D²A³D² autoadjuntos por sí solo: las D externas generan términos de
borde propios que hay que anular con BC sobre el campo físico.
**ALCANCE (importante, ChatGPT en P3.2-B2):** esta identidad es la de **CAJA COMPLETA** — los términos
internos (concomitantes de L y L³ entre auxiliares) se anulan por la BC natural en los EXTREMOS. Sobre
una **región estrictamente interior** NO se anulan (A, A³ son no locales y conectan dentro/fuera): la
forma correcta de sub-región añade `−c₁𝒲(p_w,p_u)` y `−c₂[𝒲(q1_w,q3_u)+𝒲(q2_w,q2_u)+𝒲(q3_w,q1_u)]`
(derivado y verificado a máquina en `gates/gate_P3_2B2_green.py`; ~1.3% de la señal para campos
genéricos). La forma de arriba es el límite R→caja.
**[N] Signo verificado (los segundos términos van con + ):** test directo (operador vs Wronskianos de
borde, campos reales que NO cumplen BC): LHS=−0.7629 vs 𝒲(w,p_u)+𝒲(p_w,u)=−0.7653 (dif 2e−3 = error
FD) vs la variante con − = −1.879 (descartada). El paso intermedio correcto es ⟨Dp_w,u⟩=⟨p_w,Du⟩
**−**𝒲(p_w,u) (no +); con + saldría el signo opuesto. Consistente con la SUMA de R_c₂ (§9.5).

### 9.4 Condiciones de borde y cancelación (hipótesis G'=0 en los extremos)  [verificado]
Reduciendo 𝒲₂ con la EDO de la auxiliar (q₃''=2G q₃−q₂) y Robin en las auxiliares:

    𝒲₂(u,q₃)|borde = −q₃·(∂+κ)(u''+κ²u) + q₂·(∂+κ)u

Robin de 1er orden sobre u NO basta (deja q₃·(−κu''−u''')). BC **suficiente** = DOS condiciones por
extremo sobre el campo físico (coherente con D²=∂⁴, 2 BC/extremo):

    borde derecho:  (∂+κ)u = 0  Y  (∂+κ)u'' = 0        [N sympy: ⟹ 𝒲₂=0]
    borde izquierdo:(∂−κ)u = 0  Y  (∂−κ)u'' = 0

Las de 1er orden `(∂±κ)u=0` anulan ½𝒲 y los de DAD; `D²A³D²` añade `(∂±κ)u''=0`. Las BC Robin de las
auxiliares (internas) y las de orden alto sobre u (externas) son independientes y compatibles con los
solves secuenciales — se imponen al DOMINIO, no a una solución particular.

### 9.5 Residuo para el perfil tanh finito y convergencia con L  [cuantificado]
Si G'≠0 en los extremos (caso real), la forma de borde COMPLETA del término c₂ deja un residuo
(ambas contribuciones de Green, NO solo una):

    R_c₂ = 2 c₂ G'(x_b) · ( w̄(x_b) q₃ᵤ(x_b) + q̄₃_w(x_b) u(x_b) )     [N sympy: 𝒲₂|borde = 2G'ū q₃ por contribución]

Para el fondo tanh (w=2, horizon; G(x)=n(x)), G'(x_b)=Dn·sech²(x_b/w)/(2w) ~ e^{−L/w}. Medido [N]
(dx≈0.08 fijo), residuo de borde 2G'ū q₃ (gaussiana localizada / onda plana k=0.8):

    | L  | G'(x_R)   | res_gauss | res_plane |
    |----|-----------|-----------|-----------|
    | 10 | −2.7e−05  | 1.7e−18   | 1.3e−06   |
    | 20 | −1.2e−09  | 1.5e−54   | 5.7e−11   |
    | 30 | −5.6e−14  | 1.5e−110  | 2.6e−15   |
    | 40 | <prec. num.| <prec. num.| <prec. num.|   (underflow; «compatible con cero numérico», NO G'=0 exacto)

Converge exponencialmente: localizados por doble supresión (G' y u(borde) decaen); scattering sigue a
G'~e^{−L/w}. **Solución adoptada:** bordes lejos (L≳20–30 para w=2) + cuantificar el residuo + exigir
convergencia en L. No se derivan BC corregidas en G'.

### 9.6 Simetría en caja ≠ autoadjunción ≠ scattering físico  [precisiones obligatorias]
κ=√(2G) es la tasa evanescente de la AUXILIAR, no la del campo físico (κ'≠κ, fijada por la dispersión).

**(a) Simetría, NO autoadjunción (aún).** Las BC de §9.4 definen un **dominio sobre el que la forma de
borde se anula y el operador es SIMÉTRICO**. La **autoadjunción** requiere ADEMÁS **maximalidad del
dominio** (D(K_P)=D(K_P\*)); que haya 2 BC/extremo es coherente con 4º orden pero NO demuestra
maximalidad para esta composición no local. **Maximalidad = PENDIENTE.** Lo demostrado es la anulación
de la forma de borde / simetría, no la autoadjunción.

**(b) El argumento de Fourier solo aplica al HOMOGÉNEO.** Hay que separar tres niveles:
- **Región homogénea (G const):** ε_P(−i∂) es multiplicador de Fourier de símbolo real → autoadjunto
  por Fourier. Argumento válido SOLO aquí.
- **Operador INHOMOGÉNEO (G(x)):** `D(2G(x)+D)⁻¹D` **NO es** multiplicador de Fourier; su autoadjunción
  NO se sigue del símbolo homogéneo. Requiere un argumento de **dominios** para la composición con G(x)
  variable. **PENDIENTE.**
- **Scattering discreto:** lo decide el **gate de pseudo-herm del BdG completo** `M†η=ηM` (§9.8).

**(c) scattering ≠ dominio autoadjunto en caja.** Las BC de §9.4 NO son las que satisfacen los modos
físicos propagantes/evanescentes (κ'≠κ). La forma de borde sobre estados de scattering ES el flujo
físico (no se anula; su balance = pseudo-unitaridad de S `S†ηS=η`, gate §9.8).

### 9.7 Precisión del flujo y límites de lo demostrado
Multiplicador de Fourier con símbolo real ε_P(k): corriente sobre U e^{ikx} = ε_P'(k)|U|². [N] K_P
periódico tiene símbolo ε_P(k) a 5e−13; ε_P'(k) coincide con ε_λ'(k)=(c²k+(λ/2)k³)/(ε_λ+G):

    δ_J(k) = |ε_P'(k)−ε_λ'(k)| / |ε_λ'(k)|      (RELATIVO)
    δ_J(k_*=1.6) = 0.655 %  ;  max_{k≤k_*} δ_J = 0.655 % (en k_*, monótono creciente) < 1 %

**Configuración exacta de esta cota:** λ=2, c²=G=g·n_u=1 (HOMOGÉNEO upstream), k_* ≡ 1.6 = escala de
cruce dispersiva donde ocurre la conversión del horizonte. La cota vale **solo** para esa banda
homogénea; NO se extiende a otras bandas, otros λ, ni a los **modos inhomogéneos** (donde G=G(x)).
⟹ contribución cinética a B_outer = ε_P'(k)(|U|²+|V|²) ≈ ε_λ'(k)(|U|²+|V|²) = parte cinética de J_∞
(NONLOCAL_SPEC ec.7); el término v(|U|²−|V|²) viene del convectivo. **Límite:** el símbolo homogéneo
coincide, pero que el término de borde del operador INHOMOGÉNEO DISCRETO sea exactamente J_∞ es Paso 3.

### 9.8 Gates PENDIENTES del Paso 3 (numéricos, deciden la validez para física)
1. Pseudo-herm DISCRETA del BdG completo (u,v,p,q₁,q₂,q₃): max‖M†η−ηM‖/‖M‖, convergencia en N,L,BC.
2. Prueba de Green A³: band-solve vs núcleo analítico e^{−κ|x−y|}/2κ y su triple convolución, fuentes
   gaussiana + onda plana; convergencia en L (sin postular exponente).
3. Gate 2b (regresión λ=1): reproducir n_L(ω) y T_upstream del benchmark FP v1.0 congelado.
4. Conservación: S†ηS=η (3 columnas) + balance por región con paquete de ondas.
5. Convergencia independiente en N, L, resolución de ω; coste real por frecuencia vs denso.
Solo tras 1–5 se mide n_L(λ) con G fijo (λ∈[0.5,2]) → ¿persiste el colapso?

**Archivos de verificación (fuera del repo, scratchpad de sesión):** recuento BC (sistema 1er orden,
Jordan), Green/hermiticidad (A, A³, K_P simétricos; κ-positividad), 𝒲₂ (sympy) + δ_J, residuo G' (tanh
real + convergencia L). A migrar a `residual_anatomy/` como gates reproducibles en el Paso 3.

## Referencias internas
- `NONLOCAL_SPEC.md` — operador de Weyl (referencia), conservación, identidad de balance.
- `gate0_pseudoherm.py`, `gate1_homogeneous.py` — gates que el Padé debe re-pasar.
