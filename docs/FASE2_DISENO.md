# Fase 2 — Diseño: una dispersión con dos escalas que sí aísla la física UV

**Estado:** diseño en papel + criterios de los gates técnicos F0–F3, escritos y en commit **antes** de implementarlos.
La medición física (exponente de la anchura del horizonte ensanchado) **no se ejecuta** hasta que el revisor externo
haya visto este diseño; su predicción y sus criterios están en `residual_anatomy/FASE2_PREREG.md`.
Rama `fase2` desde `main` = `3d128a1`. No se modifica v1.0 ni ningún gate existente: todo el código nuevo está en
`residual_anatomy/fase2/`. Fecha: 2026-11-18.

## 1. Por qué λ·k⁴ no servía y qué hace falta
Con una sola escala dispersiva, la dispersión (ω − vk)² = c²k² + λk⁴/4 entra en la física solo a través del número
adimensional κξ/c (o w/ξ). La fase 1 lo demostró numéricamente: (λ, w) ≡ (1, w/√λ) hasta 3,5·10⁻⁹.

En general, con (ω − vk)² = c²k² + A·k⁴ + B·k⁶, el reescalado x → x/s, ω → s·ω (que lleva una anchura w a w/s y deja c,
v y T/T_H invariantes) transforma A → A/s² y B → B/s⁴. El cociente

  **β̂ ≡ B/A²** es invariante bajo cualquier reescalado.

El modelo físico (λ = 1) tiene A = 1/4 fijo; toda la familia (λ, w) de la fase 1 tiene B = 0, es decir β̂ = 0. Un modelo
con B ≠ 0 tiene β̂ ≠ 0 y **no puede ser equivalente a ningún (λ, w′)**: hay un segundo número adimensional, que mide la
**forma** de la dispersión a k grande con la escala fija. Eso es lo que se quiere aislar.

## 2. El modelo: interacción de alcance finito
### 2.1 Por qué no una cinética de cuarto orden
Lo inmediato sería una energía cinética ε(k) = k²/2 + a·k⁴. Pero entonces ε(ε + 2c²) = c²k² + (¼ + 2ac²)k⁴ + a·k⁶ + a²·k⁸:
aparece k⁸ y el coeficiente de k⁴ cambia con c²(x). Y en un flujo, una cinética ∂⁴ rompe la invariancia galileana: el
espectro ya no es simplemente (ω − vk)². Se descarta.

### 2.2 La elección: U(k) = g(1 + γk⁴)
Con la cinética estándar y una interacción de dos cuerpos de alcance finito con transformada U(k) = g(1 + γk⁴), la
relación de Bogoliubov en un flujo uniforme es

  (ω − vk)² = ε_k(ε_k + 2nU(k)) = **c²k² + k⁴/4 + γc²·k⁶**,  con ε_k = k²/2 y c² = gn.

Es exactamente la forma pedida, con **β = γc²**. Como c depende de x, β también; el parámetro del modelo es γ, y se
usa β_H = γc_H² en el horizonte. β̂ = B/A² = 16γc², distinto de 0 para todo γ > 0. Estabilidad: U(k) > 0 para γ ≥ 0.

En el espacio real, U*f = g(f + γ f⁗) (porque (ik)⁴ = k⁴): es un operador **local** de cuarto orden, discretizable con el
esténcil (1, −4, 6, −4, 1)/dx⁴ (segundo orden).

### 2.3 Ecuaciones
GP: i∂_tψ = [−½∂² + V + (U*|ψ|²)]ψ. Con ψ₀ = √n e^{iφ}, v = φ′, las ecuaciones BdG en las variables comóviles (u, v) son
las de v1.0 más un bloque nuevo:

  M_γ = M_v1.0 + γ·g·[[P, P], [−P, −P]],  P = diag(√n)·D₄·diag(√n),

donde M_v1.0 = `fd_bdg_operator.build_MBdG(x, n, v, V₀, μ)` con el potencial de v1.0 (ver 2.4). El término de Hartree
γ·g·n⁗ de la diagonal se cancela exactamente con el que hay que añadir al potencial externo para que el fondo sea
estacionario (2.4), y por eso la diagonal efectiva es la de v1.0.

### 2.4 Fondo estacionario (la lección de la fase 1)
Se mantiene el mismo perfil n(x) en tanh, v = j/n (Mach 0,5 → 2, anchura w) y se usa el método inverso: el potencial
externo que hace estacionario ese perfil con la interacción nueva es

  V_γ = μ − ½v² − g(n + γ·n⁗) + ½(√n)″/√n = V₀ − γ·g·n⁗.

Es decir: **el fondo con γ ≠ 0 necesita un potencial externo distinto**; si se reutilizase V₀ (el error de la fase 1, en
otra forma), el fondo dejaría de ser estacionario. Con V_γ, el modo de Goldstone Φ_G = (√n, −√n) cumple
M_γΦ_G = M_v1.0Φ_G, porque el bloque nuevo actúa sobre √n·(u + v) = 0. La métrica acústica (c, v), κ y T_H son las de v1.0:
solo cambia la física UV.

### 2.5 Raíces y *matching* modal
Para ω fija, en cada región asintótica homogénea, k cumple el polinomio de grado 6

  −γc²·k⁶ − ¼k⁴ + (v² − c²)k² − 2ωv·k + ω² = 0  (grado 4 si γ = 0: el de v1.0).

- **Clasificación:** raíz real = propagante, con Ω = ω − vk, norma = signo(Ω), v_g = v + (c²k + k³/2 + 3γc²k⁵)/Ω.
  Compleja = evanescente; se acepta la que decae hacia fuera (Im k < 0 a la izquierda, > 0 a la derecha).
- **Recuento esperado en la ventana de Hawking:** a la izquierda (subsónica) 2 reales + 4 complejas → 1 entrante
  (v_g > 0), 1 saliente y 2 evanescentes aceptables (3 «away»); a la derecha (supersónica) 4 reales + 2 complejas →
  2 entrantes (una de norma −), 2 salientes y 1 evanescente (3 «away»). Se comprueba ω a ω, como el I4 de E1.
- **Modo plano:** (U, V) del símbolo [[½k² + vk + G_k, G_k], [−G_k, −(½k² − vk + G_k)]], con G_k = c²(1 + γk⁴),
  normalizado a |J| = 1 con J = v_g(|U|² − |V|²) (flujo = velocidad de grupo × densidad de norma, válido para ondas
  planas de cualquier sistema conservativo invariante por traslaciones; con γ = 0 debe coincidir con el J de v1.0:
  lo comprueba F2).
- **Ensamblado (grado 6).** Las filas u y v de M_γ tienen un esténcil de ±2 puntos. Su suma (fila u + fila v) solo
  alcanza ±1, porque el bloque P se cancela. Por eso el sistema discreto tiene orden 6 (como el continuo), y se ensambla:
  filas u y v en i = 2…N−3, la fila suma en i = 1 y en i = N−2, y el *matching* de (u, v) con la superposición de modos
  en p = 0, 1, 2 y p = N−1, N−2, N−3. Recuento: 2(N−4) + 2 + 12 = 2N + 6 ecuaciones = 2N incógnitas de campo + 3 + 3
  amplitudes. Con γ = 0 se usa exactamente la disposición de v1.0 (matching en 2 puntos por lado, 2 + 2 amplitudes).
- **Corrección de implementación (2026-11-18, antes de ejecutar ningún gate):** el *matching* de (u, v) en p = 0, 1, 2
  descrito arriba deja el sistema **casi singular**. El estado del sistema discreto de orden 6 no es «u, v en 3 nodos»,
  sino s = u + v en 4 nodos y d = u − v en 2 (la fila u alcanza ±2 en s y ±1 en d). Con u, v fijados en 0, 1, 2, la fila
  suma en i = 1 solo usa nodos ya fijados y no aporta información. Síntoma: en un flujo uniforme supersónico con γ = 1,
  inyectar un modo daba «transmisión» 1,59 y defecto 0,27. Se cambia a: *matching* de s en p = 0…3 y de d en p = 0, 1 por
  lado (mismo recuento, 2N + 6). Con eso, la transmisión en el homogéneo es exacta (defecto 10⁻⁹–10⁻¹¹). Además, la fase
  de cada modo se refiere a su borde (x_ref = x[0] o x[−1]) en lugar de a x = 0, para no escalar las columnas evanescentes
  por e^{±Im k·L}. Los criterios de los gates no cambian.
- **Solve:** equilibrado + refinamiento iterativo (`diag_conditioning.solve(..., 'eqref')`), porque el operador de cuarto
  orden empeora el condicionamiento.

## 3. La anchura del horizonte ensanchado para k⁶ (argumento de la fase)
Coutant y Parentani (arXiv:1402.2514, Ec. 7) deducen para dispersión cuártica d_br = (2κΛ²)^(−1/3) a partir de la fase
de los modos en la representación de momento cerca del horizonte. Reproduzco el argumento para una potencia general.

Cerca del horizonte, c − v ≈ κx (en el sistema de la onda que remonta). Una relación de dispersión
Ω = cp·√(1 + F(p)), con F pequeño, da Ω ≈ cp + ½cp·F(p). En la representación p, la posición del modo de frecuencia ω
es x(p) = [ω − Ω_disp(p)]/(κp) + …, y la fase ∫x(p)dp contiene un término dispersivo:

- **cuártico**, F = p²/Λ²: Ω_disp = c·p³/(2Λ²) → fase ∝ c·p³/(6κΛ²). La dispersión importa cuando ese término es de
  orden 1: p_* ~ (κΛ²/c)^{1/3}. La extensión espacial correspondiente es δx ~ c·p_*²/(κΛ²) ~ (κΛ²)^{−1/3} (con c = 1),
  que es el d_br de Coutant y Parentani salvo el factor numérico.
- **séxtico**, término β·k⁶ en Ω²: Ω ≈ cp + β·p⁵/(2c) → fase ∝ β·p⁵/(10cκ). Entonces p_* ~ (cκ/β)^{1/5} y
  δx ~ β·p_*⁴/(cκ) ~ **(β/(cκ))^{1/5}**.

Generalizando, para un término dominante k^{2n}: d ∝ κ^{−1/(2n−1)} (n = 2 → 1/3; n = 3 → 1/5).
**Esto es una deducción por análisis de escalas, no una fórmula publicada** (no hay ninguna para k⁶; ver la auditoría de
literatura). El prefactor no se predice.

**Régimen.** El término k⁶ domina sobre k⁴/4 cuando 4γc²p² ≫ 1. Con p ~ 1/d_br ~ 1 en las anchuras ensayadas
(w = 2…8, d_br ≈ 1), hace falta γc_H² ≳ 2 (γ ≳ 3) para estar en el régimen séxtico; con γ ≲ 0,1 se está en el cuártico,
y entre medias, en el cruce.

## 4. Qué se mide (la predicción, pre-registrada aparte)
Para cada w ∈ {2, 3, 4, 6, 8} y cada γ: T_up = T_0(n_L) (el observable de v1.0), y la semianchura d*(w) de la ventana
simétrica tal que κ̄(d*)/2π = T_up (la ventana que hace exacta la fórmula de Finazzi-Parentani). Ajuste
d* ∝ κ^{−p} en las cinco anchuras (κ ∝ 1/w a Machs fijos). Predicción: p = 1/3 con γ = 0 (control) y p → 1/5 con k⁶
dominante. Tolerancias, éxito y fracaso: `residual_anatomy/FASE2_PREREG.md`.

## 5. Gates técnicos (criterios fijados antes de implementar)
Malla: perfil tanh, L = 40, N = 240·L, salvo indicación. Si un gate falla, ese camino se detiene y se documenta.

| Gate | Qué comprueba | Criterio |
|---|---|---|
| **F0 — fondo y Goldstone** | (a) con V_γ, ‖M_γΦ_G‖ coincide con el residuo de v1.0 (γ = 1 y 10); (b) el residuo es O(dx²) (N = 2000, 4000, 8000; interior \|x\| ≤ L − 10); (c) control negativo: con V₀ (fondo NO consistente) el residuo es ≥ 100 veces mayor y no baja con dx | (a) dif. rel. ≤ 10⁻¹²; (b) orden observado en [1,8; 2,2]; (c) cociente ≥ 100 y orden observado < 0,5 |
| **F1 — homogéneo** | flujo uniforme subsónico (Mach 0,5) y supersónico (Mach 2), γ = 1 y 10, ω = 0,05 y 0,2·ω_max: (a) cada raíz cumple el símbolo (\|det\| relativo ≤ 10⁻¹⁰) y la onda plana cumple el operador FD con residuo O(dx²); (b) norma: signo(\|U\|² − \|V\|²) = signo(Ω) y \|J\| = 1; (c) recuento de raíces (6; reales y evanescentes según §2.5); (d) inyección en el homogéneo: transmisión = 1 y reflexión ≈ 0; (e) pseudo-unitariedad en el horizonte (w = 2): \|Σ η\|a\|² − η_in\| converge con orden ≈ 2 | (a) orden en [1,8; 2,2]; (b) exacto a 10⁻¹²; (c) exacto; (d) \|T − 1\| ≤ 10⁻⁶ y \|R\|² ≤ 10⁻⁶ con N fino; (e) orden en [1,5; 2,5] y defecto ≤ 10⁻³ en la malla fina |
| **F2 — límite γ → 0** | (a) con γ = 0, el driver nuevo (raíces del polinomio general, J = v_g(\|U\|²−\|V\|²), disposición de grado 4) reproduce n_L de v1.0 (`mode_matching.solve_inject`) en w = 2, 9 ω; (b) informativo: n_L(γ) − n_L(0) frente a γ pequeño (resuelto en la malla) | (a) dif. rel. ≤ 10⁻⁸ con solve refinado; (b) se informa el orden en γ (se espera 1) |
| **F3 — no degeneración** | w = 2, γ = 10: D(s) = máx_ω \|n_L^γ(ω) − n_L^{(0, w/s)}(s·ω)\| / n_L^γ(ω), la distancia a la familia de reescalados (que contiene toda la familia (λ, w′) de la fase 1). Control positivo: el mismo cálculo con el modelo λ = 4 debe dar D_min ≈ ruido | D_min(γ = 10) ≥ 100 × ruido, con ruido = dif. máx. de n_L^γ entre N y 1,5N; y control: D_min(λ = 4) ≤ 10 × su ruido |
