# P2 — RE-REFERENCIA DE FASE DE S · DISEÑO CONGELADO (en papel)

**Estado:** diseño **CERRADO EN PAPEL** y aprobado para congelación documental (2026-11-16) tras 4 rondas
de revisión externa (referee ChatGPT: v1 → v2 → v3 → v4 + precisiones finales). Texto = v3 consolidada
con las sustituciones de v4 y las precisiones del veredicto final.

**Esta aprobación autoriza SOLO la creación de este documento. NO autoriza implementar ni ejecutar
`gates/gate_P2_phase_reference.py`**: eso requiere autorización expresa y separada del usuario.

**Alcance:** P2 estudia **solo el solver LOCAL λ=1** (`A_at`, gates B2.1/B2.2 cerrados, importados en
solo lectura). P2 **NO** valida la S física (C1/P5), ni el operador Padé, ni n_L(λ), ni κ_out.
No modifica v1.0, gates cerrados ni documentos cerrados. κ_out sigue bloqueado hasta cerrar P2.

| Elemento | Estado |
|---|---|
| Fórmula de k_box y transformación S_ext → S_fis^box | **PROBADO algebraicamente** bajo las convenciones declaradas (§1–§2) |
| Que el re-faseo explique la deriva de B2.4 / B2.4b | **HIPÓTESIS** (P2 debe poder refutarla) |
| K1/K2/K3, estabilidad en X0 de los 9 elementos, convergencia en dx | **PENDIENTE** |
| Validez física de S_fis | **ABIERTA** (requiere P5/C1) |

---

## 0. Hipótesis bajo prueba
**H0 (HIPÓTESIS):** a dx fijo, la variación de S_ext con X0 (Δ_S O(1) en B2.4; residuo lineal ~1e−3 en
B2.4b) es **solo un re-faseo por canal** por la fase acumulada en la zona homogénea, que en la malla se
propaga con `k_box`, no con k continuo.

Predicciones: **(P-a)** módulos e invariantes de re-faseo independientes de X0 (salvo cola que decae);
**(P-b)** des-rotado con k_box → S independiente de X0; con k continuo → deriva lineal conocida;
**(P-c)** S des-rotada converge O(dx²) elemento a elemento (módulo y fase), con el mismo límite para
cualquier X0 y ambos des-rotados. **El gate se diseña para que cualquiera pueda fallar.**

## 1. Convenciones (PROBADO por inspección de `gate_B2_1_algebraic.py`, `gate_B2_2_conditioning.py`)
- Malla `xs = linspace(−X0, X0, N+1)`, dx=2X0/N. Caja: `(I − dx/2·A_{j+1})Y_{j+1} = (I + dx/2·A_j)Y_j`.
- Matching `Y_0 = Ξ_L a_L` (x=−X0), `Y_N = Ξ_R a_R` (x=+X0). Ξ = autovectores de `A_homog` (fondo
  asintótico: medias de 50 puntos extremos de [−60,60]); propagantes normalizados a |J|=1 con la
  corriente del continuo; orden por (Re k, Im k).
- Amplitudes **referidas al borde**. Canales IN={L,R1,R2}, OUT={L,R1,R2}. σ_L=−1, σ_R=+1.

### 1.1 Base modal congelada
- **Referencia por ω** calculada una vez antes de los barridos: k_m, ξ_m (fase de `la.eig`, misma
  convención que B2.x), etiquetas, orden, sha256 de Ξ_L, Ξ_R completas.
- En cada ejecución, `info["XiL"]`, `info["XiR"]` (de `g1.assemble`, NO modificado) deben ser
  **idénticas bit a bit** a la referencia, junto con etiquetas, orden y k. Cualquier diferencia → fallo I-C0.
- **Identificación de canales (solo DETECCIÓN):** min_{m≠n}|k_m−k_n| (<1e−3 → flag casi-degenerado) y
  matriz de solapamientos |ξ_ref†ξ| para confirmar correspondencia. **Nunca se reordena ni modifica Ξ**;
  una permutación o pérdida de correspondencia → fallo I-C0.
- Alineación de la convención de fase con v1.0: **ABIERTO, fuera de P2 (P5)**.

## 2. Álgebra (PROBADO algebraicamente; confirmado por el referee)
- Con Aξ=λξ, λ=ik: μ=(1+ik·dx/2)/(1−ik·dx/2)=e^{i k_box dx}, **k_box=(2/dx)·arctan(k·dx/2)** (rama
  principal). k real ⇒ |μ|=1 (Cayley); **k_box − k = −k³dx²/12 + O(k⁵dx⁴)**. Autovectores de la caja =
  autovectores de A → solo se corrige la fase acumulada, no la base.
- Amplitud referida al origen: A_m = a_m e^{−i k_box,m x_b}, x_b = σX0. Por tanto
  **S_fis^box[o,p] = e^{−i k_box,o σ_o X0} · S_ext[o,p] · e^{+i k_box,p σ_p X0}**; S_fis^cont idem con k.
- **Identidad (NO evidencia):** S_fis^cont[o,p] = S_fis^box[o,p]·exp{iX0[(k_box−k)_o σ_o − (k_box−k)_p σ_p]}
  ≡ S_fis^box·e^{iΦ_op}. Ambas salen del mismo S_ext.

## 3. Pruebas numéricas — tres convergencias separadas

### 3.1 k_box (convergencia 1)
- **K1:** fórmula → μ_K1, k_box,K1.
- **K2:** autovalores de T=(I−dx/2·A_homog)⁻¹(I+dx/2·A_homog), identificados por solapamiento. Comparar
  **μ primero**: |μ_K2−μ_K1| ≤ 1e−13. Luego k_box,K2=(arg μ + 2πn)/dx con n que deja k_box más cerca de k;
  **dx·|k_box,K2 − k_box,K1| ≤ 1e−13**. Propagantes: ||μ|−1| ≤ 1e−13. Registrar error en μ, en k_box y n.
- **Guard de Nyquist:** |k|·dx ≤ π/4 en todos los canales y mallas (estimación previa: peor caso ~0.10).
- **K3-syn (control homogéneo sintético; prueba la HERRAMIENTA; DECISIVO):** en el script de P2, sin tocar
  gates cerrados: Y_0=Ξ·a, propagar 20 celdas con el MISMO esquema de caja y A_j≡A_homog, aplicar
  **exactamente el mismo procedimiento K3** (proyección Ξ⁻¹, cociente entre celdas, filtro
  |c_m|≥1e−3·max|c|, elección de pareja). Vectores a: (i) aleatorio con semilla fija, todos los modos O(1);
  (ii) con los **cocientes de amplitud** de las columnas reales junto a cada borde (incl. modo pequeño en el
  límite del filtro). Suelo esperado ε_K3 = u·cond(Ξ)·(max|c|/|c_m|), u = épsilon de máquina.
  **Criterio:** |μ_K3syn − μ_K1| ≤ 10·ε_K3 en todos los canales/celdas/vectores (factor 10 = **margen
  empírico**; se registra medido/ε_K3).
- **K3 sobre la solución real (SIEMPRE diagnóstico, sin pass/fail):** μ_K3 = c_m(x_{j+1})/c_m(x_j),
  c=Ξ⁻¹Y_j, últimas 10 parejas junto a cada borde, columna donde el canal es mayor, filtro 1e−3.
  Registrar |μ_K3−μ_K1| por pareja, δ_b=‖A_at(x_b)−A_homog‖, **δ_10 = max_{j∈últimas 10 parejas}
  ‖A_at(x_j)−A_homog‖**, cocientes |μ_K3−μ_K1|/(dx·δ_10) y |μ_K3−μ_K1|/ε_K3; informar en qué X0 se alcanza
  el suelo de K3-syn. **Dato para A1, NO prueba independiente de H0; ninguna cota presentada.**
- **Reconstrucción:** D_out·S_fis·D_in⁻¹ = S_ext con error ≤ 1e−14.

### 3.2 Convergencia en X0 a dx fijo (convergencia 2)
dx=0.02; X0∈{8,10,12,14,16,20,24,28} (N=800…2800); ω∈{0.10,0.05}; X_ref=16.
**Métrica común:** E_X(S) = max_{i,j}|S_ij(X0) − S_ij(X_ref)|.
- **(i) Invariantes independientes de k:** E_X(|S|) (9 módulos); fases plaqueta
  I_ij = arg(S_00·S_ij·S_0j*·S_i0*), i,j∈{1,2}, sobre **S_ext cruda**; I_ij **no definida** si alguno de sus 4
  elementos tiene |S|<1e−3·‖S‖_F; E_X^I = max|wrap(ΔI_ij)| (envuelto a (−π,π]). Sin I_ij definidas → C2
  solo con módulos, y se hace constar.
- **(ii)** E_X(S_fis^box); **(iii)** E_X(S_fis^cont) (complejos).
- **(iv) Proxy heurístico de suelo numérico** ε̃_S = cond_est·η_back·max|S|. **NO es cota certificada** del
  error en S: cond_est y η_back describen la sensibilidad de la solución completa x, no de cada salida
  (P2 no reintroduce κ_out); la validación de cond_est vs SVD (B2.4c) no la convierte en cota. Registrar
  **por separado** cond_est, η_back, max|S| y ε̃_S. cond_est = estimador de B2.4 (σ_max por potencia; σ_min
  por iteración inversa con `lu_solve(trans=2)`, ruta certificada en P1), importado sin modificar.
- **(v) Diagnósticos:** η_back, residuo absoluto, ‖x‖, |a_dec|, δ_b, δ_10, verificación bit a bit de Ξ.
- **Fase vs amplitud por elemento:** |ΔS_ij|, Δlog|S_ij|, Δarg envuelto; fase interpretada solo si
  |S_ij| ≥ 1e−3·‖S‖_F. **Ningún re-faseo puede mover los módulos.**

### 3.3 Convergencia en dx (convergencia 3)
- X0=12 y X0=16, cada uno con N=500/1000/2000 (sin X0=24 de entrada; solo si 12 y 16 discrepan, como
  ampliación de alcance documentada).
- Richardson **por elemento** y **en Frobenius** (diagnósticos distintos, informados aparte), para
  S_fis^box y S_fis^cont; S_ref por X0 y por des-rotado.
- **No estimable:** |S_h−S_{h/2}| o |S_{h/2}−S_{h/4}| < 1e−10 → p no asignado.
- **Métrica de límites:** E_ref(S_a,S_b) = max_{i,j}|(S_a−S_b)_ij|, para box↔cont y X0=12↔16.
- **Conjunto material 𝓜:** elementos con |S_ij| ≥ 1e−3·‖S‖_F evaluados sobre **S_fis^box en N=2000**,
  **por separado para cada (ω, X0)**. Los no materiales se excluyen de la afirmación de orden y se listan
  por nombre. Un elemento de 𝓜 sin p estimable **no se descarta**: deja C5 **inconcluso**.

## 4. Criterios fijados antes de ejecutar

**Integridad del procedimiento — si falla cualquiera, el gate se DETIENE sin clasificación A/B/C/D:**
| # | Comprobación | Éxito |
|---|---|---|
| I-C0 | Sistema, base y solución | ventana 3×3; Ξ idéntica bit a bit a la referencia, sin permutaciones detectadas; guard de Nyquist; η_back=‖Mx−b‖/(‖M‖‖x‖+‖b‖) ≤ 1e−12 (residuo absoluto registrado aparte) |
| I-C1 | Identificación de la caja | \|μ_K2−μ_K1\| ≤ 1e−13; dx·\|k_box,K2−k_box,K1\| ≤ 1e−13; \|\|μ\|−1\| ≤ 1e−13 en propagantes; K3-syn ≤ 10·ε_K3 |
| I-C4 | Identidad entre des-rotados | \|S_cont − S_box·e^{iΦ_op}\| ≤ 1e−12·max(1,\|S_ij\|) en los 9 elementos y todas las ejecuciones (pendiente lineal solo diagnóstico) |

**Pruebas científicas — solo si pasa toda la integridad:**
| # | Comprobación | Éxito | Fallo → conclusión |
|---|---|---|---|
| C2 | Invariantes | E_X(\|S\|) ≤ 1e−6 y E_X^I ≤ 1e−6 para X0≥16; para X0=8…14 E_X decrece monótonamente (tasa informativa) | la explicación por re-faseo puro no describe los datos dentro de estos umbrales y este dominio → diagnosticar causa (cola, BC, base, solver) sin concluir nada sobre la física |
| C3 | S_fis^box | **Decisivo:** E_X(S_box) ≤ 1e−6 para X0≥16. **Indicativo (no decide):** E_28(S_cont)/E_28(S_box), informado solo si E_28(S_box) ≥ 10·ε̃_S,pair con **ε̃_S,pair = max{ε̃_S(28), ε̃_S(16)}**; si no, "no concluyente" | ver clasificación B |
| C5 | dx | todo elemento de 𝓜 con p estimable y 1.8≤p≤2.2 (por elemento; Frobenius aparte); E_ref(box,cont) ≤ 5e−5 y E_ref(S_ref(12),S_ref(16)) ≤ 5e−5 | **falla:** p fuera de rango o límites no coincidentes. **Inconcluso:** algún elemento de 𝓜 sin p estimable |

**Clasificación (solo con integridad PASS), fijada de antemano:**
- **A:** C2, C3 y C5 pasan → H0 PROBADA **dentro del alcance documentado**.
- **B:** C2 pasa, C3 falla → la estabilidad de los invariantes es compatible con un re-faseo, pero la
  des-rotación propuesta no estabiliza S_fis^box dentro de los criterios. Se investigan convención de fases,
  identificación de canales, BC y solver. **No se concluye automáticamente que k_box o σ estén mal.**
- **C:** C2 falla → diagnosticar causa.
- **D:** C2 y C3 pasan, C5 falla.  **D\*:** C2 y C3 pasan, C5 inconcluso.

Todo resultado (incluida la detención por integridad) se persiste.

## 5. Explicaciones alternativas y firma discriminante
| Alternativa | Firma |
|---|---|
| A1 Cola del perfil (A_at(x_b)≠A_homog, incl. v′ e interpolación) | afecta módulos e invariantes; decae ~e^{−X0} (2/w=1); correlacionada con δ_b, δ_10 |
| A2 Inestabilidad de la base | bloqueada por I-C0 (congelación bit a bit, detección de permutaciones) |
| A3 Orden de canales / σ mal asignados | C2 pasa, C3 falla; pendiente residual encaja con otra combinación de k |
| A4 Evanescente que decae no despreciable | deriva no de fase con tasa \|Im k_dec\|≠1; ver \|a_dec\|(X0) |
| A5 Normalización con corriente del continuo | error O(dx²) en módulos independiente de X0; afecta C5, no C2; tratamiento en C1/P5 |
| A6 Suelo numérico | η_back, ε̃_S, ‖x‖; residuo pequeño ≠ exactitud |

## 6. Trazabilidad (Control B)
- `gates/gate_P2_phase_reference.py` (NUEVO, **no autorizado aún**). Importa g1, g2 y el estimador de B2.4
  en solo lectura; ruta primal; `lu_solve(trans=2)` solo dentro del estimador de σ_min. **Sin κ_out.**
- CSV: `out/p2_modal_ref.csv` (base congelada); `out/p2_kbox_check.csv` (K1/K2/K3 en μ y k_box, rama n,
  k·dx/π, δ_b, δ_10, filas K3-syn: vector a, canal, celda, μ, ε_K3, cond(Ξ), medido/ε_K3);
  `out/p2_X0_sweep.csv` (S_ext/S_box/S_cont Re/Im, Δmódulo, Δfase, flags, η_back, residuo abs, cond_est,
  max|S|, ε̃_S, |a_dec|); `out/p2_invariants.csv` (módulos, I_ij y flags, E_X);
  `out/p2_dx_sweep.csv` (todos los S, p por elemento y Frobenius, 𝓜 y excluidos, S_ref, E_ref).
- Manifiesto `.json`: sha256 de scripts, `waterfall_state.npz` y Ξ_ref; entorno (Python, NumPy, SciPy,
  OpenBLAS, plataforma); parámetros; umbrales; comando exacto; código de salida; tiempos; pico de memoria;
  **estado de integridad (PASS / STOP + motivo)** y **clasificación A/B/C/D/D\***.

## 7. Alcance
- Barrido X0: dx=0.02, N=2X0/dx ∈ [800, 2800]. Barrido dx: X0∈{12,16}, N≤2000. Toda extensión posterior se
  documenta como **ampliación del alcance**.
- ω∈{0.05, 0.10}; operador **local λ=1**.
- P2 **no** valida la S física (C1/P5), el Padé, n_L(λ) ni κ_out.

## 8. Historial de revisión (referee externo ChatGPT)
- **v1:** álgebra y signos aprobados; pedidos: rama del log en K2, umbral de K3 por extremo, congelar base
  modal, invariantes con condición de definición, métrica común C3/C4, p no estimable; X0=24 fuera de entrada;
  "si C2 falla" ≠ refutación física.
- **v2:** pedidos: δ_10 y K3 sin cota arbitraria; C4 como identidad elemento a elemento; cociente C3
  protegido del suelo; métrica E_ref; η_back de B2.2; coherencia de alcance; solapamiento solo detecta.
- **v3:** pedidos: ε_S como proxy heurístico; K3 diagnóstico salvo control homogéneo sintético; cobertura
  mínima en C5; integridad (C0/C1/C4) separada y B más prudente.
- **v4:** aprobación científica; precisiones finales: ε̃_S,pair = max de las dos ejecuciones;
  tolerancia dx·|Δk_box| ≤ 1e−13; 𝓜 definido sobre S_fis^box en N=2000 por (ω,X0).
- **Veredicto final:** *P2 aprobado para congelación documental. La aprobación autoriza la creación de este
  documento, no la implementación ni la ejecución del gate.*

*Congelado 2026-11-16. Siguiente acción (pendiente de autorización expresa del usuario): implementar
`gates/gate_P2_phase_reference.py` conforme a este diseño, auditarlo, y solo después ejecutarlo.*
