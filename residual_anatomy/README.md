# residual_anatomy — anatomía del residuo FP (λ, w)

> **Nota de la versión pública (2026-11-17).** Este documento se escribió con todos los archivos en una sola carpeta. En el repositorio público están reorganizados: solver en `solver/`, benchmark FP en `benchmarks/fp_v1/`, cadena de validación en `methods/`, correlaciones en `correlations/`, documentos en `docs/`; `out/`, `reference/` y `residual_anatomy/` conservan su ruta. Todos los scripts se ejecutan **desde la raíz del repositorio** (p. ej. `python3 benchmarks/fp_v1/FP_check_reproduction.py`). Los SHA-256 citados aquí corresponden a la estructura original; los de la versión pública están en `docs/provenance/RELEASE_MANIFEST.txt`.

Línea de investigación **nueva y aislada**, encima del solver FP congelado v1.0.
NO modifica la entrega congelada: importa el solver raíz sin tocarlo.

> **★ ESTADO ACTUAL (2026-11-15):** la **VALIDACIÓN DEL OPERADOR Padé con G fijo está CERRADA** (alcance
> documentado). Documentos autoritativos: **`gates/GATES_REGISTRY.md`** (índice de gates) y
> **`OPERATOR_VALIDATION_CLOSURE.md`** (cierre + alcance + limitaciones). C1 (pseudo-unitaridad convergente
> de S) queda **RESERVADA** → `STABILIZED_SOLVER_RESERVATION.md`. `n_L(λ)` = fase posterior.
> *Las secciones de abajo son un LOG HISTÓRICO por orden cronológico; para el estado vigente, ver esos
> documentos, no la cabecera "PASO 0 — BLOQUEADO" (superada).*

## Objetivo
Descomponer el residuo `R(w,λ) = T_upstream/(κ_local/2π) − 1` (observable n_L, Ec.42 FP)
en contribución de **suavidad** (geométrica, sobrevive a λ→0) + **dispersión**
(UV/k⁴, muere en λ→0). Separar con dos knobs independientes: λ (dispersión) y w (suavidad).

## Dependencias congeladas importadas (SHA-256, 2026-11-15)
```
correlations.py        0e3ec0f75bcbb221cf1371186384a71defcc6db66236a24412416ae4ced7a440
fd_bdg_operator.py     083f2a822e91b14401bf1388306ddc3d67daae712ba87a06b763d107517262b8
mode_matching.py       7860da7c2ccba437d1e986f577feffa6e595f91d3f8ac99baaf3d35f68e1ac7d
FP_benchmark_table.py  4eb246a83192365c5e21796dd0d807da0ecc143135e86ea0d9411476559537e0
out/fp_benchmark.csv   47fe742af4bb4e3ee5addcf8c21b24d9ab691727251eeec3694cbac0b433d955
```
`geom()` está copiado verbatim de `FP_benchmark_table.py` (geométrico, λ-independiente).

## Estado: PASO 0 (gate) — BLOQUEADO por inconsistencia de la ocupación (causa pendiente)
`residual_gate.py`:
- **Fase A PASA**: ruta nueva (scatter_field, inyección d2) ≡ benchmark (solve_inject) en λ=1
  a 4.5e-11; `build_op_lambda(1)`≡`build_MBdG`; conservación de norma OK ∀λ.
- **Fase B BLOQUEADA**: la ocupación `n_L` con λ≠1 NO es utilizable. Para λ≠1, ω·n_L→0
  (T→0) mientras para λ=1 ω·n_L→0.0411; factor ~2.6e4 a ω=9e-5 entre λ=1 y λ=2.

**NOTA conceptual (corrección, ChatGPT acertó):** n_L SÍ puede depender de λ legítimamente
(es el efecto dispersivo que se busca medir); NO vale el argumento "λk⁴ local irrelevante en
IR ⟹ n_L λ-indep". Lo que delata ARTEFACTO (no física dispersiva) es: (1) el efecto físico
esperado es ~% y lo visto es colapso de T a 0; (2) el colapso es simétrico alrededor de λ=1
(λ=0.5 y λ=4 parecidos entre sí, lejísimos de λ=1), no monótono en λ; (3) λ=1 (único punto
sano) es justo donde √λ=1 y los factores λ de la rama `op='lambda'` desaparecen.

**Estatus honesto:** inconsistencia REPRODUCIBLE de la ocupación para λ≠1; evidencia apunta
fuerte a artefacto de la ruta de cálculo (corriente/normalización/lectura de amplitud), NO a
física; causa exacta PENDIENTE de localizar. No es "tercer bug confirmado".

**Diagnóstico (2026-11-15, pasos 1/3/5 de ChatGPT):**
- NO es bug de rama discontinuo: n_L(λ) cae CONTINUA y suave de 18.66 (λ=1) a 0.44 (λ=2) a ω fija.
- Clasificación de canales IR IDÉNTICA en λ=1 y λ=2 (mismos k, vg); solo cambian los k grandes UV.
- Pseudo-unitaridad COMPLETA de S: λ=2 BUENA (3.6e-5), λ=1 MALA (7e-3) a ω=0.02·T0 (λ=1 mal
  condicionada porque n_L~T/ω diverge → elementos S~48 que se cancelan). La S de λ≠1 es
  autoconsistente con SU corriente.
- Tensión interna: correlaciones (mismo solver, knob λ) robustas ~10-45% a λ×16; n_L colapsa 99%.
- HIPÓTESIS DE NORMALIZACIÓN **REFUTADA** (2026-11-15, paso 2 ChatGPT, verificado numéricamente):
  derivada la ecuación de continuidad de build_op_lambda → J=2t·Im(u*u'+v*v')+v(|u|²−|v|²), onda
  plana = (√λk+v)|U|²+(√λk−v)|V|² = EXACTAMENTE la J del código. Verificación: J del código se
  CONSERVA sobre el campo λ=2 (std/media 6e-4); la fórmula alternativa de ChatGPT con +v(|u|²+|v|²)
  NO se conserva (std 23%) → ChatGPT erró el signo (corriente de probabilidad vs de norma). La
  corriente/normalización del código es CORRECTA → n_L(λ) es la ocupación correcta del operador.
  Las 2 observaciones de código (doble normalización, fase de eig) no son el problema.
- **FOCO NUEVO (hipótesis, sin confirmar):** el colapso de n_L es real para ESTE operador, pero
  build_op_lambda quizá NO aísla la dispersión UV pura: para fijar c²=gn reescala t=√λ/2 (masa) y
  G=gn/√λ (ACOPLO de Bogoliubov = lo que genera la mezcla partícula-hueco = Hawking). La dispersión
  Ω(k) sale bien, pero G cambia con λ para TODOS los k (no solo UV); la ocupación depende de G y ε_k
  por separado. Si se confirma, el knob λ mezcla "dispersión UV" con "reescalado del acoplo" → afecta
  la interpretación de la robustez trans-planckiana previa (CORRELATIONS.md §7-8).
- **TEST DE MODOS HECHO (2026-11-15, 3 controles de ChatGPT) — OPERADOR NO ROTO.** |V/U|_num ==
  |V/U|_form=(ε_k+G−Ω)/(ε_k+G+Ω) a precisión de máquina (ambos λ) → extracción de modos correcta.
  IR: |V/U|→1 y vg→c=1 para λ=1 y λ=2 (k=0.001: 0.999/0.9986) → límite fonónico universal, SIN
  contaminación IR. Diferencias a k finito organizadas por k_*(λ)=2c/√λ = escala de cruce dispersiva.
  → Interpretación (ChatGPT, confirmada): "diferencias a k fijo + convergencia al bajar k = cambio
  de escala de cruce dispersiva, NO contaminación hidrodinámica".

## VEREDICTO DEL DIAGNÓSTICO: operador correcto; el knob λ no es lo que creíamos
No hay bug. Corriente correcta, S pseudo-unitaria, modos universales en IR, extracción = fórmula.
El colapso de n_L es FÍSICA REAL de esta familia de operadores. Reinterpretación clave (formulación
defendible, ChatGPT): el knob λ genera una **deformación COHERENTE que conserva la velocidad del
sonido c** — t=√λ/2 y G=gn/√λ están LIGADOS por esa exigencia (2tG=c²), NO son dos modificaciones
independientes — y que modifica CONJUNTAMENTE la dispersión, el acoplo efectivo y la escala de cruce
k_*. **n_L(λ) mide la respuesta a esa deformación COMPLETA, no a una modificación aislada del k⁴**
(los tres efectos no son separables ni sumables). En w=2 (horizonte estrecho, pleno régimen dispersivo)
la respuesta es violenta → n_L colapsa en la conversión del horizonte (k~k_*).
CAVEAT (matiz ChatGPT): la pseudo-unitaridad de S demuestra consistencia respecto a la métrica de
norma, NO que el cambio de ocupación sea libre de convención — comparar n_L(λ=1) vs n_L(λ=2) compara
DOS modelos distintos y arrastra una convención de qué se fija como base de comparación.
Para variar SÓLO la dispersión UV (c, G, densidad, flujo, perfil del horizonte fijos) haría falta un
operador NO LOCAL ε_λ(k)=√(G²+c²k²+λk⁴/4)−G (otra familia de modelos, ChatGPT §3) — cuya construcción
en fondo inhomogéneo exige especificar el ordenamiento del operador (afecta norma y dispersión cerca
del horizonte); NO se obtiene por simple sustitución k→−i∂_x.

## IMPLICACIÓN PARA EL PROGRAMA residual_anatomy
El "eje λ" NO es un eje limpio de dispersión. El plan de separar suavidad(w) de dispersión(λ) debe
reorientarse. Opciones: (a) estudiar en w GRANDE (hidrodinámico) donde la deformación λ es suave y la
separación tiene sentido; (b) usar el knob de dispersión G-fijo (operador no-local) como eje de
dispersión "limpio"; (c) redefinir el eje λ como "deformación de esta familia" y reportar n_L(λ) como
tal, sin llamarlo "dispersión UV aislada". NOTA: esto NO toca el benchmark FP v1.0 (λ=1), intacto.
Sí matiza la robustez trans-planckiana previa (CORRELATIONS §7-8): válida PARA ESTA FAMILIA, pero no
demuestra "sólo se varió la dispersión UV con todo lo demás fijo".

**Conclusión:** NO proceder al Paso 1 hasta localizar y validar. El residuo TOTAL (λ=1) sí
está detectado y congelado: R(w=2) = −4.28% vs presupuesto ~0.06%.

## GATE 0 del operador NO LOCAL — EN CURSO (2026-11-15, `gate0_pseudoherm.py`)
Primera pieza de código de la fase de implementación de NONLOCAL_SPEC.md. Ruta: cuantización de
WEYL (pseudo-espectral exacta), núcleo K(x,y)=(1/2π)∫dk e^{ik(x−y)}ε_λ((x+y)/2,k). Separa, como pidió
ChatGPT, (1) pseudo-hermiticidad de (2) consistencia/convergencia al continuo.
- **2 bugs encontrados y corregidos:** (a) faltaba la medida `dx` en el matvec (K es kernel continuo;
  matvec necesita K·dx) → daba factor 1/dx (~10) en [B]; (b) la cuadratura en k necesita Nk∝N (Nk≳N),
  si no el error de dispersión CRECE con N (converge y se degrada) — justo lo que ChatGPT advirtió.
- **[B] Consistencia PASA:** K·u/u reproduce ε_λ(k) a ~1e-5; [B2] error de dispersión converge a ~1e-7.
- **[A/C] Pseudo-herm converge O(dx²)** con N (2e-2→5.7e-3→1.5e-3 en N=150/300/600; crece con dx vía L).
  El kernel de Weyl es simétrico exacto (no contribuye); el residuo viene de la CONVECCIÓN FD en bordes
  (mismo término y comportamiento que build_op_lambda, que daba 5e-5 con N=4000). No es máquina, y es
  correcto que no lo sea.
- **Nota de coste (para decisión pseudo-espectral vs Padé):** ruta densa O(N²) en construcción y solve;
  para N grande (pseudo-herm por debajo del efecto físico) conviene evaluar Padé (operador local
  disperso). El pseudo-espectral valida el modelo; producción a N grande puede pedir Padé.
- **AUDITORÍA (gate0_audit.py, cuestiones A/B/C de ChatGPT):**
  - [A] operador inhomogéneo de Weyl CORRECTO: en λ=1 el kernel es fondo-independiente a 1.5e-16
    (ε_1=k²/2 no usa n); en fondo con horizonte reproduce el símbolo LOCAL ε_λ(k;n(x)) a 1.6e-5.
  - [C] el residuo de pseudo-herm es de BULK, no de borde: a dx≈const se mantiene (3.2e-3 ∀L=20/40/80);
    es error O(dx²) de la convección, no efecto de frontera.
  - [B] residuo = max|M†η−ηM|/max|M|, η=diag(I,−I), sobre M COMPLETO, consistencia en región central.
- **MEJORA → pseudo-herm EXACTA:** la convección `v·D1+½v'` es skew solo a O(dx²); reescrita como
  `½(v·D1+D1·v)` (D1 es skew exacto incl. bordes) → pseudo-herm = **0.00 a máquina ∀N** (verificado),
  misma física continua. El operador ahora usa esta forma; cualquier residuo futuro queda aislado al
  kernel (que ya es simétrico exacto).
- **GATE 0 SUPERADO (pendiente visto bueno de ChatGPT):** operador no local consistente (reproduce
  ε_λ(k), converge al continuo, inhomogéneo correcto) y pseudo-hermítico EXACTO. 2 bugs + 1 mejora.

## GATE 1 del operador NO LOCAL — SUPERADO (2026-11-15, `gate1_homogeneous.py`)
Validación completa en fondo homogéneo CONSISTENTE (n=1, v=0.5, V=0 → Qp=0; G=c²=1). Antes se
corrigió `bg(mode='homog')` (daba Qp=−0.125 por v=0, bug señalado por ChatGPT).
- [1+2] dispersión Ω²=c²k²+λk⁴/4 y autovalores del símbolo = vk±Ω: **max error 4.4e-15** (máquina).
- [3+4] norma |U|²−|V|²=±1 exacta; **corriente J_∞ = v_g·(|U|²−|V|²) a 1e-16** (identidad de ChatGPT).
- [5] λ=1: ε_1=k²/2 exacto → Bogoliubov estándar.
- [6] discreto: converge O(dx²) (3.35e-3→8.3e-4→1.9e-4 en N=400/800/1600); la discrepancia con el
  continuo es exactamente k̃=sin(k·dx)/dx de la convección centrada (|eig−pred(k̃)|~1e-5), k̃→k O(dx²).
- VEREDICTO: GATE 1 SUPERADO (ChatGPT OK). Dispersión, modos, normas y corrientes correctos.

## GATE 2a + DECISIÓN DE ARQUITECTURA (2026-11-15)
- **Gate 2a (regresión a nivel de operador):** no local(λ=1) → build_MBdG converge O(dx²)
  (5.87e-4→1.47e-4→3.68e-5 en N=200/400/800). A nivel de operador, el no local recupera el local.
- **Decisión (aprobada por ChatGPT con condiciones):** el operador Weyl es DENSO → scattering O(N³)
  + cola del kernel truncada en bordes = impracticable para n_L(ω). → **PADÉ** (operador racional
  LOCAL disperso) para gates 2b-4. Weyl se CONSERVA como referencia; Padé y Weyl NO son el mismo
  modelo automáticamente (equivalencia en banda a demostrar).
- **Padé [2/1] derivado y documentado** (`PADE_SPEC.md`): ε_P(s)=s/2+(λ−1)s²/(4(2G+s)), s=k².
  Ancla λ=1 exacta; sin polos en k real; precisión ~% para λ∈[0.5,2] (λ=4 necesita orden superior).
  Ordenamiento hermítico D/2+(λ−1)/4·D(2G+D)⁻¹D. Caveat: Padé_inhom≠Weyl_inhom (conmutadores ~G').
- SIGUIENTE: construir operador Padé disperso → re-pasar gates 0-1 sobre Padé → gate Padé-vs-Weyl
  inhomogéneo → regresión λ=1 de observables (gate 2b) → gates 3-4 → n_L(λ) con G fijo.

## Archivos
- `residual_gate.py` — Paso 0 del knob viejo (diagnóstico; Fase B NO fiable, ver veredicto).
- `gate0_pseudoherm.py` — Gate 0 operador NO LOCAL (Weyl): pseudo-herm + consistencia + convección skew.
- `gate0_audit.py` — auditoría Gate 0 (A/B/C de ChatGPT).
- `gate1_homogeneous.py` — Gate 1: dispersión/modos/normas/corrientes en fondo homogéneo.
- `NONLOCAL_SPEC.md` — especificación teórica (cerrada y aprobada).
- `out/residual_gate.csv` — tabla del barrido del knob viejo.
