# P3.4 — DISEÑO (papel) conservación de la matriz S (pseudo-unitaridad) — NO CONGELAR, NO IMPLEMENTADO

**Estado:** DISEÑO en papel. NO hay implementación; NO toca v1.0 ni los gates aprobados (B2, B2c, B3).
Redactado a petición de ChatGPT (orden: fijar convenciones de canales IN/OUT, métrica η, signos de norma
y orientación de S ANTES de los residuos y controles de flujo). Marcas: **[PROBADO]** aquí · **[HIPÓTESIS]**
a verificar en el gate · **[ABIERTO]** sin respuesta.

**Separación esencial (ChatGPT):** B3 fue un diagnóstico de pseudo-herm PROYECTADA del OPERADOR; NO
garantiza conservación de la norma. **P3.4 la demuestra independientemente** vía la matriz S y el balance
de flujo. P3.4 no reinterpreta B3 ni se apoya en él.

---

## P3.4.0 — Espacio de canales y convenciones (se fijan ANTES de nada)

**Geometría de scattering:** fuera de `±X0` el fondo waterfall es homogéneo (n,v,V constantes a cada lado,
L≠R). Los estados asintóticos son modos plano `e^{ikx}(U,V)ᵀ` del operador BdG homogéneo de cada lado.

**Raíces y clasificación** — SEPARAR λ=1 de Padé (corrección obligatoria ChatGPT):
- **λ=1 (BdG local):** la dispersión es `(ω−vk)² = c²k² + k⁴/4` — el término `k⁴/4` es la **presión
  cuántica** (longitud de curación), SIEMPRE presente (coef. 0.25 en `quartic_roots` de v1.0). Es
  **cuártica en k, 4 raíces** contando multiplicidades (NO `c²k²`, que sería el límite hidrodinámico
  cuadrático). Se clasifican por propagación/evanescencia, `vg` y norma.
- **λ≠1 (Padé):** el operador es una **aproximación RACIONAL** a la dispersión. Su ecuación espectral NO es
  automáticamente cuártica: al eliminar denominadores puede dar un **polinomio de grado superior con raíces
  adicionales, polos y modos auxiliares** (p,q1,q2,q3 del sistema aumentado). La dispersión IDEAL que se
  aproxima es `(ω−vk)² = c²k² + λk⁴/4` (cuártica, vg abajo), pero las raíces REALES del solver salen de la
  **ecuación espectral COMPLETA del sistema aumentado**. **NO se clasifica una raíz Padé como canal físico
  sin comprobar que es una solución PROPAGANTE del sistema aumentado y NO un modo espurio** (auxiliar o de
  polo). Esta verificación es parte de los gates de BC del solver aumentado (P3.4.5), previa a P3.4 λ≠1.
- `Ω = ω − vk` (frecuencia en el marco del fluido). Propagante si `Im(k)≈0`; si no, evanescente
  (decae/crece). Velocidad de grupo de la **dispersión ideal** `vg = v + (c²k + λk³/2)/Ω` (vale para λ=1 y
  como referencia). **Para el operador Padé, `vg` se obtiene de la RAMA ESPECTRAL REAL del sistema
  aumentado** (`vg = dω/dk` de su ecuación espectral completa), NO se reutiliza automáticamente la derivada
  de la dispersión ideal (pueden diferir donde la aproximación racional se desvía).
- **Norma de Bogoliubov del canal:** `η = sign(Ω)` [PROBADO v1.0]. (Es el signo de la norma de Krein
  `|U|²−|V|²` del modo; para modos propagantes coincide con sign(Ω).)
- **Corriente de norma asintótica** (conservada, para normalizar canales):
  `J = Im(U*U') + v|U|² + Im(V*V') − v|V|²` [forma v1.0, `norm_current`]. Cada canal se normaliza a `|J|=1`.
  Para el operador NO LOCAL (Padé) `J` es el FLUJO ASINTÓTICO de ondas plano (NO corriente local exacta en
  el interior; ver NONLOCAL_SPEC §4a). En los tramos homogéneos ±X0 es exacta.

**Convención IN/OUT** (por signo de `vg` respecto al horizonte en x≈0) [PROBADO v1.0 `smatrix`]:
- **IN (entrantes hacia el horizonte):** `Lin` (izquierda, vg>0), `Rin1,Rin2` (derecha, vg<0) → **3 IN**.
- **OUT (salientes):** `Lout` (izquierda, vg<0), `Rout1,Rout2` (derecha, vg>0) → **3 OUT**.
- Evanescentes: el que CRECE en −∞ se excluye (BC física); el que decae entra en el matching. (El recuento
  3+3 corresponde a la ventana de ω con 2 canales propagantes a la derecha; si ω sale de esa ventana el
  recuento cambia y el gate debe re-clasificar, no asumir 3×3 — **[HIPÓTESIS]** a controlar por ω.)

**Métricas de norma de canal (SEPARADAS in/out):**
`η_in = diag(sign Ω)` sobre los 3 IN ; `η_out = diag(sign Ω)` sobre los 3 OUT. En general `η_in ≠ η_out`
(distinto reparto de signos) → la pseudo-unitaridad NO es `S†ηS=η` con una sola η (ver P3.4.1).

---

## P3.4.1 — Matriz S: orientación y pseudo-unitaridad

**Orientación:** `a_out = S · a_in`. La **columna j de S = inyectar el canal IN_j con amplitud 1** (resto
0) y leer las amplitudes de los 3 OUT. "3 columnas" = los 3 canales IN. S es 3×3. [PROBADO v1.0 `smatrix`]

**Pseudo-unitaridad (conservación de la norma de Krein):** la conservación de `Q = Σ η|a|²` entre IN y OUT
da **`S† η_out S = η_in`** (NO `S†ηS=η`: in y out tienen su propia métrica de signo). [PROBADO: es la forma
verificada en v1.0, `S.conj().T @ Eout @ S - Ein`.] **Residuo:** `R_S = S† η_out S − η_in` (matriz 3×3).

Interpretación de los bloques: la diagonal = conservación de norma por canal; fuera de diagonal = ausencia
de mezcla espuria de norma. El coeficiente de Hawking `|β|²` sale de las entradas norma+ ↔ norma− de S.

**Registro de canales (convención FIJA, obligatorio — ChatGPT):** el orden de los 3 IN y 3 OUT se fija
como convención y se REGISTRA, por canal, antes de construir S (evita que una permutación o una raíz mal
clasificada produzca una S plausible con métrica incorrecta):

| Dato | Definición |
|---|---|
| Lado | izquierda / derecha |
| Dirección | IN / OUT |
| k | raíz asintótica (real si propagante) |
| vg | signo de la velocidad de grupo |
| Ω | frecuencia comóvil ω−vk |
| Norma | sign(Ω) (= η del canal) |
| Corriente | signo y módulo de J ANTES de normalizar |

La normalización `|J|=1` se aplica de forma CONSISTENTE a IN y OUT. El gate imprime esta tabla para cada ω.

---

## P3.4.2 — Controles (DOS independientes + un diagnóstico)

**(C1) Pseudo-unitaridad MATRICIAL de S:** `max|R_S| = max|S† η_out S − η_in|`, **incluyendo términos
diagonales Y no diagonales** (la conservación por canal Y la ausencia de mezcla espuria). Canales
asintóticos estacionarios.

**(C2) Balance de flujo por región con PAQUETES** (INDEPENDIENTE de S; no usa modos estacionarios).
**Identidad GENERAL (corrección ChatGPT — incluye el término de defecto de pseudo-herm, porque B3 NO
demostró `M†𝒢=𝒢M` exacto):**

    dQ_R/dt = i⟨Ψ, [M, P_R] Ψ⟩_η  +  i⟨Ψ, (M†𝒢 − 𝒢M) P_R Ψ⟩
              └─ flujo regional por la frontera ─┘   └─ término de DEFECTO de pseudo-herm ─┘
    P_R = proyector sobre R=[a,b];  [K,P_R](x,y)=K(x,y)[χ_R(y)−χ_R(x)] (≠0 solo a través de la frontera).

El 2º término se anula SOLO si la pseudo-herm es exacta; como no lo es (B3: defecto O(dx²) físico, O(1)
rejilla), **se MIDE y se REPORTA** (no se asume despreciable). Convención de producto interno `⟨·,·⟩_𝒢` y
representación matricial fijadas explícitamente (las de B3: 𝒢=diag(W,−W)).

**Protocolo operacional del paquete (obligatorio, ChatGPT) — se comprueba el balance DURANTE una
evolución, no una identidad de conmutador aislada:**
  1. Paquete inicial `Ψ(0)` y canal en que se prepara (p.ej. Lin, band-limitado baja k, como B3).
  2. Región `R=[a,b]` y su frontera (a,b interiores, lejos de ±X0).
  3. Intervalo temporal de evolución `[0,T]` (el paquete cruza el horizonte dentro de R).
  4. Variación medida `ΔQ_R = Q_R(T) − Q_R(0)`.
  5. Integral temporal del flujo regional `∫₀ᵀ (término frontera + término de defecto) dt`.
  6. **Balance:** `|ΔQ_R − ∫flujo| / escala` → 0 al refinar discretización (dx, dt) y variar el ancho del
     paquete. Reportar por separado la magnitud del término de defecto.

**(C3) Diagnóstico por canal (DERIVADO de C1, NO independiente):** `Σ_out η_out|S_ij|² = η_in,j` es
EXACTAMENTE la componente diagonal (j,j) de `S†η_out S = η_in`. Se reporta como desglose canal a canal
(balance de flujo por columna), **pero NO cuenta como evidencia matemática independiente de C1** (ChatGPT).

---

## P3.4.3 — Ensayos de convergencia

- **En dx del integrador** (RK4 del transfer, o paso del solver aumentado): `max|R_S|(dx) → 0` con el orden
  del método (RK4 → esperado alto; el límite real lo pone el matching modal y la clasificación de canales).
- **En X0** (semi-ancho de la región integrada): `R_S` → plateau (fondo homogéneo fuera de ±X0; si no se
  estabiliza, X0 es insuficiente o hay cola evanescente sin resolver).
- **En ω**: barrer ω dentro de la ventana física de 3 canales; NO extrapolar fuera de ella. Reportar
  `|β|²(ω)` y comparar con el benchmark congelado en λ=1 (enlaza con P3.3/Gate 2b).
- **Para el paquete (C2):** convergencia en resolución espacial Y temporal (dx, dt) y en el ancho del
  paquete (espectro concentrado a baja k, como en B3) → el balance regional (con el término de defecto
  incluido) debe cerrar con orden definido; reportar aparte la magnitud del término de defecto vs el de
  frontera.

---

## P3.4.4 — Criterios de aceptación / rechazo (fijados ANTES de ejecutar)

**ACEPTAR P3.4** si TODO:
- (A1) `max|R_S| = max|S†η_out S − η_in|` **medido en ≥4 resoluciones** y **→0 con pendiente registrada**
  frente a dx (el 4.36e−5 a ω=0.1 es solo un punto de referencia, NO prueba convergencia), con **meseta en
  X0**. Se reporta valor + pendiente, no solo el umbral 1e−3.
- (A2) Balance de flujo con paquete (C2): `|ΔQ_R − ∫flujo| / escala → 0` con orden definido DURANTE la
  evolución; **el término de defecto de pseudo-herm se incluye y se reporta** (no se asume cero).
- (A3) **Variar X0 y comprobar que los coeficientes de S alcanzan una MESETA**; convergencia en dx; y
  **verificar que el nº de canales propagantes NO cambia** dentro de la ventana de ω examinada (si cambia →
  re-clasificar, no forzar 3×3).
- (A4) En λ=1: `|β|²(ω)` y `T_upstream` reproducen el benchmark congelado (= P3.3/Gate 2b). Es un **ancla
  SEPARADA** de que el solver+clasificación son correctos; **NO sustituye a C1 ni a C2** (ChatGPT).

**RECHAZAR / PARAR** si:
- (R1) `R_S` no converge o queda O(1) → canal mal clasificado, BC mal puestas (evanescente que crece no
  excluido), o solver roto. Revisar antes de interpretar `|β|²`.
- (R2) El balance de flujo (C2) no cierra → la identidad regional no se cumple en la discretización →
  revisar el conmutador/BC; NO reclamar conservación.
- (R3) `η_in` o `η_out` dependen de detalles de la clasificación (modos casi degenerados, |Ω|≈0) → la
  asignación de norma es ambigua en esa ω; excluir esa ω o refinar la clasificación.

---

## P3.4.5 — Alcance y DEPENDENCIAS (crítico)

- **Baseline λ=1 (local):** el solver de transferencia de v1.0 (`bdg_smatrix.py`) ya produce S y pasa C1
  (**verificado 2026-11-15: ω=0.1 → `max|R_S| = 4.36e−5` ≪1e−3; η_in=[1,1,−1], η_out=[1,−1,1] DISTINTAS**
  → confirma que la forma es `S†η_out S = η_in`, no `S†ηS=η`). P3.4 formaliza sus controles y convergencia.
  **Enlaza con P3.3/Gate 2b** (regresión λ=1 vs benchmark). Se valida PRIMERO aquí (no toca v1.0: ancla).
- **λ≠1 (Padé, el objetivo físico):** requiere el **SCATTERING SOLVER del operador Padé** = sistema
  AUMENTADO disperso (auxiliares p,q1,q2,q3 de P3.2-A). Ese solver es un **SUB-PROYECTO APARTE** con sus
  PROPIOS gates de condiciones de contorno (orden de ChatGPT, ya registrado): (1) derivar el sistema
  aumentado desde la expresión racional; (2) justificar físicamente las BC de p,q1,q2,q3; (3) homogéneo
  reproduce raíces/modos; (4) sin modos espurios ni dependencia artificial de BC; (5) solo entonces
  scattering. **P3.4 para λ≠1 DEPENDE de ese solver** — no se puede cerrar la conservación de S dispersiva
  antes de construirlo y pasar sus gates de BC.
- **Realización en BANDA:** tarea independiente (antes de escalar a N grande), con sus propios controles.
- **NO toca v1.0 ni B3.** B3 (pseudo-herm proyectada del operador) NO implica C1/C2; la conservación se
  demuestra aquí, no se hereda.

---

## Orden
P3.4 (este) → P3.3/Gate 2b (regresión λ=1 vs benchmark congelado, usa el solver λ=1) → [construir solver
Padé aumentado con sus gates de BC] → P3.4 para λ≠1 + medir `n_L(λ)`, G fijo, λ∈[0.5,2]: **¿persiste el
colapso?** (pregunta física, ABIERTA).

**Pendiente antes de IMPLEMENTAR P3.4:** visto bueno de ChatGPT a este diseño (en especial P3.4.0
convenciones y P3.4.1 la forma `S†η_out S = η_in`, y los controles C1–C3 / criterios A1–A4). No arrancar
código hasta entonces.
