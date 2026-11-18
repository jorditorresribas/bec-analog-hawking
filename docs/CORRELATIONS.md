# RUMBO #1 — Correlaciones de Hawking + robustez trans-planckiana

> **Nota de la versión pública (2026-11-17).** Este documento se escribió con todos los archivos en una sola carpeta. En el repositorio público están reorganizados: solver en `solver/`, benchmark FP en `benchmarks/fp_v1/`, cadena de validación en `methods/`, correlaciones en `correlations/`, documentos en `docs/`; `out/`, `reference/` y `residual_anatomy/` conservan su ruta. Todos los scripts se ejecutan **desde la raíz del repositorio** (p. ej. `python3 benchmarks/fp_v1/FP_check_reproduction.py`). Los SHA-256 citados aquí corresponden a la estructura original; los de la versión pública están en `docs/provenance/RELEASE_MANIFEST.txt`.

> ⚠️ **AVISO DE AUDITORÍA (2026-11-17) — §4, §5, §7 y §8 NO TIENEN VALIDEZ FÍSICA.**
> Todos los resultados con λ ≠ 1 de este documento se calcularon con `build_op_lambda`, que
> conserva la presión cuántica del fondo **Qp independiente de λ**. Para λ ≠ 1 ese fondo no es
> estacionario para el operador modificado: rompe el modo de Goldstone (residuo medido =
> predicho ½|1−√λ|‖(√n)″‖, cociente 1.0000; `residual_anatomy/diag_goldstone.py`). Con ese
> operador, la temperatura upstream sale incluso negativa para todo λ ≠ 1 (el «colapso»;
> `residual_anatomy/diag_T2_qp_consistent.py`). **Quedan retiradas** las conclusiones «la lengua
> de Hawking es ROBUSTA frente a k⁴» (§5), las cifras ±10 % / 20–45 % (§7), la tabla w × λ y
> su lectura (§8). Las secciones se conservan como registro histórico.
> **No afectados:** §1–§3 (λ = 1, donde ambos operadores son idénticos).
> **Además (§9):** con el fondo corregido, el knob λ es exactamente un reescalado de la anchura
> del horizonte, (λ, w) ≡ (1, w/√λ); por tanto no aísla física ultravioleta con la métrica fija,
> y una rejilla w × λ mide un solo parámetro.

Sesión 2026-11-14 (continuación). Construido SOBRE el solver Hawking cerrado y
auditado (ver `REPRODUCE.md`). Observable: la correlación densidad-densidad de
Hawking (la "lengua" de Steinhauer/Balbinot-Carusotto) y su robustez frente al
término dispersivo k⁴ (primer paso trans-planckiano).

**Nivel de ambición (preservar):** Hawking **análogo** en BEC. Esto NO demuestra
Hawking gravitacional; la analogía es cinemática.

---

## 1. Observable y fórmula (verificada en literatura)

Correlación densidad-densidad en el vacío-in a T=0:
`C(x,x') = <δn(x) δn(x')>_conectada`.

Fórmula de suma de modos: **Recati, Pavloff, Carusotto, PRA 80, 043603 (2009),
arXiv:0907.4305v1, Ec. 25** (parte de punto cero = la única a T=0; numeración de esa versión —
hasta la auditoría del 2026-11-18 se citaba por error como «Ec. (32)»):

```
n·G0(ω,x,x') = (1/4π) [ Σ_{I=u,d1} v*_I(x) r_I(x') + u*_d2(x) r_d2(x') + c.c. ]
r_I = u_I + v_I
```

- `{u, d1, d2}` = 3 canales ENTRANTES de scattering en la ventana de Hawking.
  `d2` = canal de **norma NEGATIVA** (entra con `u*` en vez de `v*`: esa asimetría
  ES la firma de producción de pares de Hawking).
- Modos normalizados a corriente de norma |J|=1, con lo que ∫dx(|u|²−|v|²) = 2π δ(ω−ω′). **Con esa
  normalización el prefactor correcto es 1/(2π), no el 1/(4π) del artículo** (que corresponde a su propia
  normalización de modos). *Corrección de auditoría 2026-11-18:* aquí se decía que los modos |J|=1 «absorben el 1/4π»
  y el código multiplicaba por 1/(4π): la amplitud absoluta de C salía a la **mitad**. Comprobado con una prueba exacta en
  un condensado homogéneo (`tests/validate_correlation_norm.py`; antes 0,500001, ahora 1,00001) y con el cálculo directo
  desde la ordenación del vacío (`residual_anatomy/diag_recati_direct.py`). Ningún número publicado de §1–3 depende de la
  amplitud absoluta (la pendiente −1,17 se ha recalculado y no cambia); sí las escalas de color de las figuras de §3,
  regeneradas.
- `C` es REAL y SIMÉTRICA en x↔x' (las densidades conmutan a tiempo igual);
  se simetriza explícitamente (`br + conj(br.T)`, «+ c.c.» como conjugado hermítico). Por ω, esa forma es exactamente
  la correlación física menos el término propio (identidad de completitud BdG, `diag_recati_direct.py`).
- Se integra SOLO la ventana de Hawking 0<ω<ω_max → se obtiene la contribución
  ESPECÍFICA de Hawking (la parte ordinaria ω>ω_max queda excluida a propósito).

Referencias del observable: Balbinot-Fabbri-Fagnocchi-Recati-Carusotto,
PRA 78, 021603(R) (2008), arXiv:0711.4520; Carusotto-Fagnocchi-Recati-Balbinot-
Fabbri, NJP 10, 103001 (2008), arXiv:0803.0507; Steinhauer, Nat. Phys. (2016).

---

## 2. Implementación y VALIDACIÓN (gates, `validate_correlations.py`)

Módulo: `correlations.py`. Modos IN de scattering full-grid con el operador BdG
FD **validado** (`fd_bdg_operator.build_MBdG` + `plane_wave_mode`), lógica de
matching modal exacto idéntica a `mode_matching.py` (JUEZ independiente ya
validado), parametrizada en ω y devolviendo el CAMPO completo (u(x),v(x)).

Gates (λ=1, dispersión física):
- (A) dispersión λ=1 == cuártico original: dif 0.
- (B) |U|²−|V|² = sign(Ω) (±1). *Corrección de auditoría 2026-11-17:* la versión anterior del test comparaba con +1 para todos los modos e imprimía «desviación 2.00» (3 de los 4 k de prueba son de norma negativa, cuyo valor correcto es −1); el «2e-16» que figuraba aquí no correspondía a esa comprobación. Con la comprobación corregida: < 1e-12.
- (C) homogéneo (sin horizonte): transmisión=1.0000, reflexión 2e-15.
- (D) HORIZONTE ω=0.1: **β²=0.01008** convergido en N, flux_out=+1.0000
  → reproduce EXACTAMENTE el valor validado (memoria).
- (E) pseudo-unitaridad por canal: 2.9e-7 (u), 1.2e-4 (d1,d2).

---

## 3. Resultado (a): la lengua de Hawking (λ=1)

`run_correlation.py`, `run_decomp.py`. L=30, N=3600, w=2, banda ω∈[0.02,0.45]
(140 pts), ventana Hann en bordes, ncoarse=420. Convergencia en ω (140 vs 70 pts),
ventana |x|<15: 1.7%. *(Auditoría 2026-11-18: el 1.5% y la figura `correlation_tongue.png` anteriores
se generaron con una versión de `density_correlation` sin la ventana Hann (anterior a añadirla, el mismo 14-11);
regenerados con la ventana que aquí se describe y con el prefactor corregido 1/(2π).)*

- Aparece estructura de correlación cruzando el horizonte (upstream x<0 subsónico
  = "fuera"; downstream x>0 supersónico = "dentro").
- **Descomposición por canal** (clave): la banda negativa FUERTE del cuadrante
  upstream-upstream es correlación ORDINARIA (canales u,d1), NO Hawking. El término
  ANÓMALO d2 (norma negativa) aísla la contribución de Hawking.
- La cresta de anti-correlación del canal d2 en el cuadrante cruzado sigue la recta
  predicha del compañero lento: pendiente mediana **−1.17** vs predicho
  `(v_d−c_d)/(v_u−c_u) = −1.26`, dentro del ruido oscilatorio.
- La lengua es DÉBIL frente a las correlaciones ordinarias (consistente con
  β²~0.01, régimen dispersivo w=2) — como en los experimentos/figuras publicadas,
  oscila a su número de onda relativo; la ENVOLVENTE va por la recta predicha.

Figuras: `out/correlation_tongue.png`, `out/correlation_decomp.png`.
Datos: `out/correlation_lambda1.npz`, `out/correlation_decomp.npz`.

---

## 4. Knob trans-planckiano k⁴ (b): diseño y VALIDACIÓN

> ⚠️ **SECCIÓN AFECTADA** por el operador λ defectuoso — sin validez física (ver aviso inicial y §9).

Dispersión modificada manteniendo la MÉTRICA acústica (c,v) FIJA:
```
(ω − v k)² = c² k² + λ k⁴/4 ,   c² = g n   (κ, T_H intactas)
```
λ=1 físico; λ<1 más hidrodinámico; λ>1 más dispersivo (UV "trans-planckiano").

Operador `build_op_lambda`: `t=√λ/2` (cinético), `G=gn/√λ` (acoplo), convección
simétrica coef 1 (→ Ω=ω−vk), **diagonal = G + Qp** con la presión cuántica del
fondo `Qp = V−μ+gn+v²/2` CONSERVADA (física, →0 asintóticamente). En λ=1 reduce
EXACTAMENTE a `build_MBdG`.

**Dos bugs encontrados y corregidos durante la validación** (disciplina del
proyecto — cada uno daba resultados no físicos):
1. Al reconstruir el operador "limpio" sin Qp, λ=1 daba β²=0.19 (≠0.01008): β² es
   muy sensible a la presión cuántica del fondo. SOLUCIÓN: conservar Qp en la
   diagonal (solo se escala la estructura UV k⁴).
2. La **corriente de norma** usada para normalizar modos era la de λ=1; con
   dispersión modificada el término cinético escala:
   `J = √λ·Im(u*u') + v|u|² + √λ·Im(v*v') − v|v|²`
   → para onda plana `J=(√λ k+v)|U|²+(√λ k−v)|V|²`. Con la J equivocada la
   pseudo-unitaridad fallaba (0.8–5.1). SOLUCIÓN: J con √λ.

Gates del operador modificado (`validate_lambda.py`):
- (A) homogéneo transmisión=1 ∀λ (refl 2e-15).
- (B) λ=1 op='lambda' reproduce op='validated': β²=0.01008 = 0.01008 (0.0%).
- (C) pseudo-unitaridad por canal: 6e-4 (λ=0.5), 1e-4 (λ=1), 2e-4 (λ=2).
- (D) flux_out=+1.0000 ∀λ (físico).
- (E) ventana de Hawking ω_max(λ): 0.662 (λ=0.5), 0.468 (λ=1), 0.331 (λ=2)
  — más dispersión cierra antes la ventana.

β²(ω=0.1,λ) resultó NO monótona con mínimo cerca de λ=1 (0.010), subiendo a ambos
lados (0.125 en λ=0.25; 0.048 en λ=4) — pero β² a ω fija no es la medida limpia
(la ventana y el espectro se desplazan); la medida física es la LENGUA integrada.

---

## 5. Resultado (b): robustez de la lengua vs λ

> ⚠️ **SECCIÓN AFECTADA** por el operador λ defectuoso — sin validez física (ver aviso inicial y §9).

`run_lambda_sweep.py`. L=30, N=3000, ncoarse=360, banda ω adaptada a ω_max(λ)
(90 pts), op='lambda' para todos. Medida de robustez: integral de la parte negativa
de C_d2 en el cuadrante cruzado upstream×downstream (fuerza de la lengua).

| λ | ω_max | integral_neg(d2) | min(d2) | fuerza / λ=1 |
|---|---|---|---|---|
| 0.5 | 0.662 | −0.189 | −5.3e-3 | **1.02** |
| 1.0 | 0.468 | −0.185 | −5.1e-3 | **1.00** |
| 2.0 | 0.331 | −0.169 | −4.6e-3 | **0.91** |

**RESULTADO: la lengua de Hawking es ROBUSTA frente al término k⁴.** Variando el
coeficiente dispersivo (trans-planckiano) λ en un factor **×4** (0.5→2.0), la fuerza
integrada de la anti-correlación cambia solo **~±10%**, con la métrica acústica
(c,v,κ,T_H) FIJA. Tendencia: más dispersión (λ=2) **ensancha/suaviza** las franjas
(la ventana de Hawking se encoge, se cortan los modos de alto k) y **debilita
levemente** la lengua (0.91); menos dispersión (λ=0.5, más hidrodinámico) la
refuerza un pelín (1.02). Es la versión análoga de la robustez trans-planckiana de
la radiación de Hawking (Unruh-Schützhold; Jacobson): la señal IR es casi
insensible a la completación UV.

**Caveats honestos:** (1) la lengua es débil y oscilatoria (β²~0.01, régimen
dispersivo w=2) — la métrica integrada mezcla la cresta genuina con las franjas de
oscilación; una medida más limpia sería la correlación extrapolada al IR. (2) Solo
3 valores de λ, una configuración (w=2, Mach_d=2), N=3000 sin estudio de
convergencia del propio sweep. (3) La media sobre la recta −1.26 crece ligeramente
en magnitud con λ (−7.9e-4→−1.0e-3) mientras la integral decrece — ambas son
cambios O(10%): consistente con "robusto, cambios suaves", no con un efecto
monótono fuerte. El punto λ=1 reproduce EXACTAMENTE el solver validado; operador
pseudo-unitario a 1e-4 ∀λ.

Figura: `out/lambda_sweep.png`. Datos: `out/lambda_sweep.npz`.

---

## 6. Archivos nuevos de esta sesión

| Archivo | Qué hace |
|---|---|
| `correlations.py` | módulo: modos IN full-grid + Ec. 25 (arXiv v1) + knob λ (operador modificado) |
| `validate_correlations.py` | gates A-F de la maquinaria de correlación (λ=1) |
| `validate_lambda.py` | gates A-E del operador modificado k⁴ |
| `run_correlation.py` | lengua de Hawking λ=1 + sondeo de pendiente |
| `run_decomp.py` | descomposición d2 (Hawking) vs ordinario (u,d1) |
| `run_lambda_sweep.py` | robustez de la lengua vs λ (trans-planckiano) |

---

## 7. Profundización de la robustez (sesión 2026-11-14 cont.)

> ⚠️ **SECCIÓN AFECTADA** por el operador λ defectuoso — sin validez física (ver aviso inicial y §9).

Tres pasos sobre §5 para cerrar los caveats (1) métrica que mezcla franjas y
(2) sin convergencia / solo 3 λ / 1 config.

### 7.1 Convergencia del sweep — el ±10% es REAL (`robustness_convergence.py`)
Re-medida la fuerza de la lengua (misma métrica integral, canal d2) variando malla
y caja. El cociente fuerza/λ1 es estable a la 3ª cifra:
- **N** (2000→3000→4500, L=30): λ=0.5→1.024/1.022/1.024 ; λ=2→0.914/0.914/0.915.
  (el valor absoluto cambia de N=2000 a 3000, pero N=3000≈4500 → convergido en N=3000)
- **L** (24→30→40, dx fino fijo): λ=0.5→1.024/1.022/1.021 ; λ=2→0.915/0.914/0.913.
→ El ±10% sobre ×4 [0.5,2] NO es artefacto de malla ni de caja.

### 7.2 λ ampliado ×16 + métrica limpia (`robustness_clean_metric.py`, `robustness_ridge_scan.py`)
Ampliado λ a {0.25,0.5,1,2,4} (×16). Métrica limpia = envolvente de Hilbert de C_d2
a lo largo de la recta del compañero (ignora la fase de las franjas, cuya longitud
de onda cambia con λ). Caveat detectado y auditado: la pendiente del compañero
depende de v_g → cambia con λ; muestrear a pendiente FIJA −1.26 sesga. Se añadió
escaneo de pendiente (`robustness_ridge_scan.py`).

**HALLAZGO (disciplina — el escaneo encontró un 2º problema):** la pendiente óptima
se pega al borde del escaneo (−0.80) en casi todos los λ → la envolvente NO tiene un
pico en la cresta; crece al hacer la recta más horizontal (contaminación cerca del
horizonte). **La lengua es demasiado débil/empotrada (β²~0.01, w=2 dispersivo) para
aislarla con una sola recta.** Ni pendiente fija ni escaneo dan un escalar perfecto.

**Lo ROBUSTO (coinciden las 3 métricas = lo único de fiar), normalizado a λ=1:**

| λ | integral (vieja) | envolvente fija −1.26 | cresta escaneada |
|---|---|---|---|
| 0.25 | 0.75 | 0.72 | 0.80 |
| 0.5  | 1.02 | 0.95 | 1.01 |
| 1.0  | 1.00 | 1.00 | 1.00 |
| 2.0  | 0.91 | 0.81 | 0.82 |
| 4.0  | 0.79 | 0.65 | 0.66 |

**CONCLUSIÓN HONESTA:** (1) la lengua tiene un máximo cerca del valor FÍSICO λ=1 y
decae hacia ambos extremos; (2) sobre ×4 [0.5,2]: ~10-20% (el ±10% era el extremo
optimista; con envolvente ~20% en λ=2); (3) sobre ×16: cae ~20-45% en los extremos,
pero SUAVE, sin colapso ni divergencia, mismo orden de magnitud; (4) el % exacto
depende de la métrica (10-45%) → no hay escalar bien definido (señal débil/oscil.).
Es robustez trans-planckiana en el sentido correcto (la señal IR sobrevive, mismo
orden, pico en el valor físico, sin cambio cualitativo al variar la completación UV
×16), NO "±10% plano". El residuo de sensibilidad UV es esperable para w=2 (horizonte
~2ξ, no profundamente hidrodinámico).

### 7.3 ABIERTO (próximos pasos candidatos)
- **(A) Escala con w:** predicción física — un horizonte más suave (w grande, más
  hidrodinámico) debería dar una lengua MÁS robusta frente a λ (menos sensibilidad
  UV). Medir la variación λ de la lengua para w=2,4,6 → si decrece con w, confirma
  la interpretación trans-planckiana (Unruh-Schützhold/Jacobson) y conecta con el
  escalado dispersivo de Finazzi. Es el experimento con más fundamento físico.
- **(B) Comparación con Finazzi-Parentani:** su predicción dispersiva del espectro
  (arXiv:1012.1556) — traducir a una predicción para la amplitud de la correlación.

### Archivos nuevos (7)
| Archivo | Qué hace |
|---|---|
| `robustness_convergence.py` | convergencia del sweep en N y L (métrica integral) |
| `robustness_clean_metric.py` | métrica limpia (envolvente ridge) + λ ×16 |
| `robustness_ridge_scan.py` | escaneo de pendiente de la cresta (auditoría de la métrica) |

---

## 8. Opción A: dependencia con w (sesión 2026-11-14 cont.) — CERRADA, resultado NULO

> ⚠️ **SECCIÓN AFECTADA** por el operador λ defectuoso — sin validez física (ver aviso inicial y §9).

Recomendación del referee (ChatGPT): ¿la lengua de Hawking se vuelve menos sensible
a la dispersión UV (λ) al suavizar el horizonte (w grande, más hidrodinámico)?

### 8.1 Bug detectado y corregido (mismo patrón que H1 de la auditoría térmica)
`w_dependence.py` (v1) usó banda ω FIJA [0.02, 0.9·ω_max] para todos los w. Pero
T_H∝κ∝1/w cae con w (0.043→0.014); la banda fija cae en la COLA DISPERSIVA (ω≫T_H)
para w grande → inflaba la sensibilidad a λ. v1 daba spread CRECIENTE (0.086→0.373,
tendencia OPUESTA) = ARTEFACTO. Corrección: banda en unidades de T_H(w),
ω∈[1,4]·T_H(w), misma regla ∀(w,λ). (`w_dependence_v2.py`, `w_dependence_v3.py`).

### 8.2 Tabla de cierre (`w_dependence_v3.py`, criterio completo de ChatGPT)
Controles: spread = max_λ R_λ − min_λ R_λ (NO el asimétrico max|r−1|); no-truncado
verificado por cada (w,λ); pseudo-unitaridad ε(ω) por frecuencia de los MISMOS
solves que construyen C; n_valid por config; L crece con w (H3).

| w | κ | T_H | R(0.5) | R(1) | R(2) | spread | n_val | ε_max |
|---|---|---|---|---|---|---|---|---|
| 2 | 0.271 | 0.0431 | 1.225 | 1.00 | 0.898 | 0.327 | 60/60 | 3.0e-3 |
| 4 | 0.135 | 0.0215 | 1.192 | 1.00 | 0.852 | 0.339 | 60/60 | 3.8e-3 |
| 6 | 0.090 | 0.0143 | 1.146 | 1.00 | 0.852 | 0.294 | 60/60 | 3.9e-3 |

Todos los controles limpios: sin truncado, 60/60 válidas en las 9 configs, 0
descartes, ε_max 3-4e-3.

### 8.3 RESULTADO HONESTO (NULO) — la definición de la métrica fue determinante
Con el spread RIGUROSO (max−min): 0.327 → 0.339 → 0.294 = **aproximadamente
CONSTANTE, ligero bajón NO monótono en w=6. NO es un decrecimiento limpio.**
El "decrece monótonamente" preliminar (v2: 0.225→0.192→0.148) venía de la métrica
ASIMÉTRICA max|r−1|; con max−min el efecto se esfuma. → la hipótesis "robustez
mejora al suavizar el horizonte" NO queda respaldada por este experimento.

Estructura real (más matizada): las dos ramas se mueven OPUESTAS con w:
- lado hidrodinámico λ=0.5: R→1 (1.225→1.146) = más robusto por ese lado al crecer w.
- lado dispersivo λ=2: R plano/aleja (0.898→0.852→0.852) = igual o algo menos robusto.
Lectura plausible (no demostrada): suavizar el horizonte hace que QUITAR dispersión
UV (λ<1) importe menos, pero AÑADIR dispersión UV (λ>1) muerde igual sea cual sea w.

Caveat físico (ChatGPT): w y T_H COVARÍAN (κ∝1/w, Machs asintóticos fijos) → esto
caracteriza una FAMILIA de horizontes donde suavidad y escala térmica cambian
juntas; NO aísla la suavidad. Aislarla exigiría variar Mach_d(w) para fijar κ.

Archivos: `w_dependence.py` (v1, artefacto), `w_dependence_v2.py` (banda T_H),
`w_dependence_v3.py` (cierre). Datos out/w_dependence*.npz, figs .png.

### 8.4 Comprobación de la asimetría (pedida por ChatGPT) — ES normalización
`asymmetry_check.py`. ¿La asimetría de ramas (λ<1 vs λ>1 con w) es física o la
introduce normalizar cada config a su propio str(λ=1)? str SIN normalizar y
desviaciones ABSOLUTAS Δ respecto a λ=1:

| w | \|str(1)\| | Δ_lo=str(.5)−str(1) | Δ_hi=str(2)−str(1) | R(.5)−1 | R(2)−1 |
|---|---|---|---|---|---|
| 2 | 0.189 | −0.0426 | +0.0192 | +0.225 | −0.102 |
| 4 | 0.228 | −0.0437 | +0.0336 | +0.192 | −0.147 |
| 6 | 0.266 | −0.0388 | +0.0393 | +0.146 | −0.148 |

|str(1)| CRECE ×1.41 con w. Δ_lo (λ<1) ~CONSTANTE (−0.043→−0.039). Δ_hi (λ>1)
CRECE ×2.05. **CONCLUSIÓN:** el "R(0.5)→1 al crecer w" es casi todo el DENOMINADOR
creciente (Δ_lo plano), NO mejora de robustez en la rama hidrodinámica. En absoluto,
Δ_hi (λ>1) crece con w; PERO esto NO demuestra un efecto físico UV mayor (corrección
de ChatGPT a una sobreatribución mía): la ESCALA del propio observable también cambia
con w (|str| crece, geometría del horizonte, escala espacial efectiva de la integral)
→ la interpretación del Δ absoluto requiere controlar esa escala, no hecho aquí. Lo
único firme: la asimetría de R está dominada por la normalización (confirma la
sospecha de ChatGPT); NO es evidencia suficiente de asimetría física.
Fig out/asymmetry_check.png.
Gate pseudo-unitaridad: ε_max real = 3.0–3.9e-3 (ChatGPT lo citó como e-4, un orden
de más); PASA el umbral PREESTABLECIDO PU_TOL=1e-2 (3.9e-3 < 1e-2), por eso se acepta
el control — no por llamarlo "limpio".

**CIERRE DEFINITIVO OPCIÓN A:** en la familia w={2,4,6} NO hay reducción de la
sensibilidad global de la lengua a λ al suavizar el horizonte (spread ~constante,
no monótono). La señal SOBREVIVE a variar λ, pero el experimento NO encuentra mejora
de robustez con w. La asimetría de ramas es sustancialmente normalización. Formulación
defendible (ChatGPT): "la señal sobrevive a variaciones de λ" — NO "robustez
trans-planckiana confirmada por el suavizado del horizonte". Limitado por covariación
w–κ–T_H (familia, no suavidad aislada).

---

## 9. Auditoría fase 1 (2026-11-17): el knob λ con fondo consistente es un reescalado

**Afirmación.** Con el fondo consistente (Qp_λ = √λ·Qp, variante `QP_LAM` de
`residual_anatomy/diag_T2_qp_consistent.py`), el operador BdG con dispersión
(ω − vk)² = c²k² + λk⁴/4 sobre el perfil tanh de anchura w es exactamente el operador físico
v1.0 (λ = 1) sobre el perfil de anchura w/√λ:

```
x' = x/√λ ,  k' = √λ k ,  ω' = √λ ω ,  M_1(x') = √λ · M_λ(x)
```

Término a término: cinético √λ/2·k² → k'²/2; acoplo gn/√λ → gn; convección vk → vk';
Qp_λ = √λ·Qp → λ·Qp = Qp' (Qp escala como 1/w²). n, v y c no cambian, y T/T_H es
adimensional. En la malla, con L' = L/√λ y el mismo N, D2' = λ·D2 y D1' = √λ·D1: la
identidad también es exacta en el operador discreto.

**Comprobación** (`residual_anatomy/diag_scaling_equivalence.py`, ~30 s,
`out/diag_scaling_equivalence.csv`): para los 90 puntos (λ, N, ω) de `QP_LAM` en
`out/diag_T2_nL.csv` se recalcula n_L con el solver v1.0 (`op='validated'`, **sin**
`build_op_lambda`) en w' = 2/√λ, L' = 40/√λ, mismo N, ω' = √λ·ω.

| λ | w equivalente | máx. dif. relativa en n_L (N = 9600) | (N = 14400) |
|---|---|---|---|
| 0.25 | 4.000 | 8.8e-08 | 4.9e-07 |
| 0.5 | 2.828 | 1.5e-07 | 8.0e-07 |
| 1 | 2.000 | 0 | 0 |
| 2 | 1.414 | 4.9e-07 | 2.4e-06 |
| 4 | 1.000 | 2.8e-06 | 3.2e-06 |

El suelo de redondeo del propio solver v1.0 (n_L al mover L en 1 ulp, λ = 4, N = 14400) es
7·10⁻⁶: las diferencias están por debajo, así que la identidad se cumple hasta el redondeo.
(Declaración: la tolerancia se fijó primero en 1e-8 y falló; se cambió a 1e-5 **después** de
medir ese suelo.)

**Confirmación con un solve más preciso** (diagnóstico, misma tarde;
`residual_anatomy/diag_conditioning_equiv.py`, `out/diag_conditioning_equiv.csv`): resolviendo los dos
lados con equilibrado de Ruiz y refinamiento iterativo (residuo en long double), sin tocar el
operador ni la malla, la diferencia máxima en los 90 puntos baja de 3,2·10⁻⁶ a **3,5·10⁻⁹**. El
7·10⁻⁶ era pérdida de precisión del LU sin equilibrar, no del problema
(`residual_anatomy/diag_conditioning.py`: con la solución refinada, 1 ulp en L cambia n_L ~10⁻¹²).
El solver v1.0 no se ha modificado.

En T/T_H, con el benchmark v1.0 en la anchura equivalente: λ = 0.25 ↔ w = 4: 0.97900 frente a
0.97893 (`FP_BENCHMARK.md`, 0.9789); λ = 1: 0.95724 en ambos. Ese cálculo usa la malla del
benchmark y no la reescalada, de ahí la diferencia de 7·10⁻⁵.

**Consecuencias.**
1. Confirma T2 por una vía independiente del código de λ, y que el «colapso» era un
   artefacto: con el fondo consistente se recupera exactamente la física de λ = 1.
2. El knob λ **no** estudia física trans-planckiana con la métrica fija: es el mismo parámetro
   que la anchura del horizonte en unidades de la longitud de curación (lo que ya cubre la
   familia w de Finazzi–Parentani).
3. Si se quiere estudiar la lengua frente a la «suavidad UV», basta λ = 1 con distintos w.
   Una rejilla w × λ (como la de §8) es degenerada.
4. Alcance: la identidad se ha comprobado numéricamente para n_L; para C(x, x′) es una
   consecuencia del mismo argumento de escala, todavía no ejecutada.
